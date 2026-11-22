def analyze_drug_lineage_v4(drug_name, lineage, verbose=False):
    """
    Phase 4: Per (drug, lineage) GSEA against split male/female libraries.
    - Uses z-scores DIRECTLY (no Mann-Whitney)
    - Mean z-score per gene across all signatures of (drug, lineage)
    - GSEA against 88 split sets
    - Direction comes from library name (Male_X or Female_X)
    """
    sigs = sig_known[
        (sig_known['cmap_name'] == drug_name) &
        (sig_known['cell_lineage'] == lineage)
    ]
    if len(sigs) < 10:
        return None

    valid_ids = [s for s in sigs['sig_id'] if s in gctx_sigs]
    if len(valid_ids) < 10:
        return None

    # Pull z-scores
    data = parse(gctx_file, cid=valid_ids)
    expr = data.data_df  # 12,328 genes × N sigs

    # Mean z-score per gene across signatures (collapse sigs to gene-level summary)
    mean_z = expr.mean(axis=1)
    mean_z.index = mean_z.index.astype(str)

    # Entrez → Ensembl translation
    ensg_scores = {}
    for gid, z in mean_z.items():
        ensg = gid_to_ensg.get(gid)
        if ensg and not pd.isna(ensg):
            ensg_scores[ensg] = z
    if not ensg_scores:
        return None

    # Sort by z-score (descending) — this IS the ranked list
    ranked = pd.Series(ensg_scores).sort_values(ascending=False)

    # Optional info: how many genes pass |z| > 2
    n_strong = (ranked.abs() > 2.0).sum()
    if verbose:
        print(f"  {drug_name} × {lineage}: {len(valid_ids)} sigs, {n_strong} genes |z|>2")

    # Run GSEA against 88 split libraries
    try:
        result = gseapy.prerank(
            rnk=ranked.reset_index(),
            gene_sets=sb_libs_split,
            permutation_num=200,   # 200 for speed; can bump to 1000 later
            min_size=15,
            max_size=10000,
            no_plot=True,
            verbose=False,
            seed=42
        )
        gsea_df = result.res2d
    except Exception as e:
        if verbose:
            print(f"  GSEA error: {e}")
        return None

    return {
        'drug': drug_name,
        'lineage': lineage,
        'n_sigs': len(valid_ids),
        'n_M': sigs[sigs['donor_sex']=='M']['cell_iname'].nunique(),
        'n_F': sigs[sigs['donor_sex']=='F']['cell_iname'].nunique(),
        'n_strong_z': int(n_strong),
        'gsea': gsea_df
    }

def classify_drug_lineage(result, cmap_to_gtex):
    """
    Convert GSEA result to A/B/C classification + sex direction.
    Each significant library (FDR < 0.05) → one row of classification.
    """
    if result is None:
        return []

    drug = result['drug']
    lineage = result['lineage']
    home_tissues = set(cmap_to_gtex.get(lineage, []))

    gsea = result['gsea']
    sig_rows = gsea[gsea['FDR q-val'] < 0.05]

    classifications = []
    for _, row in sig_rows.iterrows():
        term = row['Term']  # e.g., 'Female_Lung' or 'Male_Brain_Cortex'

        # Parse direction (sex) and tissue from library name
        if term.startswith('Female_'):
            sex_direction = 'Female'
            tissue = term.replace('Female_', '')
        elif term.startswith('Male_'):
            sex_direction = 'Male'
            tissue = term.replace('Male_', '')
        else:
            continue

        # Classify A/B/C
        if tissue in home_tissues:
            category = 'A'  # same-tissue match
        else:
            category = 'B'  # cross-tissue

        classifications.append({
            'drug': drug,
            'lineage': lineage,
            'tissue': tissue,
            'sex_direction': sex_direction,
            'category': category,
            'NES': row['NES'],
            'FDR': row['FDR q-val'],
            'n_sigs': result['n_sigs']
        })

    # If no significant libraries → Cat C
    if not classifications:
        classifications.append({
            'drug': drug, 'lineage': lineage, 'tissue': None,
            'sex_direction': None, 'category': 'C',
            'NES': None, 'FDR': None, 'n_sigs': result['n_sigs']
        })

    return classifications
