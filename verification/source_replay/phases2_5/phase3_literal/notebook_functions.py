def per_lineage_DE(drug_name, lineage, min_per_group=5):
    """
    Run M vs F differential expression for ONE drug × ONE lineage.

    Args:
      drug_name: drug name (e.g., 'vorinostat')
      lineage: CMAP cell_lineage (e.g., 'lung')
      min_per_group: minimum cell lines needed in EACH of M and F (default 5)

    Returns:
      DataFrame with gene_id, gene_symbol, ensembl_id, difference, p_value, p_adj, score
      OR None if not enough samples in this lineage
    """
    # Step 1: Pick this drug × this lineage's experiments
    sigs = sig_known[
        (sig_known['cmap_name'] == drug_name) &
        (sig_known['cell_lineage'] == lineage)
    ]
    m_sigs = sigs[sigs['donor_sex'] == 'M']
    f_sigs = sigs[sigs['donor_sex'] == 'F']

    # Step 2: Check if enough samples
    n_m_cells = m_sigs['cell_iname'].nunique()
    n_f_cells = f_sigs['cell_iname'].nunique()
    if n_m_cells < min_per_group or n_f_cells < min_per_group:
        return None

    # Step 3: Pull only these experiments from gctx (small extract)
    all_ids = list(sigs['sig_id'])
    valid_ids = [s for s in all_ids if s in gctx_sigs]
    if len(valid_ids) < 20:
        return None

    drug_data = parse(gctx_file, cid=valid_ids)
    expr = drug_data.data_df  # 12,328 genes × N experiments

    # Step 4: Separate M vs F columns
    sex_map = sigs.set_index('sig_id')['donor_sex'].to_dict()
    m_cols = [s for s in expr.columns if sex_map.get(s) == 'M']
    f_cols = [s for s in expr.columns if sex_map.get(s) == 'F']

    if len(m_cols) < 10 or len(f_cols) < 10:
        return None  # too few experiments

    m_expr = expr[m_cols].values
    f_expr = expr[f_cols].values

    # Step 5: Mann-Whitney U test per gene
    p_vals, diffs = [], []
    for i in range(expr.shape[0]):
        try:
            _, p = mannwhitneyu(m_expr[i], f_expr[i], alternative='two-sided')
        except:
            p = 1.0
        p_vals.append(p)
        diffs.append(m_expr[i].mean() - f_expr[i].mean())

    # Step 6: Build results table
    de = pd.DataFrame({
        'gene_id': expr.index.astype(str),
        'difference': diffs,
        'p_value': p_vals
    })
    de['p_adj'] = multipletests(de['p_value'], method='fdr_bh')[1]

    # GSEA needs a ranked score: signed -log10(p)
    de['score'] = np.sign(de['difference']) * -np.log10(de['p_value'] + 1e-300)

    # Add Ensembl IDs + gene symbols
    de['ensembl_id'] = de['gene_id'].map(gid_to_ensg)
    de = de.merge(
        genes_info[['gene_id', 'gene_symbol']].astype({'gene_id': str}),
        on='gene_id', how='left'
    )

    return de.sort_values('score', ascending=False)

def run_gsea(de_results, gene_sets_dict, min_genes=15):
    """
    Run GSEA preranked on DE results against multiple gene sets.

    Args:
      de_results: DataFrame from per_lineage_DE() with 'ensembl_id' and 'score'
      gene_sets_dict: dict of {set_name: set of Ensembl IDs}
      min_genes: minimum genes a set needs (default 15, gseapy requires this)

    Returns:
      DataFrame with NES, FDR, lead genes per gene set
    """
    # Step 1: Build ranked list (sort by score, drop NaN)
    ranked = de_results.dropna(subset=['ensembl_id', 'score']).copy()
    ranked = ranked.drop_duplicates(subset='ensembl_id', keep='first')
    ranked = ranked.set_index('ensembl_id')['score'].sort_values(ascending=False)

    # Step 2: Convert sets to lists (gseapy expects lists)
    gene_sets_for_gsea = {
        name: list(genes & set(ranked.index))
        for name, genes in gene_sets_dict.items()
    }
    gene_sets_for_gsea = {
        name: genes for name, genes in gene_sets_for_gsea.items()
        if len(genes) >= min_genes
    }

    if not gene_sets_for_gsea:
        return pd.DataFrame()

    # Step 3: Run GSEA prerank
    pre_res = gseapy.prerank(
        rnk=ranked.reset_index(),
        gene_sets=gene_sets_for_gsea,
        outdir=None,           # don't save to disk
        permutation_num=1000,  # for FDR
        seed=42,
        verbose=False,
        no_plot=True,
        min_size=min_genes,
        max_size=10000
    )

    return pre_res.res2d

def analyze_drug_full(drug_name, verbose=True):
    """
    Full per-lineage analysis for one drug.
    Returns dict with per-lineage DE/GSEA results and final A/B/C category.
    """
    # Step 1: Find all eligible lineages for this drug
    drug_sigs = sig_known[sig_known['cmap_name'] == drug_name]
    drug_lineages = set(drug_sigs['cell_lineage'].dropna().unique()) - {'unknown'}

    if verbose:
        print(f"\n=== {drug_name} ===")
        print(f"Lineages drug was tested in: {len(drug_lineages)}")

    # Step 2: For each lineage with enough samples, run DE + GSEA
    lineage_results = {}
    for lineage in drug_lineages:
        de = per_lineage_DE(drug_name, lineage)
        if de is None:
            continue

        if verbose:
            n_strong = ((de['p_adj']<0.05) & (de['difference'].abs()>0.5)).sum()
            print(f"  Lineage '{lineage}': {n_strong} strong DE genes")

        # Run GSEA
        try:
            gsea = run_gsea(de, sb_library_sets)
        except Exception as e:
            if verbose:
                print(f"    GSEA error: {e}")
            continue

        lineage_results[lineage] = {
            'de': de,
            'gsea': gsea,
            'sig_tissues': gsea[gsea['FDR q-val'] < 0.05]['Term'].tolist()
        }

    # Step 3: Decide A/B/C based on COMBINED signal
    category = 'C'
    cat_a_evidence = []
    cat_b_evidence = []

    for lineage, res in lineage_results.items():
        sig_tissues = res['sig_tissues']
        if not sig_tissues:
            continue

        # Drug's home tissues for this lineage
        home_tissues = set(cmap_to_gtex.get(lineage, []))

        for t in sig_tissues:
            nes = res['gsea'][res['gsea']['Term']==t]['NES'].iloc[0]
            fdr = res['gsea'][res['gsea']['Term']==t]['FDR q-val'].iloc[0]
            if t in home_tissues:
                cat_a_evidence.append((lineage, t, nes, fdr))
            else:
                cat_b_evidence.append((lineage, t, nes, fdr))

    if cat_a_evidence:
        category = 'A'
    elif cat_b_evidence:
        category = 'B'

    if verbose:
        print(f"\nCategory: {category}")
        if cat_a_evidence:
            print(f"  Cat A evidence (lineage, tissue, NES, FDR):")
            for ev in cat_a_evidence[:5]:
                print(f"    {ev}")
        elif cat_b_evidence:
            print(f"  Cat B evidence (lineage, tissue, NES, FDR):")
            for ev in cat_b_evidence[:5]:
                print(f"    {ev}")

    return {
        'drug': drug_name,
        'category': category,
        'n_lineages_tested': len(lineage_results),
        'cat_a_evidence': cat_a_evidence,
        'cat_b_evidence': cat_b_evidence,
        'lineage_results': lineage_results
    }
