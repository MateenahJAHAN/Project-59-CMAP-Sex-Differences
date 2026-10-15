"""
Phase 6 - step 1: build the A / B / C / D table.

Marouen's instruction for this phase: run the pipeline with NO sample-size
filter. The only filter on the data is |z| > 2.

Output: results/per_signature_v4_ABCD.csv  (39,996 rows = drug x lineage x tissue)

Verbatim from Phase6_NoSampleFilter/notebook/phase6_no_sample_filter.ipynb
(cells 2, 5 and 7). Paths assume Google Colab with Drive mounted.
"""

# ---------------------------------------------------------------
# load reference objects  (sets, NOT lists - a list here silently
# breaks the set intersection inside run_gsea)
# ---------------------------------------------------------------
import json, pandas as pd

BASE = '/content/drive/MyDrive/Project59'
REF  = f'{BASE}/data/reference'

# ---- reference objects (sets, lists nahi — yeh bug tha pehle) ----
cmap_to_gtex    = {k: set(v) for k,v in json.load(open(f'{REF}/cmap_to_gtex.json')).items()}
sb_library_sets = {k: set(v) for k,v in json.load(open(f'{REF}/sb_library_sets.json')).items()}
gid_to_ensg     = json.load(open(f'{REF}/gid_to_ensg.json'))
top50_names     = json.load(open(f'{REF}/top50_names.json'))

# ---- dataframes ----
genes_info = pd.read_csv(f'{BASE}/data/cmap/geneinfo.csv')
sig_known  = pd.read_csv(f'{BASE}/results/sig_known/sig_known.csv', low_memory=False)

# ---- verify ----
print(f"cmap_to_gtex    : {len(cmap_to_gtex)} lineages   (expected 20)")
print(f"sb_library_sets : {len(sb_library_sets)} tissues  (expected 44)")
print(f"  type check    : {type(list(sb_library_sets.values())[0]).__name__}  (MUST be 'set')")
print(f"gid_to_ensg     : {len(gid_to_ensg):,} genes")
print(f"top50_names     : {len(top50_names)} drugs")
print(f"genes_info      : {len(genes_info):,} rows")
print(f"sig_known       : {len(sig_known):,} rows")
print(f"  columns       : {list(sig_known.columns)}")

# ---------------------------------------------------------------
# the two core functions
# ---------------------------------------------------------------
from cmapPy.pandasGEXpress.parse import parse
import pandas as pd, numpy as np, gseapy

gctx_file = '/content/cmap/level5_beta_trt_cp_n720216x12328.gctx'

# gctx ke andar kaunse signatures hain
col_meta  = parse(gctx_file, col_meta_only=True)
gctx_sigs = set(col_meta.index)
print(f"gctx signatures: {len(gctx_sigs):,}")


def per_lineage_DE_v3(drug_name, lineage, min_sigs=3):

    # check: gctx have minimum signatures .

    sigs = sig_known[(sig_known['cmap_name']   == drug_name) &
                     (sig_known['cell_lineage'] == lineage)]

    valid_ids = [s for s in sigs['sig_id'] if s in gctx_sigs]
    if len(valid_ids) < min_sigs:
        return None

    expr   = parse(gctx_file, cid=valid_ids).data_df
    mean_z = expr.mean(axis=1)# avg of every gene
    strong = mean_z[mean_z.abs() > 2] # this is a filter

    if len(strong) == 0:# if in a signature there is no any gene with z+_ 2 return none
        return None

    de = pd.DataFrame({'gene_id': strong.index.astype(str), # create a table and make datatype same so that we can use after
                       'z_score': strong.values, # we have a dataframe with gene_id index , another col is z_score and score  both z score and score are same
                       'score'  : strong.values})
    de['ensembl_id'] = de['gene_id'].map(gid_to_ensg)# create another col ensemble id to match cmap and gtex gene kind of translator
    de = de.merge(genes_info[['gene_id','gene_symbol']].astype({'gene_id':str}),
                  on='gene_id', how='left')# primary value merge
    return de.sort_values('score', ascending=False)


def run_gsea(de_results, gene_sets_dict, min_genes=15):
    ranked = de_results.dropna(subset=['ensembl_id']).copy()
    ranked = ranked.groupby('ensembl_id')['score'].mean().sort_values(ascending=False)

    gene_sets = {name: list(genes & set(ranked.index))
                 for name, genes in gene_sets_dict.items()}
    gene_sets = {n: g for n, g in gene_sets.items() if len(g) >= min_genes}
    if not gene_sets:
        return pd.DataFrame()

    pre = gseapy.prerank(rnk=ranked, gene_sets=gene_sets,
                         permutation_num=1000, min_size=min_genes,
                         max_size=5000, seed=42, no_plot=True,
                         outdir=None, threads=4, verbose=False)
    return pre.res2d


print(" per_lineage_DE_v3 + run_gsea ready")

# ---------------------------------------------------------------
# the batch loop  (~24 min on Colab)
# ---------------------------------------------------------------
import pandas as pd, time

all_tissues = sorted(sb_library_sets.keys())

print(f"Drugs: {len(top50_names)} | Tissues: {len(all_tissues)}")
print("Sample-size filter REMOVED (Marouen's instruction)")
print("Categories:  A = home significant | B = cross significant")
print("             C = tested, NOT significant | D = could NOT be tested")
print("Time: ~30-50 min\n")

rows = []
t0 = time.time()

for i, drug in enumerate(top50_names, 1):
    lineages = sorted(set(sig_known[sig_known['cmap_name'] == drug]['cell_lineage']
                          .dropna().unique()) - {'unknown'})

    for lin in lineages:
        home = set(cmap_to_gtex.get(lin, []))

        s = sig_known[(sig_known['cmap_name'] == drug) &
                      (sig_known['cell_lineage'] == lin)]
        n_m   = s[s['donor_sex'] == 'M']['cell_iname'].nunique()
        n_f   = s[s['donor_sex'] == 'F']['cell_iname'].nunique()
        n_sig = s['sig_id'].nunique()

        base = dict(drug=drug, lineage=lin, n_M=n_m, n_F=n_f, n_signatures=n_sig)

        # D = could NOT be tested
        def all_untested(reason):
            for t in all_tissues:
                rows.append({**base, 'tissue': t, 'is_home_tissue': t in home,
                             'NES': None, 'FDR': None,
                             'category': 'D', 'skip_reason': reason})

        if not home:
            all_untested('no_GTEx_home_tissue')
            continue

        de = per_lineage_DE_v3(drug, lin)
        if de is None:
            all_untested('too_few_signatures_or_DE_genes')
            continue

        try:
            g = run_gsea(de, sb_library_sets)
        except Exception as e:
            print(f"  !! GSEA error: {drug} x {lin} -> {type(e).__name__}: {e}")
            all_untested('gsea_error')
            continue

        if g is None or len(g) == 0 or 'FDR q-val' not in g.columns:
            print(f"  !! GSEA empty: {drug} x {lin}")
            all_untested('gsea_empty')
            continue

        gres = g.drop_duplicates(subset='Term').set_index('Term')

        for t in all_tissues:
            # D = this tissue could not be tested (too few matching genes)
            if t not in gres.index:
                rows.append({**base, 'tissue': t, 'is_home_tissue': t in home,
                             'NES': None, 'FDR': None,
                             'category': 'D', 'skip_reason': 'gene_set_too_small'})
                continue

            fdr = float(gres.loc[t, 'FDR q-val'])
            nes = float(gres.loc[t, 'NES'])

            # A / B = significant   |   C = tested but not significant
            cat = ('A' if t in home else 'B') if fdr < 0.05 else 'C'

            rows.append({**base, 'tissue': t, 'is_home_tissue': t in home,
                         'NES': nes, 'FDR': fdr,
                         'category': cat, 'skip_reason': None})

    print(f"[{i:2}/50] {drug:24} · {len(lineages)} lin · {time.time()-t0:.0f}s")

sig_table = pd.DataFrame(rows)[
    ['drug', 'lineage', 'tissue', 'is_home_tissue', 'n_M', 'n_F',
     'n_signatures', 'NES', 'FDR', 'category', 'skip_reason']]

out = f'{BASE}/results/per_signature_v4_ABCD.csv'
sig_table.to_csv(out, index=False)

print(f"\n{'='*58}")
print(f"Rows: {len(sig_table):,}   Time: {(time.time()-t0)/60:.1f} min")
print('='*58)
print("A = home significant   B = cross significant")
print("C = tested, NOT significant   D = could NOT be tested\n")
print(sig_table['category'].value_counts().reindex(['A','B','C','D']).to_string())
print("\nD ke reasons (why not tested):")
print(sig_table[sig_table['category'] == 'D']['skip_reason']
      .value_counts().to_string())
print(f"\nDrugs with hits (A or B): "
      f"{sig_table[sig_table['category'].isin(['A','B'])]['drug'].nunique()} / 50")
print(f"\nSaved: {out}")
