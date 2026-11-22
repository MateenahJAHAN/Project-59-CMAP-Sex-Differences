"""Exploratory (not pre-registered): nearest other line by deep-cell profile correlation (all conditions pooled,
no condition labels), to help interpret line-level disagreements. Same gene set as the identity check (plan 7d)."""
import sys, os, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(__file__))
from common import NGROUP, NBINS, Y_GENES
ROOT = sys.argv[1]
acc = None; cnt = None
for f in sorted(glob.glob(os.path.join(ROOT, "data/F1/emerald_proc/acc_w*.npz"))):
    z = np.load(f); acc = z["acc"].astype(float) if acc is None else acc + z["acc"]; cnt = z["cnt"] if cnt is None else cnt + z["cnt"]
cl = pd.read_parquet(os.path.join(ROOT, "data/F1/meta/cell_line_metadata.parquet")).drop_duplicates("Cell_ID_Cellosaur")
lines = sorted(cl.Cell_ID_Cellosaur); name = cl.set_index("Cell_ID_Cellosaur").cell_name
g = pd.read_parquet(os.path.join(ROOT, "data/F1/meta/gene_metadata.parquet"))
P = np.vstack([acc[[i * NGROUP + k for k in range(0, 1 + NBINS)]].sum(0) for i in range(len(lines))])
n = np.array([cnt[[i * NGROUP + k for k in range(0, 1 + NBINS)]].sum() for i in range(len(lines))])
tot = P.sum(0); keep = (1e6 * tot / tot.sum() >= 10)
keep[g[g.gene_symbol.isin(["XIST"] + Y_GENES)].token_id.values] = False
valid = np.zeros(P.shape[1], bool); valid[g.token_id.values] = True; keep &= valid
L = np.log2(1e6 * P[:, keep] / P[:, keep].sum(1, keepdims=True) + 1)
R = np.corrcoef(L); np.fill_diagonal(R, np.nan)
rows = []
for i, l in enumerate(lines):
    if n[i] < 50:
        continue
    j = int(np.nanargmax(R[i])); r2 = np.sort(R[i][~np.isnan(R[i])])[-2]
    rows.append(dict(cellosaurus_id=l, cell_name=name[l], n_deep_cells=int(n[i]), nearest_line=name[lines[j]],
                     r_nearest=R[i, j], r_second=r2, gap_nearest_second=R[i, j] - r2))
out = pd.DataFrame(rows).sort_values("r_nearest", ascending=False)
allr = R[np.triu_indices(len(lines), 1)]
out.attrs = {}
out.to_csv(os.path.join(ROOT, "results/F1/F1_exploratory_line_similarity.tsv"), sep="\t", index=False)
pd.DataFrame([dict(n_pairs=int(np.isfinite(allr).sum()), r_median=float(np.nanmedian(allr)), r_p95=float(np.nanpercentile(allr, 95)),
                   r_max=float(np.nanmax(allr)), n_genes=int(keep.sum()))]).to_csv(
    os.path.join(ROOT, "results/F1/F1_exploratory_line_similarity_summary.tsv"), sep="\t", index=False)
print(out.head(8).round(3).to_string()); print(pd.read_csv(os.path.join(ROOT, "results/F1/F1_exploratory_line_similarity_summary.tsv"), sep="\t"))
