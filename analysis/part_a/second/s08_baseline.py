"""Steps 7 to 9: baseline description of marker classes (vehicle cells, descriptive), Y dose, ATCC known-answer note.
Arm-level ratios are skipped (no permitted source of gene positions); F2_arm_ratios.tsv records that."""
import html, re, sys, time
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import runlog, Y8

t0 = time.time()
R, D = "results/F2", "data/F2"
a = pd.read_parquet(f"{D}/cells_analysis.parquet")
st = pd.read_csv(f"{R}/F2_states.tsv", sep="\t").set_index("cell_name")
M = pd.read_csv(f"{R}/F2_mosaic.tsv", sep="\t")

# ---------- step 7 ----------
targets = [("A549", "Y"), ("SW48", "XIST"), ("NCI-H460", "Y"), ("LoVo", "Y")]
for _, r in M[M.flagged_b_d_e].iterrows():
    if (r.cell_name, r.marker) not in targets:
        targets.append((r.cell_name, r.marker))
rows = []
for line, mk in targets:
    g = a[a.line == line].copy()
    if not len(g):
        continue
    g["marker"] = g.Y8 if mk == "Y" else g.XIST
    g["cls"] = np.where(g.marker == 0, f"{mk}_zero", f"{mk}_positive")
    g["quintile"] = pd.qcut(g.umis, 5, labels=False) + 1
    g["genes_per_umi"] = g.n_genes / g.umis
    g["mito_frac"] = g.mito_umis / g.umis
    for (q, c), h in g.groupby(["quintile", "cls"]):
        row = dict(cell_name=line, marker=mk, depth_quintile=int(q), marker_class=c, n_cells=len(h),
                   umis_min_quintile=int(g[g.quintile == q].umis.min()), umis_max_quintile=int(g[g.quintile == q].umis.max()),
                   median_umis=float(h.umis.median()), mean_umis=float(h.umis.mean()), mean_genes_per_umi=float(h.genes_per_umi.mean()),
                   mean_mito_frac=float(h.mito_frac.mean()))
        for ph in ["G1", "S", "G2M"]:
            row[f"phase_{ph}_share"] = float((h.phase == ph).mean())
        rows.append(row)
pd.DataFrame(rows).to_csv(f"{R}/F2_subpop_baseline.tsv", sep="\t", index=False)
pd.DataFrame([dict(status="skipped", reason="gene_metadata holds only gene_symbol, ensembl_id and token_id; no permitted source gives gene chromosome positions or arm boundaries (plan Part 3)")]).to_csv(
    f"{R}/F2_arm_ratios.tsv", sep="\t", index=False)

# ---------- step 8 ----------
yd = []
for line in st.index[st.state == "XY"]:
    g = a[(a.line == line) & ~a.foreign]
    p = g[g.Y8 >= 1]
    yd.append(dict(cell_name=line, n_cells=len(g), n_y_positive=len(p), y_positive_share=len(p) / len(g),
                   y_per_1000_umis_y_positive_cells=1000 * p.Y8.sum() / p.umis.sum(), y_per_1000_umis_all_cells=1000 * g.Y8.sum() / g.umis.sum(),
                   arm_state=st.loc[line, "arm_state"]))
yd = pd.DataFrame(yd).sort_values("y_per_1000_umis_y_positive_cells")
yd["relative_dose_vs_median_xy_line"] = yd.y_per_1000_umis_y_positive_cells / yd.y_per_1000_umis_y_positive_cells.median()
yd.to_csv(f"{R}/F2_y_dose.tsv", sep="\t", index=False)

# The independent product-sheet comparison is documented in the supplement.
runlog("s08_baseline", t0)
