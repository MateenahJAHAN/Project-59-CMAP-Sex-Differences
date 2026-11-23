"""Step 5: ambient Y and XIST per plate; foreign-marker cells per line against plate background and abundance."""
import json, sys, time
import numpy as np, pandas as pd
from scipy.stats import poisson, pearsonr, spearmanr
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import runlog

t0 = time.time()
R, D = "results/F2", "data/F2"
a = pd.read_parquet(f"{D}/cells_analysis.parquet")
st = pd.read_csv(f"{R}/F2_states.tsv", sep="\t").set_index("cell_name")
cut = pd.read_csv(f"{R}/F2_cutoffs.tsv", sep="\t")
cf = cut[cut.cutoff_set == "final_line_x_plate"].set_index("score").cutoff
yneg_final = set(st.index[st.y_score <= cf["y_score"]])
xneg_final = set(st.index[st.xist <= cf["xist"]])
yneg_p1 = set(st.index[st.y_negative_pass1])
xneg_p1 = set(st.index[st.xist_negative_pass1])
rows = []
for plate, g in a.groupby("plate"):
    gy = g[g.line.isin(yneg_final)]
    gx = g[g.line.isin(xneg_final)]
    gy0, gx0 = gy[~gy.foreign], gx[~gx.foreign]
    rows.append(dict(plate=plate, y_neg_lines=len(set(gy.line)), y_neg_cells_clean=len(gy0), y_neg_umis_clean=int(gy0.umis.sum()),
                     y_umis_clean=int(gy0.Y8.sum()), y_per_1000_umis=1000 * gy0.Y8.sum() / gy0.umis.sum(),
                     y_per_1000_umis_incl_foreign=1000 * gy.Y8.sum() / gy.umis.sum(),
                     y_cells_with_any_y_share=float((gy0.Y8 > 0).mean()),
                     xist_neg_lines=len(set(gx.line)), xist_neg_cells_clean=len(gx0), xist_umis_clean=int(gx0.XIST.sum()),
                     xist_per_1000_umis=1000 * gx0.XIST.sum() / gx0.umis.sum(),
                     xist_per_1000_umis_incl_foreign=1000 * gx.XIST.sum() / gx.umis.sum()))
amb = pd.DataFrame(rows)
amb["plate_n"] = amb.plate.str.replace("plate", "").astype(int)
amb = amb.sort_values("plate_n")
amb.to_csv(f"{R}/F2_ambient.tsv", sep="\t", index=False)
print(amb[["plate", "y_neg_lines", "y_neg_cells_clean", "y_per_1000_umis", "y_per_1000_umis_incl_foreign", "xist_neg_lines", "xist_per_1000_umis"]].to_string())
bgY = dict(zip(amb.plate, amb.y_per_1000_umis / 1000))
bgX = dict(zip(amb.plate, amb.xist_per_1000_umis / 1000))

fr = []
for line, g in a.groupby("line"):
    r = dict(cell_name=line, n_vehicle_cells=len(g), log10_cells=np.log10(len(g)),
             y_negative_pass1=line in yneg_p1, xist_negative_pass1=line in xneg_p1)
    if line in yneg_p1:
        lam = g.plate.map(bgY).to_numpy() * g.umis.to_numpy()
        r.update(foreign_Y_cells=int(g.foreign_Y.sum()), foreign_Y_share=float(g.foreign_Y.mean()),
                 expected_foreign_Y_cells_from_background=float(poisson.sf(2, lam).sum()))
    if line in xneg_p1:
        lam = g.plate.map(bgX).to_numpy() * g.umis.to_numpy()
        r.update(foreign_XIST_cells=int(g.foreign_XIST.sum()), foreign_XIST_share=float(g.foreign_XIST.mean()),
                 expected_foreign_XIST_cells_from_background=float(poisson.sf(2, lam).sum()))
    fr.append(r)
fr = pd.DataFrame(fr)
fr["foreign_Y_obs_over_expected"] = fr.foreign_Y_cells / fr.expected_foreign_Y_cells_from_background
fr["foreign_XIST_obs_over_expected"] = fr.foreign_XIST_cells / fr.expected_foreign_XIST_cells_from_background
fr.to_csv(f"{R}/F2_foreign_cells.tsv", sep="\t", index=False)
summ = []
for mk in ["Y", "XIST"]:
    s = fr[fr[f"foreign_{mk}_share"].notna()]
    pr = pearsonr(s.log10_cells, s[f"foreign_{mk}_share"])
    sr = spearmanr(s.log10_cells, s[f"foreign_{mk}_share"])
    summ.append(dict(marker=mk, n_lines=len(s), total_foreign_cells=int(s[f"foreign_{mk}_cells"].sum()),
                     total_expected_from_background=float(s[f"expected_foreign_{mk}_cells_from_background"].sum()),
                     median_share=float(s[f"foreign_{mk}_share"].median()), min_share=float(s[f"foreign_{mk}_share"].min()),
                     max_share=float(s[f"foreign_{mk}_share"].max()),
                     pearson_r_share_vs_log10_cells=pr.statistic, pearson_p=pr.pvalue,
                     spearman_rho_share_vs_log10_cells=sr.statistic, spearman_p=sr.pvalue))
# exploratory (not pre-registered): Y foreign share without lines whose arm state shows partial Y (Yq lost or Yp lost)
part = set(st.index[st.arm_state.isin(["Yq lost", "Yp lost"])])
s = fr[fr.foreign_Y_share.notna() & ~fr.cell_name.isin(part)]
pr = pearsonr(s.log10_cells, s.foreign_Y_share); sr = spearmanr(s.log10_cells, s.foreign_Y_share)
summ.append(dict(marker="Y, exploratory: without partial-Y lines " + ",".join(sorted(part & set(fr[fr.foreign_Y_share.notna()].cell_name))), n_lines=len(s),
                 total_foreign_cells=int(s.foreign_Y_cells.sum()), total_expected_from_background=float(s.expected_foreign_Y_cells_from_background.sum()),
                 median_share=float(s.foreign_Y_share.median()), min_share=float(s.foreign_Y_share.min()), max_share=float(s.foreign_Y_share.max()),
                 pearson_r_share_vs_log10_cells=pr.statistic, pearson_p=pr.pvalue, spearman_rho_share_vs_log10_cells=sr.statistic, spearman_p=sr.pvalue))
pd.DataFrame(summ).to_csv(f"{R}/F2_foreign_summary.tsv", sep="\t", index=False)
print(pd.DataFrame(summ).to_string())
runlog("s06_ambient", t0)
