"""Step 6: summary counts quoted in the report, written to results/F1/F1_report_numbers.tsv (key, value, source)."""
import sys, os
import numpy as np, pandas as pd

ROOT = sys.argv[1]
RES = os.path.join(ROOT, "results/F1")
rows = []


def put(key, value, source, note=""):
    rows.append(dict(key=key, value=value, source=source, note=note))


st = pd.read_csv(os.path.join(RES, "F1_emerald_states.tsv"), sep="\t")
S = "F1_emerald_states.tsv"
put("emerald_n_lines", len(st), S)
for k, v in st.state.value_counts().items():
    put(f"emerald_state_{k}", v, S, "column state")
for k, v in st.confidence.value_counts().items():
    put(f"emerald_confidence_{k}", v, S, "column confidence")
for k, v in st.arm_state.value_counts().items():
    put(f"emerald_arm_{k}", v, S, "column arm_state")
for k, v in st.annotation.value_counts().items():
    put(f"emerald_annotation_{k}", v, S, "column annotation")
mf = st[st.annotation.isin(["Male", "Female"])]
usual = ((mf.annotation == "Male") & (mf.state == "XY")) | ((mf.annotation == "Female") & (mf.state == "XX-like"))
put("emerald_male_female_lines", len(mf), S, "annotation Male or Female")
put("emerald_male_female_usual_pattern", int(usual.sum()), S, "Male&XY or Female&XX-like (descriptive, not pre-registered)")
put("emerald_male_female_not_usual_pattern", int((~usual).sum()), S, "descriptive, not pre-registered")
hm = mf[mf.confidence.isin(["high", "medium"])]
u2 = ((hm.annotation == "Male") & (hm.state == "XY")) | ((hm.annotation == "Female") & (hm.state == "XX-like"))
put("emerald_male_female_highmed_lines", len(hm), S); put("emerald_male_female_highmed_not_usual", int((~u2).sum()), S)
put("emerald_lines_ge1000_vehicle", int((st.n_vehicle_cells >= 1000).sum()), S)
put("emerald_lines_lt10_vehicle_per_sample_all", int((st.n_plates == 0).sum()), S, "n_plates == 0: no sample with >= 10 vehicle cells")
yneg = st[st.y_negative_line]
put("emerald_y_negative_lines", len(yneg), S, "y_negative_line")
put("emerald_x_negative_lines", int(st.x_negative_line.sum()), S, "x_negative_line")
fy = yneg[yneg.n_vehicle_cells >= 100].foreign_share_y
put("emerald_foreign_share_y_median_lines_ge100", float(fy.median()), S, "Y-negative lines with >= 100 vehicle cells")
fy2 = yneg[(yneg.n_vehicle_cells >= 100) & (yneg.cell_name != "NCI-H2030")].foreign_share_y
put("emerald_foreign_share_y_max_excl_NCI-H2030", float(fy2.max()), S, "Y-negative lines with >= 100 vehicle cells, NCI-H2030 excluded")
put("emerald_foreign_share_y_expected_median", float(yneg[yneg.n_vehicle_cells >= 100].foreign_share_y_expected_ambient.median()), S)
xn = st[st.x_negative_line & (st.n_vehicle_cells >= 100)]
put("emerald_foreign_share_xist_median_lines_ge100", float(xn.foreign_share_xist.median()), S)
put("emerald_foreign_share_xist_max_lines_ge100", float(xn.foreign_share_xist.max()), S)

m = pd.read_csv(os.path.join(RES, "F1_emerald_mosaic.tsv"), sep="\t")
M = "F1_emerald_mosaic.tsv"
t = m[m.role == "test"]
put("emerald_y_tested_lines", int((t.marker == "Y").sum()), M, "role == test, marker Y")
put("emerald_xist_tested_lines", int((t.marker == "XIST").sum()), M)
put("emerald_mosaic_flagged", ";".join(t[t.mosaic_flag == "mosaic"].cell_name + "(" + t[t.mosaic_flag == "mosaic"].marker + ")"), M)
put("emerald_mosaic_flagged_n", int((t.mosaic_flag == "mosaic").sum()), M)
for c in ["a_holds", "b_holds", "c_holds", "d_holds"]:
    tt = t[t[c].notna()]
    put(f"emerald_tests_{c}_true", f"{int((tt[c] == True).sum())}/{len(tt)}", M, "among tested lines where evaluable")
ty = t[t.marker == "Y"].sort_values("pi_neg", ascending=False)
put("emerald_y_tested_pi_ranked", ";".join(f"{a}:{b:.3f}" for a, b in zip(ty.cell_name, ty.pi_neg)), M)
ctrl = m[(m.role == "control") & (m.marker == "Y")]
ctrl_neg = ctrl.merge(st[["cell_name", "y_negative_line"]], on="cell_name")
ctrl_neg = ctrl_neg[ctrl_neg.y_negative_line]
put("emerald_y_control_ynegative_lines_n", len(ctrl_neg), M, "Y EM on Y-negative lines with >= 1,000 vehicle cells")
put("emerald_y_control_pi_min", float(ctrl_neg.pi_neg.min()), M); put("emerald_y_control_pi_median", float(ctrl_neg.pi_neg.median()), M)
put("emerald_y_control_pi_max", float(ctrl_neg.pi_neg.max()), M)
put("emerald_y_control_pi_in_15_85", int(ctrl_neg.pi_neg.between(0.15, 0.85).sum()), M)
put("emerald_y_control_rate_per_1k_min", float(ctrl_neg.marker_umi_per_1k.min()), M)
put("emerald_y_control_rate_per_1k_max", float(ctrl_neg.marker_umi_per_1k.max()), M)
ms = pd.read_csv(os.path.join(RES, "F1_emerald_mosaic_sensitivity_background.tsv"), sep="\t")
MS = "F1_emerald_mosaic_sensitivity_background.tsv"
cs = ms[ms.role == "control"].merge(st[["cell_name", "y_negative_line"]], on="cell_name")
cs = cs[cs.y_negative_line & (cs.cell_name != "NCI-H2030")]
put("emerald_sens_y_control_n", len(cs), MS, "Y-negative control lines, NCI-H2030 excluded")
put("emerald_sens_y_control_pi_min", float(cs.pi_neg.min()), MS); put("emerald_sens_y_control_pi_median", float(cs.pi_neg.median()), MS)
put("emerald_sens_y_control_pi_in_15_85", int(cs.pi_neg.between(0.15, 0.85).sum()), MS)
ts = ms[ms.role == "test"]
put("emerald_sens_mosaic_flagged", ";".join(ts[ts.mosaic_flag == "mosaic"].cell_name), MS)

ty2 = t[t.marker == "Y"]
ratio = ty2.lam_per_1k / ty2.background_per_1k
put("emerald_y_tests_lambda_over_background_range", f"{ratio.min():.1f}-{ratio.max():.1f}", M, "Y test rows: lam_per_1k / background_per_1k")
put("emerald_y_tests_or_x_matched_range", f"{ty2.or_x_matched.min():.3f}-{ty2.or_x_matched.max():.3f}", M)
put("emerald_y_tests_or_x_unthinned_range", f"{ty2.or_x_unthinned_5050.min():.3f}-{ty2.or_x_unthinned_5050.max():.3f}", M)
put("emerald_y_control_rate_max_excl_NCI-H2030", float(ctrl_neg[ctrl_neg.cell_name != "NCI-H2030"].marker_umi_per_1k.max()), M)
bgmin = m[m.marker == "Y"].background_per_1k.min()
put("emerald_y_control_lines_rate_below_min_background", f"{int((ctrl_neg.marker_umi_per_1k < bgmin).sum())}/{len(ctrl_neg)}", M,
    f"control rate < smallest primary Y background ({bgmin:.4f} per 1k)")
sm = pd.read_csv(os.path.join(RES, "F1_sciplex_mosaic.tsv"), sep="\t")
spt = sm[sm.role == "test"]
put("sciplex_tests_or_x_matched_range", f"{spt[spt.marker == 'Y'].or_x_matched.min():.3f}-{spt[spt.marker == 'Y'].or_x_matched.max():.3f}", "F1_sciplex_mosaic.tsv")
put("sciplex_mean_depth_range", f"{sm.mean_depth.min():.0f}-{sm.mean_depth.max():.0f}", "F1_sciplex_mosaic.tsv", "mean UMIs per vehicle cell, all rows")
xt = spt[spt.marker == "XIST"]
put("sciplex_xist_tests_expected_xist_umis_per_cell", ";".join(f"{a}/{b}:{r * d / 1000:.2f}" for a, b, r, d in zip(xt.dataset, xt.line, xt.marker_umi_per_1k, xt.mean_depth)),
    "F1_sciplex_mosaic.tsv", "marker_umi_per_1k x mean_depth / 1000")
for tag in ("sciPlex2", "sciPlex3", "sciPlex4"):
    r = sm[(sm.dataset == tag) & (sm.line == "A549") & (sm.marker == "Y")].iloc[0]
    put(f"{tag}_A549_pi", float(r.pi_neg), "F1_sciplex_mosaic.tsv"); put(f"{tag}_A549_flag", r.mosaic_flag, "F1_sciplex_mosaic.tsv")

amb = pd.read_csv(os.path.join(RES, "F1_ambient.tsv"), sep="\t")
b = amb[(amb.level == "batch_background")]
for (ds, mk), g in b.groupby(["dataset", "marker"]):
    put(f"ambient_{ds}_{mk}_per_1k_range", f"{g.rate_per_1k.min():.4f}-{g.rate_per_1k.max():.4f}", "F1_ambient.tsv", "level == batch_background")
eb = b[(b.dataset == "EmeraldBay") & (b.marker == "Y")].rate_per_1k
spb = b[(b.dataset.isin(["sciPlex3", "sciPlex4"])) & (b.marker == "Y")].rate_per_1k
put("ambient_y_emerald_over_sciplex_fold_range", f"{eb.min() / spb.max():.1f}-{eb.max() / spb.min():.1f}", "F1_ambient.tsv",
    "EmeraldBay primary Y background / sciPlex3-4 Y background")
pd.DataFrame(rows).to_csv(os.path.join(RES, "F1_report_numbers.tsv"), sep="\t", index=False)
print(pd.DataFrame(rows).to_string())
