"""K9 cross-check: local Part A caller against cloud task F2, Tahoe-100M vehicle cells.

Rule (analysis_plan/PartA_frozen_rules.md, line 9, committed 7cac9a4 before either caller ran):
  the two tables must give the same state for at least 48 of 50 lines; every disagreement is
  reported and resolved by a written rule ("ambiguous" if not resolvable), never by choosing
  the preferred call.
Stricter wording in FINAL_RECOMMENDATION.md (K9): same state AND confidence for >= 48 of 50.
Both are evaluated and reported.

The report records the checksums of both regenerated caller outputs.
"""
import hashlib, json, sys
from pathlib import Path
import numpy as np
import pandas as pd

work = Path(sys.argv[1]).resolve()
LOCAL_DIR = work / "results/tahoe/latest"
F2_DIR = work / "results/F2"
OUT = work / "results/consensus"
OUT.mkdir(parents=True, exist_ok=True)
log_lines = []
def log(s=""):
    print(s); log_lines.append(str(s))

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()

inputs = [LOCAL_DIR / "partA_local_states.tsv", F2_DIR / "F2_states.tsv"]
input_hashes = {p.name: sha(p) for p in inputs}

loc = pd.read_csv(LOCAL_DIR / "partA_local_states.tsv", sep="\t")
locL = pd.read_csv(LOCAL_DIR / "partA_local_states_linelevel.tsv", sep="\t")
f2 = pd.read_csv(F2_DIR / "F2_states.tsv", sep="\t")
for name, d, key in [("local", loc, "cellosaurus_id"), ("local_linelevel", locL, "cellosaurus_id"), ("F2", f2, "cell_line_id")]:
    assert len(d) == 50, (name, len(d))
    assert d[key].is_unique, name
    assert d["cell_name"].is_unique, name
log(f"rows: local {len(loc)}, local line-level {len(locL)}, F2 {len(f2)}")

# F2 has two Cellosaurus columns; check they agree with each other first.
assert (f2["cell_line_id"] == f2["Cell_ID_Cellosaur"]).all(), "F2 cell_line_id != Cell_ID_Cellosaur"

m = loc.merge(f2, left_on="cellosaurus_id", right_on="cell_line_id", how="outer",
              suffixes=("_loc", "_f2"), indicator=True, validate="one_to_one")
assert (m["_merge"] == "both").all(), m.loc[m["_merge"] != "both", ["cell_name_loc", "cell_name_f2"]]
assert (m["cell_name_loc"] == m["cell_name_f2"]).all(), "cell names differ after join on Cellosaurus id"
m = m.drop(columns="_merge")
mL = m.merge(locL[["cellosaurus_id", "state", "confidence", "arm_state"]].rename(
    columns={"state": "state_linelevel_loc", "confidence": "confidence_linelevel_loc", "arm_state": "arm_state_linelevel_loc"}),
    on="cellosaurus_id", how="left", validate="one_to_one")
assert len(mL) == 50 and mL["state_linelevel_loc"].notna().all()
m = mL

# plates_agree: local stores a fraction and n_plates_agree; F2 stores a count.
assert (m["n_plates_loc"] == m["n_plates_f2"]).all(), "plate counts differ"

cmp = pd.DataFrame({
    "cell_name": m["cell_name_loc"],
    "cellosaurus_id": m["cellosaurus_id"],
    "label_loc": m["annotation"],
    "label_f2": m["cellosaurus_sex"],
    "state_loc": m["state_loc"], "state_f2": m["state_f2"],
    "confidence_loc": m["confidence_loc"], "confidence_f2": m["confidence_f2"],
    "arm_state_loc": m["arm_state_loc"], "arm_state_f2": m["arm_state_f2"],
    "state_linelevel_loc": m["state_linelevel_loc"], "state_linelevel_f2": m["state_linelevel"],
    "plates_agree_loc": m["n_plates_agree"], "plates_agree_f2": m["plates_agree_f2"],
    "grid_stable_loc": m["grid_stable_loc"], "grid_stable_f2": m["grid_stable_f2"],
    "n_cells_before_filter_loc": m["n_vehicle_cells_before_filter"], "n_cells_before_filter_f2": m["n_vehicle_cells_f2"],
    "n_cells_called_loc": m["n_vehicle_cells_loc"], "n_cells_called_f2": m["n_cells_called"],
    "XIST_loc": m["XIST"], "XIST_f2": m["xist"],
    "y_score_loc": m["y_score_loc"], "y_score_f2": m["y_score_f2"],
    "yp_score_loc": m["yp_score_loc"], "yp_score_f2": m["yp_score_f2"],
    "yq_score_loc": m["yq_score_loc"], "yq_score_f2": m["yq_score_f2"],
})
for c in ["state", "confidence", "arm_state", "state_linelevel", "grid_stable"]:
    cmp[f"{c}_same"] = cmp[f"{c}_loc"].astype(str) == cmp[f"{c}_f2"].astype(str)
cmp["plates_agree_same"] = cmp["plates_agree_loc"].astype(int) == cmp["plates_agree_f2"].astype(int)
cmp["label_same"] = cmp["label_loc"].astype(str) == cmp["label_f2"].astype(str)
cmp["n_cells_called_same"] = cmp["n_cells_called_loc"] == cmp["n_cells_called_f2"]
cmp["n_cells_before_same"] = cmp["n_cells_before_filter_loc"] == cmp["n_cells_before_filter_f2"]
for s in ["XIST", "y_score", "yp_score", "yq_score"]:
    cmp[f"{s}_absdiff"] = (cmp[f"{s}_loc"] - cmp[f"{s}_f2"]).abs()

n_state = int(cmp["state_same"].sum())
n_state_conf = int((cmp["state_same"] & cmp["confidence_same"]).sum())
log()
log(f"K9 primary (frozen file): same state for {n_state} of 50 lines -> {'PASS' if n_state >= 48 else 'FAIL'} (needs >= 48)")
log(f"K9 strict (FINAL_RECOMMENDATION): same state and confidence for {n_state_conf} of 50 -> {'PASS' if n_state_conf >= 48 else 'FAIL'}")
for c in ["arm_state", "state_linelevel", "plates_agree", "grid_stable", "label", "n_cells_called", "n_cells_before"]:
    log(f"  same {c}: {int(cmp[c + '_same'].sum())} of 50")
for s in ["XIST", "y_score", "yp_score", "yq_score"]:
    log(f"  max |{s} local - F2|: {cmp[s + '_absdiff'].max():.6f}  (values rounded to 4 dp in the local table)")
dis = cmp.loc[~(cmp["state_same"] & cmp["confidence_same"] & cmp["arm_state_same"])]
log(f"lines with any difference in state, confidence or arm state: {len(dis)}")
if len(dis):
    log(dis[["cell_name", "state_loc", "state_f2", "confidence_loc", "confidence_f2", "arm_state_loc", "arm_state_f2"]].to_string(index=False))
lab_dis = cmp.loc[~cmp["label_same"], ["cell_name", "label_loc", "label_f2"]]
if len(lab_dis):
    log("label differences (label source, not a state call):")
    log(lab_dis.to_string(index=False))

# Counts by state for both callers
cnt = pd.DataFrame({"local": cmp["state_loc"].value_counts(), "F2": cmp["state_f2"].value_counts()}).fillna(0).astype(int)
log(); log("state counts:"); log(cnt.to_string())
cc = pd.DataFrame({"local": cmp["confidence_loc"].value_counts(), "F2": cmp["confidence_f2"].value_counts()}).fillna(0).astype(int)
log("confidence counts:"); log(cc.to_string())

# Consensus table. Written rule (frozen file line 9): where the callers differ, the state is "ambiguous".
cons = cmp[["cell_name", "cellosaurus_id", "label_f2", "state_loc", "state_f2", "confidence_loc", "confidence_f2",
            "arm_state_loc", "arm_state_f2", "n_cells_called_f2", "XIST_f2", "y_score_f2", "yp_score_f2", "yq_score_f2"]].copy()
cons = cons.rename(columns={"label_f2": "cellosaurus_sex", "n_cells_called_f2": "n_vehicle_cells_called",
                            "XIST_f2": "xist", "y_score_f2": "y_score", "yp_score_f2": "yp_score", "yq_score_f2": "yq_score"})
cons["state_final"] = np.where(cmp["state_same"], cmp["state_f2"], "ambiguous")
lowest = {"high": 2, "medium": 1, "low": 0}
inv = {v: k for k, v in lowest.items()}
cons["confidence_final"] = [cl if cl == cf else inv[min(lowest[cl], lowest[cf])] for cl, cf in zip(cmp["confidence_loc"], cmp["confidence_f2"])]
cons["arm_state_final"] = np.where(cmp["arm_state_same"], cmp["arm_state_f2"], "ambiguous")
cons["callers_agree"] = cmp["state_same"] & cmp["confidence_same"] & cmp["arm_state_same"]
cons = cons.sort_values("cell_name").reset_index(drop=True)

cmp.to_csv(OUT / "k9_line_comparison.tsv", sep="\t", index=False)
cons.to_csv(OUT / "tahoe_partA_states_consensus.tsv", sep="\t", index=False, float_format="%.6f")

# Mosaic: compare the two callers' verdicts (no pre-set agreement rule exists for this; report both).
ml = pd.read_csv(LOCAL_DIR / "partA_local_mosaic.tsv", sep="\t")
mf = pd.read_csv(F2_DIR / "F2_mosaic.tsv", sep="\t")
mm = ml.merge(mf, on=["cell_name", "marker"], how="outer", suffixes=("_loc", "_f2"), indicator=True)
log(); log(f"mosaic tests: local {len(ml)}, F2 {len(mf)}, joined both {int((mm['_merge']=='both').sum())}")
if (mm["_merge"] != "both").any():
    log(mm.loc[mm["_merge"] != "both", ["cell_name", "marker", "_merge"]].to_string(index=False))
both = mm[mm["_merge"] == "both"].copy()
both["mosaic_loc"] = both["mosaic_verdict"].str.startswith("mosaic")
both["mosaic_f2"] = both["mosaic"].astype(bool)
mos = both[["cell_name", "marker", "n_cells_loc", "n_cells_f2", "b_pi_2000", "b_share_2000", "c_pi", "c_share",
            "d_diff_pi_2000", "d_abs_diff_b2000", "e_or_y", "e_or_ctrl", "e_or_yq0_yp0", "e_or_xhalves_control",
            "r_neg_pos", "f_r_neg_vs_pos", "r_split_half", "f_r_split_half_pos", "best_other_line", "f_best_other_line",
            "r_neg_best_other", "f_r_neg_vs_best_other", "f_margin_over_best_other", "pass_b", "rule_b", "pass_d", "rule_d",
            "pass_e", "rule_e", "pass_f", "rule_f", "f_single_neg_share_own_line", "nc_neg_own_share",
            "share_range_lo", "share_range_hi", "share_range_min", "share_range_max", "mosaic_loc", "mosaic_f2"]].copy()
mos["f_margin_f2"] = mos["f_r_neg_vs_pos"] - mos["f_r_neg_vs_best_other"]
mos.to_csv(OUT / "k9_mosaic_comparison.tsv", sep="\t", index=False)
log(f"mosaic verdict same in {int((mos['mosaic_loc'] == mos['mosaic_f2']).sum())} of {len(mos)} tests")
show = mos[(mos["mosaic_loc"]) | (mos["mosaic_f2"]) | (mos["cell_name"].isin(["A549", "SW48", "NCI-H460", "LoVo", "HCT15"]))]
log(show[["cell_name", "marker", "b_pi_2000", "b_share_2000", "c_pi", "c_share", "f_margin_over_best_other", "f_margin_f2",
          "best_other_line", "f_best_other_line", "mosaic_loc", "mosaic_f2"]].to_string(index=False))

summary = {
    "k9_primary_same_state": n_state, "k9_primary_pass": n_state >= 48,
    "k9_strict_same_state_and_confidence": n_state_conf, "k9_strict_pass": n_state_conf >= 48,
    "same_arm_state": int(cmp["arm_state_same"].sum()),
    "same_state_linelevel": int(cmp["state_linelevel_same"].sum()),
    "state_counts_final": cons["state_final"].value_counts().to_dict(),
    "confidence_counts_final": cons["confidence_final"].value_counts().to_dict(),
    "mosaic_same_verdict": int((mos["mosaic_loc"] == mos["mosaic_f2"]).sum()), "mosaic_tests": int(len(mos)),
    "inputs_sha256": input_hashes,
    "consensus_sha256": sha(OUT / "tahoe_partA_states_consensus.tsv"),
}
(OUT / "k9_summary.json").write_text(json.dumps(summary, indent=1, default=str))
log(); log(json.dumps(summary, indent=1, default=str))
(OUT / "k9_log.txt").write_text("\n".join(log_lines) + "\n")
