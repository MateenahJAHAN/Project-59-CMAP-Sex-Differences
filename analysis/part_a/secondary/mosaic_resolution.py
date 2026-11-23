"""Resolve the mosaic calls of the two Tahoe callers by the written rules; report line-level sensitivity differences.

Written rules applied (analysis_plan/PartA_frozen_rules.md, commit 7cac9a4):
  line 9:  disagreements are resolved by a written rule ("ambiguous" if not resolvable), never by choosing the preferred call.
  line 52: report every share as the range across (a) to (c), never as one number.
  line 81: if A549 fails the mosaic rule on full data, say so, give the full-data share, and drop the mosaic claim.
Resolution used here: a line is "mosaic" only if both callers call it mosaic; if they differ, the call is
"not established" and the mosaic claim is dropped. Shares are the union range over both callers' estimators.
"""
import json
import sys
from pathlib import Path
import pandas as pd

OUT = Path(sys.argv[1]).resolve() / "results/consensus"
mos = pd.read_csv(OUT / "k9_mosaic_comparison.tsv", sep="\t")
cmp = pd.read_csv(OUT / "k9_line_comparison.tsv", sep="\t")

mos["verdict_final"] = [
    "mosaic" if (a and b) else ("not mosaic" if (not a and not b) else "not established (callers differ)")
    for a, b in zip(mos["mosaic_loc"], mos["mosaic_f2"])]
mos["share_lo_both"] = mos[["share_range_lo", "share_range_min"]].min(axis=1)
mos["share_hi_both"] = mos[["share_range_hi", "share_range_max"]].max(axis=1)
keep = ["cell_name", "marker", "mosaic_loc", "mosaic_f2", "verdict_final",
        "f_margin_over_best_other", "f_margin_f2", "best_other_line", "f_best_other_line",
        "pass_b", "rule_b", "pass_d", "rule_d", "pass_e", "rule_e", "pass_f", "rule_f",
        "share_range_lo", "share_range_hi", "share_range_min", "share_range_max", "share_lo_both", "share_hi_both",
        "f_single_neg_share_own_line", "nc_neg_own_share"]
res = mos[keep].rename(columns={
    "f_margin_over_best_other": "f_margin_local", "best_other_line": "best_other_line_local",
    "f_best_other_line": "best_other_line_f2", "share_range_lo": "share_lo_local", "share_range_hi": "share_hi_local",
    "share_range_min": "share_lo_f2", "share_range_max": "share_hi_f2",
    "f_single_neg_share_own_line": "single_neg_nearest_own_f2", "nc_neg_own_share": "single_neg_nearest_own_local"})
res.to_csv(OUT / "mosaic_resolution.tsv", sep="\t", index=False, float_format="%.6f")

ll = cmp.loc[~cmp["state_linelevel_same"], ["cell_name", "state_loc", "state_linelevel_loc", "state_linelevel_f2"]]
ll.to_csv(OUT / "linelevel_sensitivity_differences.tsv", sep="\t", index=False)

lines = []
lines.append("Mosaic verdicts after resolution:")
lines.append(res["verdict_final"].value_counts().to_string())
lines.append("")
lines.append(res[res["verdict_final"] != "not mosaic"][["cell_name", "marker", "verdict_final", "f_margin_local", "f_margin_f2",
            "share_lo_both", "share_hi_both", "single_neg_nearest_own_local", "single_neg_nearest_own_f2"]].to_string(index=False))
for ctrl in ["NCI-H460", "LoVo"]:
    r = res[(res.cell_name == ctrl) & (res.marker == "Y")].iloc[0]
    lines.append(f"control {ctrl}: verdict {r.verdict_final}; share {r.share_lo_both:.3f} to {r.share_hi_both:.3f}; "
                 f"f margin local {r.f_margin_local:.4f}, F2 {r.f_margin_f2:.4f}")
lines.append("")
lines.append("Line-level (sensitivity) cut-offs: lines whose line-level state differs between callers:")
lines.append(ll.to_string(index=False) if len(ll) else "none")
txt = "\n".join(lines)
print(txt)
(OUT / "mosaic_resolution_log.txt").write_text(txt + "\n")
