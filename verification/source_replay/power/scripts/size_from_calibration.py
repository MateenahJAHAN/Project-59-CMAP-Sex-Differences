"""Analysis 1a (../PLAN.md): size of the panel-level test under pseudo-labels, from the saved calibration rows.

Reads results/G/G_calibration.tsv (one row per pseudo-labelling, configured runs, B = 2,000) and writes
results/size_panel_test.tsv: per panel and term, the number and share of pseudo-labellings with panel p below 0.05 and
below 0.025, the binomial standard error, the exact (Clopper-Pearson) 95% interval, the mean and median panel p, and the
number of pseudo-labellings whose panel p lies within 2 permutation standard errors of each threshold.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as K  # noqa: E402

import numpy as np
import pandas as pd


def main():
    src = os.path.join(K.GRES, "G_calibration.tsv")
    cal = pd.read_csv(src, sep="\t", float_precision="round_trip")
    cal = cal[cal.run == "configured"]
    summ = pd.read_csv(os.path.join(K.GRES, "G_calibration_summary.tsv"), sep="\t", float_precision="round_trip")
    rows = []
    for (pn, term), g in cal.groupby(["panel", "term"], sort=False):
        p = g.panel_p.values.astype(float)
        p = p[np.isfinite(p)]
        n = len(p)
        B = int(g.B.iloc[0])
        r = dict(panel=pn, term=term, n_pseudo=n, B=B, source=g.source.iloc[0], mean_panel_p=float(np.mean(p)),
                 median_panel_p=float(np.median(p)))
        for a, lab in ((0.05, "005"), (0.025, "0025")):
            x = int((p < a).sum())
            lo, hi = K.clopper_pearson(x, n)
            se_perm = float(np.sqrt(a * (1 - a) / B))
            r.update({f"n_p_lt_{lab}": x, f"share_p_lt_{lab}": x / n, f"mc_se_{lab}": float(np.sqrt((x / n) * (1 - x / n) / n)),
                      f"cp95_lo_{lab}": lo, f"cp95_hi_{lab}": hi,
                      f"n_within_2se_perm_{lab}": int((np.abs(p - a) <= 2 * se_perm).sum()),
                      f"reading_{lab}": ("conservative" if hi < a else ("liberal" if lo > a else "consistent with nominal"))})
        s = summ[(summ.panel == pn) & (summ.term == term)]
        r["share_pseudo_panel_p_lt_005_in_summary"] = float(s.share_pseudo_panel_p_lt_005.iloc[0]) if len(s) else np.nan
        r["mean_share_compounds_p_lt_005"] = float(s.mean_share_p_lt_005.iloc[0]) if len(s) else np.nan
        rows.append(r)
    out = pd.DataFrame(rows)
    assert np.allclose(out.share_p_lt_005, out.share_pseudo_panel_p_lt_005_in_summary, rtol=0, atol=1e-12), \
        "recount differs from G_calibration_summary.tsv"
    out["source_file"] = src
    out.to_csv(os.path.join(K.RES, "size_panel_test.tsv"), sep="\t", index=False)
    print(out[["panel", "term", "n_pseudo", "n_p_lt_005", "n_p_lt_0025", "cp95_hi_005", "cp95_hi_0025", "mean_panel_p"]].to_string(index=False))


if __name__ == "__main__":
    main()
