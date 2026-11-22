"""Analysis 5 (../PLAN.md): E3 for every panel, term and scale, read from results/G/G_H3.tsv (supplementary table).

Adds the HC3 standard error implied by each interval (95% half-width / t quantile on the residual df, in percent of
mean S), the 95% half-width, and which rows make claims. The interval method is the frozen one (gcore.scalar_fit2,
ganalyses.express): OLS of each line's mean response size S on lineage fixed effects, TP53 and the tested terms; HC3
standard error of the tested coefficient; estimate +- t(df, 0.95 or 0.975) x SE; all multiplied by 100 / mean S of the
model lines, with mean S treated as fixed. Writes results/e3_all_panels_scales.tsv.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as K  # noqa: E402

import numpy as np
import pandas as pd
from scipy import stats

TERM_NAME = {"factY": "Y term (lines that express Y genes vs neither marker)",
             "factX": "X term (lines that express XIST vs neither marker)",
             "annotated": "annotated sex (male vs female label)"}


def claims(r):
    if r.role == "viability-adjusted":
        return "description (PRISM AUC as a covariate)"
    if r.panel != "P1":
        return "no claim (P1b and P2 make no claims of their own)"
    if r.role == "S12":
        return "sensitivity analysis S12 (level 5), no claim"
    if r.term == "annotated":
        return "H4, reported beside H3, no claim of its own"
    return "H3 (primary, P1)"


def main():
    src = os.path.join(K.GRES, "G_H3.tsv")
    h = pd.read_csv(src, sep="\t", float_precision="round_trip")
    rows = []
    for _, r in h.iterrows():
        df = int(r["df"])
        hw95 = (float(r.hi95) - float(r.lo95)) / 2
        hw90 = (float(r.hi90) - float(r.lo90)) / 2
        se = hw95 / stats.t.ppf(0.975, df)
        se90 = hw90 / stats.t.ppf(0.95, df)
        assert abs(se - se90) <= 1e-9 * max(1.0, abs(se)), "90% and 95% intervals imply different SEs"
        rows.append(dict(panel=r.panel, term=r.term, term_name=TERM_NAME[r.term], scale=r.scale, tensor=r.tensor, role=r.role,
                         claim=claims(r), n_lines=int(r.n_lines), df=df, mean_S=float(r.mean_S), estimate=float(r.estimate),
                         lo90=float(r.lo90), hi90=float(r.hi90), lo95=float(r.lo95), hi95=float(r.hi95),
                         halfwidth95=hw95, se_hc3_pct=se, difference_declared=bool(r.difference_declared),
                         sign=("negative" if r.estimate < 0 else "positive"), unit=r.unit, interval=r.interval))
    out = pd.DataFrame(rows)
    out["interval_method"] = ("OLS of each line's mean response size S on lineage fixed effects, TP53 and the tested terms; "
                              "HC3 standard error of the tested coefficient; estimate +- t(df) quantile x SE (90%: 0.95; "
                              "95%: 0.975); x 100 / mean S of the model lines (mean S treated as fixed)")
    out["source_file"] = src
    out.to_csv(os.path.join(K.RES, "e3_all_panels_scales.tsv"), sep="\t", index=False)
    print(out[["panel", "term", "scale", "role", "n_lines", "df", "estimate", "lo95", "hi95", "se_hc3_pct"]].to_string(index=False))


if __name__ == "__main__":
    main()
