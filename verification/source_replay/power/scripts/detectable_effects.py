"""Analysis 6 (../PLAN.md): effects detectable with 80% power, from the observed interval widths of the real run.

For each interval, SE = half-width / t quantile on its residual df (HC3, as the frozen code computes it). The effect
detectable with 80% power by a two-sided test at 0.05 is (t(df, 0.975) + t(df, 0.80)) x SE, and the approximate power
at the plan's margin m is P(|T| > t(df, 0.975)) with T noncentral t (df, m / SE). E3: results/G/G_H3.tsv (P1, both
scales). E2 and E1: results/G/G_H2_per_compound.tsv (P1, tested compounds), per compound, then the median over
compounds (E1 uses each compound's median-gene HC3 half-width in SD units). These are design approximations; the
planted-effect simulation (power.py) is the direct estimate. Writes results/detectable_effects.tsv.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as K  # noqa: E402

import numpy as np
import pandas as pd
from scipy import stats


def mde(se, df):
    return (stats.t.ppf(0.975, df) + stats.t.ppf(0.80, df)) * se


def power_at(m, se, df):
    tc = stats.t.ppf(0.975, df)
    ncp = m / se
    return 1 - stats.nct.cdf(tc, df, ncp) + stats.nct.cdf(-tc, df, ncp)


def main():
    cfg = K.gconfig.load(os.path.join(K.REPO, "config", "G_defaults.json"))
    m3, m2, m1 = cfg.const("E3_margin_pct"), cfg.const("E2_sesoi_pct"), cfg.const("E1_sesoi_sd")
    h3 = pd.read_csv(os.path.join(K.GRES, "G_H3.tsv"), sep="\t", float_precision="round_trip")
    h2 = pd.read_csv(os.path.join(K.GRES, "G_H2_per_compound.tsv"), sep="\t", float_precision="round_trip")
    rows = []
    for term in ("factY", "factX", "annotated"):
        for scale, role in (("level3", "primary"), ("level5", "S12")):
            r = h3[(h3.panel == "P1") & (h3.term == term) & (h3.scale == scale) & (h3.role == role)].iloc[0]
            df = int(r["df"])
            se = (r.hi95 - r.lo95) / 2 / stats.t.ppf(0.975, df)
            rows.append(dict(estimand="E3", term=term, scale=scale, unit="% of mean S", summary="single interval", n=1, df=df,
                             se=se, mde80=mde(se, df), margin=m3, power_at_margin=power_at(m3, se, df),
                             source="results/G/G_H3.tsv"))
        t = h2[(h2.panel == "P1") & (h2.term == term) & h2.in_ranking & h2.E2_est.notna()]
        df = t["df"].values.astype(float)
        se2 = (t.E2_hi90.values - t.E2_lo90.values) / 2 / stats.t.ppf(0.95, df)
        rows.append(dict(estimand="E2", term=term, scale="level5", unit="% of the compound's mean s",
                         summary="median over tested compounds", n=len(t), df=float(np.median(df)), se=float(np.median(se2)),
                         mde80=float(np.median(mde(se2, df))), margin=m2, power_at_margin=float(np.median(power_at(m2, se2, df))),
                         mean_power_at_margin=float(np.mean(power_at(m2, se2, df))), source="results/G/G_H2_per_compound.tsv"))
        se1 = t.E1_hw90_hc3_sd_median.values / stats.t.ppf(0.95, df)
        rows.append(dict(estimand="E1", term=term, scale="level5", unit="SD per gene",
                         summary="median over tested compounds of the median-gene HC3 half-width", n=len(t), df=float(np.median(df)),
                         se=float(np.median(se1)), mde80=float(np.median(mde(se1, df))), margin=m1,
                         power_at_margin=float(np.median(power_at(m1, se1, df))), mean_power_at_margin=float(np.mean(power_at(m1, se1, df))),
                         source="results/G/G_H2_per_compound.tsv"))
    out = pd.DataFrame(rows)
    out["method"] = ("SE = interval half-width / t(df) quantile; detectable effect at 80% power, two-sided 0.05: "
                     "(t(df, 0.975) + t(df, 0.80)) x SE; power at the margin from a noncentral t; per-compound values summarised by the median")
    out.to_csv(os.path.join(K.RES, "detectable_effects.tsv"), sep="\t", index=False)
    print(out[["estimand", "term", "scale", "n", "df", "se", "mde80", "margin", "power_at_margin"]].to_string(index=False))


if __name__ == "__main__":
    main()
