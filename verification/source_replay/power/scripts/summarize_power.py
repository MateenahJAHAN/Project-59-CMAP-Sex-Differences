"""Summaries of the planted-effect run (../PLAN.md, analyses 1b and 2): detection rates, coverage, attenuation and the
plotted values of Figure S3, each with its Monte Carlo error.

Reads raw/<tag>_runs.tsv.gz, raw/<tag>_planted.tsv.gz and raw/<tag>_e2.npz. Writes, in results/:
size_null_scenario.tsv, power_summary.tsv, coverage_e3.tsv, coverage_e2.tsv, attenuation_e3.tsv,
e3_power_interpolated.tsv and figS3_power_values.tsv.
Monte Carlo error: for a share over replicates, the binomial SE and the Wilson 95% interval; for a share over
compound-replicate pairs (E2) or a mean of per-replicate shares, the SE of the per-replicate values over replicates and a
normal 95% interval.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as K  # noqa: E402
import argparse

import numpy as np
import pandas as pd

E3_EST = ("e3l3", "e3l5", "e3l3w", "e3l3s", "e3l5w", "e3l5s")
E3_NAME = {"e3l3": "E3 level 3 (primary scale)", "e3l5": "E3 level 5 (S12 scale)",
           "e3l3w": "E3 level 3, weak-template compounds only", "e3l3s": "E3 level 3, strong-template compounds only",
           "e3l5w": "E3 level 5, weak-template compounds only", "e3l5s": "E3 level 5, strong-template compounds only"}


def rate(x, n):
    x, n = int(x), int(n)
    lo, hi = K.wilson(x, n)
    p = x / n if n else np.nan
    return dict(n=n, x=x, rate=p, mc_se=float(np.sqrt(p * (1 - p) / n)) if n else np.nan, lo95=lo, hi95=hi)


def mean_rate(v):
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    n = len(v)
    m = float(np.mean(v)) if n else np.nan
    se = float(np.std(v, ddof=1) / np.sqrt(n)) if n > 1 else np.nan
    return dict(n=n, rate=m, mc_se=se, lo95=max(0.0, m - 1.959964 * se), hi95=min(1.0, m + 1.959964 * se))


def reading(lo, hi, nominal):
    return "consistent with nominal" if lo <= nominal <= hi else ("below nominal" if hi < nominal else "above nominal")


def main():
    import warnings
    warnings.simplefilter("ignore", category=RuntimeWarning)   # all-NaN slices: the compound that is never tested
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="power")
    ap.add_argument("--raw-dir", default=K.RAW)
    ap.add_argument("--res-dir", default=K.RES)
    a = ap.parse_args()
    runs = pd.read_csv(os.path.join(a.raw_dir, f"{a.tag}_runs.tsv.gz"), sep="\t", float_precision="round_trip")
    planted = pd.read_csv(os.path.join(a.raw_dir, f"{a.tag}_planted.tsv.gz"), sep="\t", float_precision="round_trip")
    z = np.load(os.path.join(a.raw_dir, f"{a.tag}_e2.npz"))
    R = int(runs.rep.nunique())
    B = int(runs.B.iloc[0])
    out = lambda name, df: df.to_csv(os.path.join(a.res_dir, name), sep="\t", index=False)
    for c in [c for c in runs.columns if c.endswith("difference_declared")]:
        runs[c] = runs[c].map(lambda v: True if str(v) == "True" else (False if str(v) == "False" else np.nan))

    # ---------------- 1b: size in the null scenario
    rows = []
    for term in K.TERMS:
        g = runs[(runs.term == term) & (runs.scenario == "N0")]
        r = dict(term=term, B=B, n_rep=len(g))
        for al, lab in ((0.05, "005"), (0.025, "0025")):
            x = int((g.panel_p < al).sum())
            lo, hi = K.clopper_pearson(x, len(g))
            r.update({f"n_panel_p_lt_{lab}": x, f"share_panel_p_lt_{lab}": x / len(g), f"mc_se_{lab}": np.sqrt((x / len(g)) * (1 - x / len(g)) / len(g)),
                      f"cp95_lo_{lab}": lo, f"cp95_hi_{lab}": hi, f"reading_{lab}": "conservative" if hi < al else ("liberal" if lo > al else "consistent with nominal")})
        for col, lab in (("n_p_lt_005", "005"), ("n_p_lt_0025", "0025")):
            m = mean_rate(g[col] / g.n_tested)
            r.update({f"compound_share_p_lt_{lab}": m["rate"], f"compound_share_mc_se_{lab}": m["mc_se"]})
        r.update(mean_panel_p=float(g.panel_p.mean()), median_panel_p=float(g.panel_p.median()),
                 share_rep_any_q_lt_01=float((g.n_q_lt_01 > 0).mean()), mean_n_q_lt_01=float(g.n_q_lt_01.mean()),
                 share_rep_E4_met=float(g.rule3_E4.astype(bool).mean()))
        rows.append(r)
    out("size_null_scenario.tsv", pd.DataFrame(rows))

    # ---------------- planted values of E3 (mean shift from N0, same pseudo-labels)
    n0 = runs[runs.scenario == "N0"].set_index(["term", "rep"])
    theta = {}
    for (term, scen), g in runs.groupby(["term", "scenario"]):
        g = g.set_index(["term", "rep"])
        for e in E3_EST:
            sh = (g[f"{e}_estimate"] - n0.loc[g.index, f"{e}_estimate"]).values.astype(float)
            sh = sh[np.isfinite(sh)]
            theta[(term, scen, e)] = (0.0 if scen == "N0" else float(np.mean(sh)),
                                      float(np.std(sh, ddof=1)) if scen != "N0" and len(sh) > 1 else 0.0,
                                      float(np.std(sh, ddof=1) / np.sqrt(len(sh))) if scen != "N0" and len(sh) > 1 else 0.0)

    # ---------------- power summary
    rows = []
    for (term, scen), g in runs.groupby(["term", "scenario"], sort=False):
        base = dict(term=term, scenario=scen, delta=g.delta.iloc[0] if "delta" in g else np.nan,
                    k=g.k.iloc[0] if "k" in g else np.nan, u=g.u.iloc[0] if "u" in g else np.nan, B=B, n_rep=len(g))
        for e, lab in (("e3l3", "E3 level 3: 95% CI excludes 0"), ("e3l5", "E3 level 5: 95% CI excludes 0")):
            v = g[f"{e}_difference_declared"].dropna().astype(bool)
            th = theta[(term, scen, e)]
            rows.append(dict(base, measure=lab, **rate(v.sum(), len(v)), planted_value=th[0], planted_value_sd=th[1],
                             mean_estimate=float(g[f"{e}_estimate"].mean())))
        if scen.startswith("A"):
            continue
        for al, lab in ((0.025, "panel-level test: p < 0.025"), (0.05, "panel-level test: p < 0.05")):
            rows.append(dict(base, measure=lab, **rate((g.panel_p < al).sum(), len(g))))
        rows.append(dict(base, measure="E4 rule 3: more than 10 compounds at p < 0.05", **rate(g.rule3_E4.astype(bool).sum(), len(g))))
        if scen == "N0" or scen.startswith("M"):
            rows.append(dict(base, measure="any compound at q < 0.1", **rate((g.n_q_lt_01 > 0).sum(), len(g)),
                             mean_count=float(g.n_q_lt_01.mean())))
        if scen[0] in "CG":
            rows.append(dict(base, measure="planted compounds at q < 0.1 (share of planted)", **mean_rate(g.n_planted_q_lt_01 / g.n_planted),
                             mean_count=float(g.n_planted_q_lt_01.mean())))
            rows.append(dict(base, measure="at least one planted compound at q < 0.1", **rate((g.n_planted_q_lt_01 > 0).sum(), len(g))))
            rows.append(dict(base, measure="planted compounds meeting H2 rules 1 and 2 (share of planted)",
                             **mean_rate(g.n_planted_rule12 / g.n_planted)))
            rows.append(dict(base, measure="planted compounds at p < 0.05 (share of planted)", **mean_rate(g.n_planted_p_lt_005 / g.n_planted)))
            rows.append(dict(base, measure="other compounds at q < 0.1 (count per replicate)", n=len(g), rate=np.nan,
                             mean_count=float(g.n_other_q_lt_01.mean()), mc_se=float(g.n_other_q_lt_01.std(ddof=1) / np.sqrt(len(g)))))
            pl = planted[(planted.term == term) & (planted.scenario == scen)]
            dec = ((pl.e2_lo95 > 0) | (pl.e2_hi95 < 0)).groupby(pl.rep).mean()
            rows.append(dict(base, measure="planted compounds with E2 95% CI excluding 0 (share of planted)", **mean_rate(dec.values),
                             mean_estimate=float(pl.e2_est.mean()), planted_value=float((pl.e2_est - pl.e2_N0_est).mean())))
            if scen[0] == "G":
                rows.append(dict(base, measure="planting check: mean E1 over module genes of planted compounds (SD units)",
                                 n=len(pl), rate=np.nan, mean_estimate=float(pl.e1_module_mean.mean()),
                                 mc_se=float(pl.e1_module_mean.std(ddof=1) / np.sqrt(len(pl)))))
    ps = pd.DataFrame(rows)
    out("power_summary.tsv", ps)

    # ---------------- coverage of the E3 intervals
    rows = []
    for (term, scen), g in runs.groupby(["term", "scenario"], sort=False):
        if not (scen == "N0" or scen.startswith("M")):
            continue
        for e in E3_EST:
            th, sd, se_th = theta[(term, scen, e)]
            r = dict(term=term, scenario=scen, delta=0.0 if scen == "N0" else K.M_DELTA[scen], estimand=e, estimand_name=E3_NAME[e],
                     B=B, n_rep=len(g), planted_value=th, planted_value_sd=sd, planted_value_mc_se=se_th,
                     mean_estimate=float(g[f"{e}_estimate"].mean()), sd_estimate=float(g[f"{e}_estimate"].std(ddof=1)))
            for lv in (90, 95):
                lo, hi = g[f"{e}_lo{lv}"].values, g[f"{e}_hi{lv}"].values
                ok = np.isfinite(lo) & np.isfinite(hi)
                x = int(((lo <= th) & (th <= hi))[ok].sum())
                rr = rate(x, ok.sum())
                r.update({f"coverage{lv}": rr["rate"], f"coverage{lv}_mc_se": rr["mc_se"], f"coverage{lv}_lo95": rr["lo95"],
                          f"coverage{lv}_hi95": rr["hi95"], f"reading{lv}": reading(rr["lo95"], rr["hi95"], lv / 100),
                          f"median_halfwidth{lv}": float(np.nanmedian((hi - lo) / 2))})
            from scipy import stats as _st
            se_hc3 = (g[f"{e}_hi95"].values - g[f"{e}_lo95"].values) / 2 / _st.t.ppf(0.975, g[f"{e}_df"].values.astype(float))
            r["median_se_hc3"] = float(np.nanmedian(se_hc3))
            r["se_to_sd_ratio"] = float(np.nanmedian(se_hc3) / g[f"{e}_estimate"].std(ddof=1))   # PLAN.md deviation 8
            rows.append(r)
    out("coverage_e3.tsv", pd.DataFrame(rows))

    # ---------------- coverage of the E2 intervals (per compound, level 5)
    E2 = z["E2"]
    scen_e2 = list(z["scenarios"])
    strata = {"all compounds": np.arange(E2.shape[3]), "weak-template compounds (level 5)": z["weak5"],
              "strong-template compounds (level 5)": z["strong5"]}
    rows = []
    for ti, term in enumerate(K.TERMS):
        est0 = E2[:, ti, 0, :, 0]
        for si, scen in enumerate(scen_e2):
            est = E2[:, ti, si, :, 0]
            shift = est - est0
            with np.errstate(invalid="ignore"):
                th_c = np.zeros(E2.shape[3]) if scen == "N0" else np.nanmean(shift, axis=0)
                sd_c = np.zeros(E2.shape[3]) if scen == "N0" else np.nanstd(shift, axis=0, ddof=1)
            for sname, idx in strata.items():
                r = dict(term=term, scenario=scen, delta=0.0 if scen == "N0" else K.M_DELTA[scen], stratum=sname,
                         n_compounds=len(idx), B=B, n_rep=E2.shape[0],
                         median_planted_value=float(np.nanmedian(th_c[idx])), median_planted_value_sd=float(np.nanmedian(sd_c[idx])))
                for lv, (ilo, ihi) in ((90, (1, 2)), (95, (3, 4))):
                    lo, hi = E2[:, ti, si][:, idx, ilo], E2[:, ti, si][:, idx, ihi]
                    ok = np.isfinite(lo) & np.isfinite(hi)
                    cov = (lo <= th_c[idx][None, :]) & (th_c[idx][None, :] <= hi)
                    per_rep = np.array([cov[i][ok[i]].mean() if ok[i].any() else np.nan for i in range(len(cov))])
                    m = mean_rate(per_rep)
                    r.update({f"n_pairs{lv}": int(ok.sum()), f"coverage{lv}": float(cov[ok].mean()), f"coverage{lv}_mc_se": m["mc_se"],
                              f"coverage{lv}_lo95": float(cov[ok].mean()) - 1.959964 * m["mc_se"],
                              f"coverage{lv}_hi95": float(cov[ok].mean()) + 1.959964 * m["mc_se"]})
                    r[f"reading{lv}"] = reading(r[f"coverage{lv}_lo95"], r[f"coverage{lv}_hi95"], lv / 100)
                rows.append(r)
    out("coverage_e2.tsv", pd.DataFrame(rows))

    # ---------------- attenuation: labelled group-1 lines that carry no effect
    rows = []
    for term in K.TERMS:
        for e in ("e3l3", "e3l5"):
            ref = theta[(term, "M20", e)][0]
            for scen in ("M20", "A1", "A2"):
                g = runs[(runs.term == term) & (runs.scenario == scen)]
                v = g[f"{e}_difference_declared"].dropna().astype(bool)
                th = theta[(term, scen, e)]
                rows.append(dict(term=term, estimand=E3_NAME[e], scenario=scen,
                                 unaffected_lines=0 if scen == "M20" else K.A_U[scen], n_g1=int(g.n_g1.iloc[0]),
                                 planted_value=th[0], planted_value_mc_se=th[2], ratio_to_M20=th[0] / ref,
                                 **{f"detect_{k}": v for k, v in rate(v.sum(), len(v)).items()}))
    out("attenuation_e3.tsv", pd.DataFrame(rows))

    # ---------------- E3 detection interpolated at a planted value of 20% (linear in the planted value)
    rows = []
    for term in K.TERMS:
        for e, lab in (("e3l3", "E3 level 3: 95% CI excludes 0"), ("e3l5", "E3 level 5: 95% CI excludes 0")):
            g = ps[(ps.term == term) & (ps.measure == lab) & ps.scenario.isin(["N0"] + list(K.M_DELTA))].copy()
            g = g.sort_values("planted_value")
            x, y = g.planted_value.values, g.rate.values
            i = int(np.searchsorted(x, 20.0))
            if 0 < i < len(x):
                f = (20.0 - x[i - 1]) / (x[i] - x[i - 1])
                rows.append(dict(term=term, measure=lab, target_planted_value=20.0, scenario_lo=g.scenario.values[i - 1],
                                 scenario_hi=g.scenario.values[i], planted_lo=x[i - 1], planted_hi=x[i], rate_lo=y[i - 1],
                                 rate_hi=y[i], interpolated_rate=y[i - 1] + f * (y[i] - y[i - 1])))
    out("e3_power_interpolated.tsv", pd.DataFrame(rows))

    # ---------------- Figure S3 plotted values
    rows = []
    def add(panel, series, term, scen, x, xlab, r):
        th = theta.get((term, scen, "e3l5" if "level 5" in series else "e3l3"), (np.nan,))[0] if panel.startswith("A") else np.nan
        rows.append(dict(figure="S3", panel=panel, series=series, term=term, scenario=scen, x=x, x_label=xlab,
                         planted_e3_pct=th, y=r["rate"], y_lo95=r["lo95"], y_hi95=r["hi95"], y_mc_se=r["mc_se"], n_rep=R, B=B))
    for term in K.TERMS:
        for scen in ["N0"] + list(K.M_DELTA):
            for e, ser in (("e3l3", "E3 on level 3: 95% CI excludes 0"), ("e3l5", "E3 on level 5: 95% CI excludes 0")):
                r = ps[(ps.term == term) & (ps.scenario == scen) & (ps.measure == ser.replace("E3 on", "E3"))].iloc[0]
                add("A: pooled magnitude", ser, term, scen, 0.0 if scen == "N0" else K.M_DELTA[scen],
                    "planted pooled response-size difference (fraction)", r)
            r = ps[(ps.term == term) & (ps.scenario == scen) & (ps.measure == "panel-level test: p < 0.025")].iloc[0]
            add("A: pooled magnitude", "panel-level test: p < 0.025", term, scen, 0.0 if scen == "N0" else K.M_DELTA[scen],
                "planted pooled response-size difference (fraction)", r)
        for pan, keys in (("B: 50% larger response in k compounds", K.C_K), ("C: +1 SD in a 49-gene module in k compounds", K.G_K)):
            for scen, k in keys.items():
                for ser in ("planted compounds at q < 0.1 (share of planted)", "panel-level test: p < 0.025",
                            "at least one planted compound at q < 0.1"):
                    r = ps[(ps.term == term) & (ps.scenario == scen) & (ps.measure == ser)].iloc[0]
                    add(pan, ser, term, scen, k, "number of compounds with the planted effect (k)", r)
    out("figS3_power_values.tsv", pd.DataFrame(rows))
    print(f"summaries written: R {R}, B {B}")


if __name__ == "__main__":
    main()
