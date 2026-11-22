"""Analysis 3 (../PLAN.md): TP53 status as a yardstick on the same tests and scales.

From the run's outputs: the TP53 comparison exposure (results/G/G_controls.tsv) and, from the capture rerun's saved
table and permutation statistics of that call, AMG-232's p, BH q and rank in the factorial line set
(gcore.finish_table, as ganalyses.exposure_run applies it). The recomputed panel p must equal G_controls.tsv.
Computed here with the frozen functions: E3 for TP53 status (mutant minus wild type) on level 3 and level 5, on the
factorial set's TP53 analysis set (as exposure_run builds it, rule 7i A) and on the PC1 gate's line set (TP53_cmp,
task D's rule); and AMG-232's E2 under TP53 status (level 5). No state label is read except the factorial line-set
membership, which the frozen design builds from the line table.
Writes results/tp53_yardstick.tsv and results/tp53_compounds_q_lt_01.tsv.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as K  # noqa: E402

import numpy as np
import pandas as pd

C, GS, GA = K.C, K.GS, K.GA
CAP_TAB = "0019_exposure__TP53_status_comparison_exposure__P1_run_term_table.tsv"
CAP_NPZ = "0019_exposure__TP53_status_comparison_exposure__P1_run_term.npz"
CAP_NULL = "0020_exposure__TP53_status_comparison_exposure__P1_panel_null.npz"


def main():
    D = K.Data(with_level3=True)
    P, cfg = D.P, D.cfg
    rows = []
    # ---------------- the run's outputs
    ctl = pd.read_csv(os.path.join(K.GRES, "G_controls.tsv"), sep="\t", float_precision="round_trip")
    for pn in ("P1", "P1b", "P2"):
        r = ctl[(ctl.panel == pn) & (ctl.exposure == "TP53 status (comparison exposure)")].iloc[0]
        rows.append(dict(quantity="panel-level test, TP53 status as exposure", panel=pn, line_set="factorial set (TP53 analysis set)",
                         scale="level5", n_lines=int(r.n_lines), n_g1=int(r.n_g1), n_g0=int(r.n_g0), n_lineages=int(r.n_lineages),
                         value=float(r.panel_stat), p=float(r.panel_p), n_tested=int(r.n_tested), n_p_lt_005=int(r.n_p_lt_005),
                         n_q_lt_01=int(r.n_q_lt_01), source="results/G/G_controls.tsv"))
    gates = pd.read_csv(os.path.join(K.GRES, "G_gates.tsv"), sep="\t", float_precision="round_trip")
    pc1 = gates[gates.gate == "PC1"].iloc[0]
    rows.append(dict(quantity="PC1 gate: AMG-232 per-compound test", panel="P1", line_set="PC1 gate set (TP53_cmp, task D's rule)",
                     scale="level5", n_lines=int(pc1.n_lines), n_g1=int(pc1.n_g1), n_g0=int(pc1.n_g0), n_lineages=int(pc1.n_lineages),
                     value=float(pc1.T_obs), p=float(pc1.p), rank=int(pc1["rank"]), n_ranked=int(pc1.n_ranked),
                     source="results/G/G_gates.tsv"))
    rows.append(dict(quantity="PC1 gate: panel-level test", panel="P1", line_set="PC1 gate set (TP53_cmp, task D's rule)",
                     scale="level5", n_lines=int(pc1.n_lines), value=float(pc1.panel_stat), p=float(pc1.panel_p),
                     source="results/G/G_gates.tsv"))
    # ---------------- AMG-232 in the factorial set, from the capture
    tab = pd.read_csv(os.path.join(K.CAPTURE, CAP_TAB), sep="\t", float_precision="round_trip")
    z = np.load(os.path.join(K.CAPTURE, CAP_NPZ))
    assert list(z["compounds"]) == list(tab.compound) == list(P.comps)
    permT = z["permT"].astype(np.float64)
    dfF, obs, pp, null = C.finish_table(tab, permT, compounds_only=P.comps)
    zn = np.load(os.path.join(K.CAPTURE, CAP_NULL))
    r0 = ctl[(ctl.panel == "P1") & (ctl.exposure == "TP53 status (comparison exposure)")].iloc[0]
    same = bool(pp == float(zn["panel_p"]) == float(r0.panel_p) and abs(obs - float(r0.panel_stat)) <= 1e-12 * abs(obs))
    assert same, (pp, float(zn["panel_p"]), float(r0.panel_p))
    a = dfF[dfF.compound == "AMG-232"].iloc[0]
    rows.append(dict(quantity="AMG-232 per-compound test, TP53 status as exposure", panel="P1",
                     line_set="factorial set (TP53 analysis set)", scale="level5", n_lines=int(a.n_lines), n_g1=int(a.n_g1),
                     n_g0=int(a.n_g0), n_lineages=int(a.n_lineages), value=float(a.T_obs), p=float(a.p), q=float(a.q),
                     rank=int(a["rank"]), n_ranked=int(a.n_ranked),
                     source=f"capture/{CAP_TAB} and {CAP_NPZ}; gcore.finish_table; panel p reproduces G_controls.tsv: {same}"))
    qq = dfF[dfF.in_ranking & (dfF.q < cfg.const("H2_bh_q"))].sort_values("rank")[["compound", "T_obs", "p", "q", "rank", "n_lines"]]
    qq.to_csv(os.path.join(K.RES, "tp53_compounds_q_lt_01.tsv"), sep="\t", index=False)
    # ---------------- E3 and AMG-232's E2 for TP53 status (computed here)
    lt = pd.read_csv(D.I.real_line_table_path(), sep="\t")
    al = GS.Aligned(P, lt)
    prim = GS.primary_settings(cfg)
    dY = al.design("factY", prim)
    sl = pd.read_csv(os.path.join(K.GRES, "G_sets_lines.tsv"), sep="\t")
    f41 = set(sl[(sl.panel == "P1") & (sl.family == "factorial")].cell_iname)
    assert set(np.array(P.lines)[dY["mask"]]) == f41
    tcode = GS.tp53_code(cfg)
    tp1 = GS.tp53_columns(P.cov, "7gB" if tcode == "7gB" else "cmp")[0]
    m_fact = C.analysis_set(dY["lin"], tp1, [0], dY["mask"])
    v = np.array(P.cov[cfg.const("pc1")["genotype_column"]].astype(str).tolist(), dtype=object)
    x = np.where(v == "mutant", 1.0, np.where(v == "wild_type", 0.0, np.nan))[:, None]
    m_pc1 = C.analysis_set(P.lin, x, [0], np.ones(len(P.lines), bool))
    s5, s3 = D.s0()
    S = {"level3": np.nanmean(s3, axis=1), "level5": np.nanmean(s5, axis=1)}
    sets = (("factorial set (TP53 analysis set)", dict(term="TP53", lin=dY["lin"], Xo=tp1, k=0, mask=m_fact, rule=prim["rule"])),
            ("PC1 gate set (TP53_cmp, task D's rule)", dict(term="TP53", lin=P.lin, Xo=x, k=0, mask=m_pc1, rule=None)))
    for name, d in sets:
        for scale, tn in (("level3", "R3_P1_d111"), ("level5", "P1_resp")):
            e = GA.e3_row(D.ctx, "P1", "TP53", d, S[scale], scale, tn, "yardstick")
            t = d["Xo"][d["mask"], 0]
            rows.append(dict(quantity="E3, TP53 status (mutant minus wild type), % of mean S", panel="P1", line_set=name, scale=scale,
                             n_lines=int(e["n_lines"]), n_g1=int(np.sum(t == 1)), n_g0=int(np.sum(t == 0)),
                             n_lineages=len(set(np.asarray(d["lin"])[d["mask"]])), df=int(e["df"]), mean_S=e["mean_S"],
                             value=e["estimate"], lo90=e["lo90"], hi90=e["hi90"], lo95=e["lo95"], hi95=e["hi95"],
                             difference_declared=e["difference_declared"], source="computed here: ganalyses.e3_row"))
    d = sets[0][1]
    j = P.comps.index("AMG-232")
    f = C.scalar_fit2(d["lin"], d["Xo"], 0, np.flatnonzero(d["mask"]), s5[:, j], rule=d["rule"])
    est, lo90, hi90 = GA.express(f, cfg, 0.90)
    _, lo95, hi95 = GA.express(f, cfg, 0.95)
    rows.append(dict(quantity="E2 of AMG-232, TP53 status (mutant minus wild type), % of the compound's mean s", panel="P1",
                     line_set=sets[0][0], scale="level5", n_lines=int(f["n"]), df=int(f["df"]), mean_S=f["mean_s"], value=est,
                     lo90=lo90, hi90=hi90, lo95=lo95, hi95=hi95, difference_declared=bool(lo95 > 0 or hi95 < 0),
                     source="computed here: gcore.scalar_fit2 + ganalyses.express"))
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(K.RES, "tp53_yardstick.tsv"), sep="\t", index=False)
    print(out[["quantity", "line_set", "scale", "n_lines", "value", "p", "q", "rank", "lo95", "hi95"]].to_string(index=False))
    print(qq.to_string(index=False))


if __name__ == "__main__":
    main()
