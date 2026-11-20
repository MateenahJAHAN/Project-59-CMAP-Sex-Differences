"""Merge steps of task G: the only code that joins a state or label to a response.

It is called by run_g.py only after the guard and the gates. Functions here never print a statistic; they
return tables that run_g.py writes under results/G/.

Plan items (draft section in brackets): H1 [4, 6], H2 rules 1 to 4 [4], E1, E2, E4, E5 [4], H3 / E3 [4, 6],
H4 [3, 4], P2 agreement (7a) [6], leads L1 to L4 (7l) [1, 4], NC2 to NC5 and comparison exposures [5c].
"""
import numpy as np, pandas as pd
from scipy import stats
import gcore as C
import engine as E

DMSO = "DMSO (pseudo-compound)"


# ---------------------------------------------------------------- helpers
def e_unit(cfg):
    return "percent of mean S" if cfg.decide("7j")["unit"] == "A" else "100 x coefficient"


def e_interval(cfg):
    return "hc3" if cfg.decide("7j")["interval"] == "A" else "classical"


def express(f, cfg, level=0.90):
    """Estimate and CI of a scalar fit in the 7j unit and interval."""
    if f is None:
        return np.nan, np.nan, np.nan
    scale = 100.0 / f["mean_s"] if cfg.decide("7j")["unit"] == "A" else 100.0
    hw = f[f"hw{int(round(level*100))}_{e_interval(cfg)}"]
    return scale * f["beta"], scale * (f["beta"] - hw), scale * (f["beta"] + hw)


def joint_null_count(permT, cols, alpha, q=0.95):
    """7k B: per permutation, the number of compounds whose permuted statistic has permutation p < alpha."""
    X = permT[:, cols]
    B = X.shape[0]
    cnt = np.zeros(B, int)
    for j in range(X.shape[1]):
        r = stats.rankdata(-X[:, j], method="max")  # number of permuted values >= this one (itself included)
        cnt += (r / (B + 1) < alpha).astype(int)
    return float(np.quantile(cnt, q)), cnt


def beta_with_cov(lin, Xo, k, rows, extra_col, rule):
    """Design with one extra covariate placed first (the tested column moves to k + 1), on the rows where the
    covariate is known. Returns (Design or None, rows, reason)."""
    Xc = np.c_[extra_col, Xo]
    rows = rows[np.isfinite(extra_col[rows])]
    rows, d, why = C.compound_design(lin, Xc, k + 1, rows, rule["min_group"], rule["min_relab"])
    return d, rows, why


# ---------------------------------------------------------------- one term (H1, H2, E1, E2, E4, NC3, NC4)
def run_main_term(ctx, pn, term, d, B, keep_perm=False, label="primary"):
    cfg, P = ctx.cfg, ctx.P[pn]
    s5 = ctx.s_level5[pn]
    head = dict(panel=pn, term=term, analysis=label, n_lines=d["n_lines"], n_g1=d["n_g1"], n_g0=d["n_g0"],
                n_lineages=d["n_lineages"], n_relabellings=float(d["n_relabellings"]), contrast_ok=d["contrast_ok"],
                contrast_reason=d["contrast_reason"], B=B)
    if not d["contrast_ok"]:
        return dict(panel_row=dict(head, n_tested=0), per_cpd=None, permT=None, recs=None)
    extra = {DMSO: P.dmso} if (P.dmso is not None and label == "primary") else None
    df, permT, recs = C.run_term(P, d["Xo"], d["k"], d["mask"], B, extra=extra, lin=d["lin"], rule=d["rule"])
    if "p" not in df.columns or not df.tested.any():
        return dict(panel_row=dict(head, n_tested=0, untested="; ".join(f"{a}: {b}" for a, b in df.reason.value_counts().items())),
                    per_cpd=None, permT=None, recs=None)
    df, obs_all, p_all, null_all = C.finish_table(df, permT, compounds_only=P.comps)
    mask_lines = np.flatnonzero(d["mask"])
    complete = np.array([bool(r.get("_j", -1) >= 0 and P.avail[mask_lines, r["_j"]].all()) if r["tested"] else
                         (c in P.comps and P.avail[mask_lines, P.comps.index(c)].all()) for r, c in zip(recs, df.compound)])
    df["complete_over_analysis_set"] = complete
    act = ctx.activity[pn]
    df["above_median_activity"] = df.compound.map(act).fillna(False).astype(bool)
    in_rank = df.in_ranking.values
    # H1 compound set (7h)
    h1_complete = cfg.decide("7h") == "A"
    sub_complete = in_rank & complete
    o_c, p_c, null_c, n_c = C.panel_stat(df, permT, sub_complete)
    if h1_complete:
        H1 = dict(H1_compound_set="complete compounds (7h A)", panel_stat=o_c, panel_p=p_c, n_in_H1=n_c, null=null_c)
        SENS = dict(sens_compound_set="all tested compounds", sens_panel_stat=obs_all, sens_panel_p=p_all, sens_n=int(in_rank.sum()))
    else:
        H1 = dict(H1_compound_set="all tested compounds (7h B)", panel_stat=obs_all, panel_p=p_all, n_in_H1=int(in_rank.sum()), null=null_all)
        SENS = dict(sens_compound_set="complete compounds", sens_panel_stat=o_c, sens_panel_p=p_c, sens_n=n_c)
    # activity subset (label-free activity, D_activity_scores.tsv)
    sub_act = in_rank & df.above_median_activity.values
    o_a, p_a, _, n_a = C.panel_stat(df, permT, sub_act)
    df["q_activity_subset"] = np.nan
    df.loc[sub_act, "q_activity_subset"] = C.bh(df.loc[sub_act, "p"].values)
    # E1 and E2 per compound
    e1, e1hw, e1d, e2 = [], [], [], []
    for r in recs:
        if not r["tested"] or r["compound"] == DMSO:
            e1.append(np.nan); e1hw.append(np.nan); e1d.append(np.nan); e2.append((np.nan,) * 4); continue
        dd = r["_d"]
        tq = stats.t.ppf(0.95, dd.df)
        with np.errstate(invalid="ignore", divide="ignore"):
            z = np.where(r["_sig"] > 0, r["_b"] / r["_sig"], 0.0)
            hw = np.where(r["_sig"] > 0, tq * r["_se"] / r["_sig"], np.nan)
        e1.append(float(np.sqrt(np.mean(z ** 2)))); e1hw.append(float(np.nanmedian(hw))); e1d.append(float(tq * np.sqrt(dd.xtxi_kk)))
        f = C.scalar_fit2(d["lin"], d["Xo"], d["k"], r["_rows"], s5[:, r["_j"]], rule=d["rule"])
        est, lo, hi = express(f, cfg)
        e2.append((est, lo, hi, f["n"] if f else np.nan))
    df["E1_rms_sd"] = e1; df["E1_hw90_hc3_sd_median"] = e1hw; df["E1_hw90_design_sd"] = e1d
    df["E2_est"] = [x[0] for x in e2]; df["E2_lo90"] = [x[1] for x in e2]; df["E2_hi90"] = [x[2] for x in e2]
    df["E2_n"] = [x[3] for x in e2]
    df["E2_unit"] = e_unit(cfg); df["E2_interval"] = e_interval(cfg)
    # E4 (H2 rule 3)
    t = df[in_rank]
    n_p = int((t.p < cfg.const("E4_count_p")).sum())
    n_q = int((t.q < cfg.const("H2_bh_q")).sum())
    if cfg.decide("7k") == "A":
        thr, thr_src = float(cfg.const("E4_threshold_7kA")), "fixed: task D pseudo-label 95th percentile (7k A)"
    else:
        thr, _ = joint_null_count(permT, np.flatnonzero(in_rank), cfg.const("E4_count_p"))
        thr_src = "joint permutation null at G (7k B)"
    rule3 = bool(n_p > thr)
    # H2 rules 1 and 2 (NC4 is the rank)
    cut = int(np.floor(cfg.const("H2_top_fraction") * int(in_rank.sum())))
    df["top5pct_cut"] = cut
    df["rule1_q"] = df.q < cfg.const("H2_bh_q")
    df["rule2_top5pct"] = df["rank"] <= cut
    df["rule3_E4"] = rule3
    # NC3 (DMSO pseudo-compound under the state term)
    dm = df.compound == DMSO
    dmso = dict(dmso_T=np.nan, dmso_p=np.nan, dmso_rank_if_inserted=np.nan)
    if dm.any() and df.loc[dm, "tested"].iloc[0]:
        dmso = dict(dmso_T=float(df.loc[dm, "T_obs"].iloc[0]), dmso_p=float(df.loc[dm, "p"].iloc[0]),
                    dmso_rank_if_inserted=float(1 + (df.loc[in_rank, "T_obs"] > df.loc[dm, "T_obs"].iloc[0]).sum()))
        df.loc[dm, "rank"] = dmso["dmso_rank_if_inserted"]
    null = H1.pop("null")
    prow = dict(head, n_tested=int(in_rank.sum()), **H1, null_median=float(np.median(null)) if null is not None else np.nan,
                E5_excess=float(H1["panel_stat"] - np.median(null)) if null is not None else np.nan,
                **SENS, activity_subset_panel_stat=o_a, activity_subset_panel_p=p_a, activity_subset_n=n_a,
                n_p_lt_005=n_p, n_q_lt_01=n_q, E4_threshold=thr, E4_threshold_source=thr_src, E4_met=rule3, **dmso)
    df.insert(0, "panel", pn); df.insert(1, "term", term); df.insert(2, "analysis", label)
    out = dict(panel_row=prow, per_cpd=df, recs=recs, design=d, permT=permT if keep_perm else None)
    return out


# ---------------------------------------------------------------- H2 rule 4 (ploidy and PRISM viability)
def rule4(ctx, pn, res):
    """Ratio = <b_adj, b_primary> / |b_primary|^2 for the gene-coefficient vector of the tested term, after adding
    ploidy, and after adding the compound's PRISM AUC, as a covariate. Rule 4 holds if both ratios are at least
    plan constant H2_rule4_min_retained (same sign and at least half the estimate). Not testable counts as not met."""
    cfg, P, d = ctx.cfg, ctx.P[pn], res["design"]
    ploidy = P.cov.ploidy.to_numpy(dtype=float)
    auc = ctx.auc.get(pn)
    out = []
    for r in res["recs"]:
        if not r["tested"] or r["compound"] == DMSO:
            continue
        b0 = r["_b"]; nb = float(b0 @ b0)
        row = dict(panel=pn, term=d["term"], compound=r["compound"])
        for nm, col in (("ploidy", ploidy), ("prism_auc", auc[:, r["_j"]] if auc is not None else None)):
            if col is None:
                row.update({f"{nm}_ratio": np.nan, f"{nm}_ratio_same_rows": np.nan, f"{nm}_n": 0, f"{nm}_note": "no PRISM data in this run"})
                continue
            dd, rows, why = beta_with_cov(d["lin"], d["Xo"], d["k"], r["_rows"], col, d["rule"])
            if dd is None:
                row.update({f"{nm}_ratio": np.nan, f"{nm}_ratio_same_rows": np.nan, f"{nm}_n": len(rows),
                            f"{nm}_note": f"not testable: {why}" if np.isfinite(col[r["_rows"]]).any() else "no covariate value"})
                continue
            Y = P.R[rows, r["_j"], :]
            ba, _, _, _ = dd.observed(Y)
            rows2, d0, _ = C.compound_design(d["lin"], d["Xo"], d["k"], rows, d["rule"]["min_group"], d["rule"]["min_relab"])
            b_same = d0.observed(Y)[0] if d0 is not None and len(rows2) == len(rows) else None
            row.update({f"{nm}_ratio": float(ba @ b0 / nb) if nb > 0 else np.nan,
                        f"{nm}_ratio_same_rows": float(ba @ b_same / (b_same @ b_same)) if b_same is not None and b_same @ b_same > 0 else np.nan,
                        f"{nm}_n": len(rows), f"{nm}_note": ""})
        thr = cfg.const("H2_rule4_min_retained")
        row["rule4_met"] = bool(np.isfinite(row["ploidy_ratio"]) and np.isfinite(row["prism_auc_ratio"])
                                and row["ploidy_ratio"] >= thr and row["prism_auc_ratio"] >= thr)
        out.append(row)
    return pd.DataFrame(out)


def apply_findings(ctx, pn, df, r4, old):
    """H2: a compound is a finding only if rules 1 to 4 hold; old-list compounds cannot be findings from P1 alone."""
    df = df.merge(r4, on=["panel", "term", "compound"], how="left", validate="1:1")
    df["rule4_met"] = df.rule4_met.fillna(False).astype(bool)
    df["old_list"] = df.compound.isin(old)
    allr = df.rule1_q.fillna(False) & df.rule2_top5pct.fillna(False) & df.rule3_E4.fillna(False) & df.rule4_met
    df["meets_rules_1_to_4"] = allr & df.in_ranking
    df["finding"] = df.meets_rules_1_to_4 & ~df.old_list & (pn == ctx.cfg.const("primary_panel"))
    df["finding_note"] = np.where(df.finding, "finding (candidate for follow-up, not a claim about patients)",
                         np.where(df.meets_rules_1_to_4 & df.old_list, "meets rules 1 to 4; on the old 50-drug list, so not a finding from P1 alone",
                         np.where(df.meets_rules_1_to_4, "meets rules 1 to 4 outside the primary panel (no claim)", "")))
    return df


# ---------------------------------------------------------------- gene level for findings
def gene_level(ctx, pn, res, compounds, B):
    P, d = ctx.P[pn], res["design"]
    js = [P.comps.index(c) for c in compounds]
    recs = C.run_term_signed(P, d["Xo"], d["k"], d["mask"], B, lambda nm, rec: (lambda t, b: np.abs(t)), comps=js,
                             lin=d["lin"], rule=d["rule"])
    rows = []
    for r in recs:
        if not r["tested"]:
            continue
        obs = np.abs(r["_t"])
        p = (1 + (r["perm"] >= obs[None, :] * (1 - C.TIE_TOL)).sum(0)) / (B + 1)
        q = C.bh(p)
        for g, gene in enumerate(P.genes):
            rows.append(dict(panel=pn, term=d["term"], compound=r["compound"], gene=gene, beta=float(r["_b"][g]),
                             se_hc3=float(r["_se"][g]), t=float(r["_t"][g]), p_perm=float(p[g]), q_bh_genes=float(q[g])))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- H3 (E3) and S12
def e3_row(ctx, pn, term, d, S, scale, tensor, role):
    cfg = ctx.cfg
    f = C.scalar_fit2(d["lin"], d["Xo"], d["k"], np.flatnonzero(d["mask"]), S, rule=d["rule"], levels=cfg.const("ci_levels"))
    row = dict(panel=pn, term=term, scale=scale, tensor=tensor, role=role, unit=e_unit(cfg), interval=e_interval(cfg))
    if f is None:
        row.update(n_lines=0, note="not estimable")
        return row
    est, lo90, hi90 = express(f, cfg, 0.90)
    _, lo95, hi95 = express(f, cfg, 0.95)
    row.update(n_lines=f["n"], df=f["df"], mean_S=f["mean_s"], estimate=est, lo90=lo90, hi90=hi90, lo95=lo95, hi95=hi95,
               difference_declared=bool(lo95 > 0 or hi95 < 0))
    return row


def equivalence(row, p_labelfree, cfg):
    m, need = cfg.const("E3_margin_pct"), cfg.const("E3_equivalence_min_probability")
    if not np.isfinite(row.get("estimate", np.nan)):
        return "not estimable"
    if not np.isfinite(p_labelfree):
        return "no label-free probability for this setting: no equivalence claim; the 90% CI is a bound"
    if p_labelfree < need:
        return f"label-free P(90% CI inside +-{m:g}%) = {p_labelfree:.3f} < {need}: no equivalence claim; the 90% CI is a bound"
    inside = row["lo90"] > -m and row["hi90"] < m
    return "equivalence declared (90% CI inside the margin)" if inside else "equivalence not shown (90% CI not inside the margin)"


# ---------------------------------------------------------------- H4 attenuation (labels only)
def h4_attenuation(al, d_ann):
    """Expected dilution of a pure Y (or X) effect under annotated labels, from the line table only:
    P(Y=1 | M) - P(Y=1 | F) and P(X=1 | F) - P(X=1 | M) among annotated-set lines with a factorial state."""
    m = d_ann["mask"]
    M = d_ann["Xo"][:, d_ann["k"]]
    Y = al.a.fact_Y.to_numpy(dtype=float); X = al.a.fact_X.to_numpy(dtype=float)
    ok = m & np.isfinite(Y) & np.isfinite(X)
    f = lambda v, g: float(np.mean(v[ok & (M == g)])) if (ok & (M == g)).any() else np.nan
    return dict(n_annotated_set=int(m.sum()), n_with_state=int(ok.sum()), P_Y1_given_M=f(Y, 1), P_Y1_given_F=f(Y, 0),
                P_X1_given_F=f(X, 0), P_X1_given_M=f(X, 1), dilution_Y=f(Y, 1) - f(Y, 0), dilution_X=f(X, 0) - f(X, 1),
                n_M_Y_lost=int(np.sum(ok & (M == 1) & (Y == 0))), n_M=int(np.sum(ok & (M == 1))))


# ---------------------------------------------------------------- P2 agreement (7a)
def p2_agreement(ctx, res1, res2, B):
    """For compounds tested in both P1 and P2: Spearman correlation over genes of the tested coefficient vectors,
    with a one-sided permutation p (P2 permutations), and the sign agreement of E2."""
    P2, d2 = ctx.P["P2"], res2["design"]
    b1 = {r["compound"]: r["_b"] for r in res1["recs"] if r["tested"] and r["compound"] != DMSO}
    shared = [r["compound"] for r in res2["recs"] if r["tested"] and r["compound"] in b1]
    if not shared:
        return pd.DataFrame()
    js = [P2.comps.index(c) for c in shared]

    def red_for(name, rec):
        rk = stats.rankdata(b1[name]); rk = (rk - rk.mean()) / np.sqrt(((rk - rk.mean()) ** 2).sum())

        def red(t, b):
            R = stats.rankdata(b, axis=1)
            R = R - R.mean(axis=1, keepdims=True)
            R /= np.sqrt((R ** 2).sum(axis=1, keepdims=True))
            return R @ rk
        return red
    recs = C.run_term_signed(P2, d2["Xo"], d2["k"], d2["mask"], B, red_for, comps=js, lin=d2["lin"], rule=d2["rule"])
    e2_1 = res1["per_cpd"].set_index("compound").E2_est
    e2_2 = res2["per_cpd"].set_index("compound").E2_est
    rows = []
    for r in recs:
        if not r["tested"]:
            continue
        rho = float(r["obs"]); null = r["perm"]
        p = float((1 + np.sum(null >= rho - 1e-12)) / (B + 1))
        same = bool(np.sign(e2_1[r["compound"]]) == np.sign(e2_2[r["compound"]])) if np.isfinite(e2_1[r["compound"]]) and np.isfinite(e2_2[r["compound"]]) else False
        rows.append(dict(term=d2["term"], compound=r["compound"], spearman_P1_P2=rho, p_perm_one_sided=p,
                         E2_P1=float(e2_1[r["compound"]]), E2_P2=float(e2_2[r["compound"]]), E2_same_sign=same,
                         shape_and_E2_agree=bool(p < ctx.cfg.const("replication_p") and rho > 0 and same)))
    return pd.DataFrame(rows)


def p2_labels(cfg, agree):
    """Wording under 7a for panel- and compound-level agreement."""
    opt = cfg.decide("7a")
    if opt == "A":
        return "replicated" if agree else "not replicated in P2 (a P2 null does not count against P1)"
    if opt == "B":
        return "P2 reported as description only (7a B)"
    return "also seen in P2" if agree else "not seen in P2 (a P2 null does not count against P1)"


# ---------------------------------------------------------------- leads (7l)
def lead_test(ctx, pn, lead, d, B, scale_s, use_shape=False, h2_table=None):
    P = ctx.P[pn]
    c = lead["compound"]
    row = dict(lead=lead["id"], compound=c, panel=pn, term=lead["term"], alternative=lead["alternative"])
    if c not in P.comps:
        row.update(note="compound not in panel"); return row
    j = P.comps.index(c)
    if use_shape:  # 7l C for L4: the H2 shape test
        t = h2_table[(h2_table.compound == c)]
        row.update(statistic="T_c (shape)", p=float(t.p.iloc[0]) if len(t) and t.tested.iloc[0] else np.nan, scale="level 5")
        return row
    s = scale_s[:, j]
    Yov = {j: np.where(np.isfinite(s), s, np.nan)[:, None]}
    recs = C.run_term_signed(P, d["Xo"], d["k"], d["mask"], B, lambda nm, rec: (lambda t, b: t[:, 0]), comps=[j],
                             Y_override=Yov, lin=d["lin"], rule=d["rule"])
    r = recs[0]
    if not r["tested"]:
        row.update(note=f"not tested: {r['reason']}"); return row
    t_obs, null = float(r["obs"]), r["perm"]
    tol = 1e-9 * max(1.0, abs(t_obs))
    if lead["alternative"] == "greater":
        p = (1 + np.sum(null >= t_obs - tol)) / (B + 1)
    elif lead["alternative"] == "less":
        p = (1 + np.sum(null <= t_obs + tol)) / (B + 1)
    else:
        p = (1 + np.sum(np.abs(null) >= abs(t_obs) - tol)) / (B + 1)
    f = C.scalar_fit2(d["lin"], d["Xo"], d["k"], r["_rows"], s, rule=d["rule"])
    est, lo, hi = express(f, ctx.cfg)
    row.update(statistic="E2 HC3 t, Freedman-Lane", t_obs=t_obs, p=float(p), E2_est=est, E2_lo90=lo, E2_hi90=hi,
               n_lines=r["n_lines"], note="")
    return row


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    m = len(p)
    adj = np.empty(m)
    run = 0.0
    for i, idx in enumerate(o):
        run = max(run, min(1.0, (m - i) * p[idx]))
        adj[idx] = run
    return adj


# ---------------------------------------------------------------- controls and comparison exposures
def exposure_run(ctx, pn, name, Xo, k, tcols, mask, lin, rule, B):
    P = ctx.P[pn]
    m = C.analysis_set(lin, Xo, tcols, mask)
    t = Xo[m, k]
    base = dict(panel=pn, exposure=name, n_lines=int(m.sum()), n_g1=int(np.sum(t == 1)), n_g0=int(np.sum(t == 0)),
                n_lineages=len(set(lin[m])), B=B)
    if m.sum() < 4 or np.unique(t).size < 2:
        base.update(n_tested=0, note="exposure constant or too few lines: not testable"); return base
    df, permT, recs = C.run_term(P, Xo, k, m, B, lin=lin, rule=rule)
    if "p" not in df.columns or not df.tested.any():
        base.update(n_tested=0, note="no compound testable"); return base
    df, obs, pp, null = C.finish_table(df, permT, compounds_only=P.comps)
    base.update(n_tested=int(df.in_ranking.sum()), panel_stat=obs, panel_p=pp, null_median=float(np.median(null)),
                E5_excess=float(obs - np.median(null)), n_p_lt_005=int((df[df.in_ranking].p < 0.05).sum()),
                n_q_lt_01=int((df[df.in_ranking].q < 0.1).sum()), note="")
    return base


def ploidy_split(P, mask, lin):
    v = P.cov.ploidy.to_numpy(dtype=float)
    hi = np.full(len(v), np.nan)
    for l in np.unique(lin):
        idx = np.flatnonzero((lin == l) & mask & np.isfinite(v))
        if len(idx) >= 2:
            hi[idx] = (v[idx] > np.median(v[idx])).astype(float)
    return hi
