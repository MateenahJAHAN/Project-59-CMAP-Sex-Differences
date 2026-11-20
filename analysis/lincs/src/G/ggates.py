"""Gates of task G (draft section 5), run before any state or label is joined to a response.

Everything here is label-free: genotype labels (TP53) from task D's state-free covariate table, pseudo-labels
from per-lineage counts, and line sets given as membership lists without any group column.

- PC1 (level 5, task D's control): AMG-232 under lineage + TP53 (CMP). Must reproduce task D's
  D_pc1_summary.tsv row exactly (rank, p, T_obs), and pass the rule (p < 0.05, top 5%).
- K4 (level 3, task E's engine): PC1 on the level-3 tensor of 7d, plus the CDKN1A sign. Must reproduce task D2's
  D2_K4_summary.tsv (7d A) or task E's E_pc1_summary.tsv (7d B).
- Calibration (K6): pseudo-labels within lineage; mean share of compounds with p < 0.05 within [0.03, 0.07] for
  every term analysed. Prior runs (tasks D and D2) are replayed and must be reproduced exactly.
- Label-free E3 precision: task D2's draw scheme; gives P(90% CI inside +-20% under no effect) for H3.
- DMSO gate: only under 7b C.
"""
import os, time, math
import numpy as np, pandas as pd
from scipy import stats
from concurrent.futures import ProcessPoolExecutor
import gcore as C
import engine_E

CODES_FACT = {"Y": (1, 4), "X": (2, 4)}


# ======================= PC1 (level 5) =======================
def pc1(P, cfg, B):
    pc = cfg.const("pc1")
    v = np.array(P.cov[pc["genotype_column"]].astype(str).tolist(), dtype=object)
    x = np.where(v == "mutant", 1.0, np.where(v == "wild_type", 0.0, np.nan))
    Xo = x[:, None]
    mask = C.analysis_set(P.lin, Xo, [0], np.ones(len(P.lines), bool))
    extra = {"DMSO (pseudo-compound)": P.dmso} if P.dmso is not None else None
    df, permT, recs = C.run_term(P, Xo, 0, mask, B, extra=extra, dtype=np.float64)
    df, obs, pp, null = C.finish_table(df, permT, compounds_only=P.comps)
    gi = P.genes.index(pc["gene"])
    df["beta_gene"] = [r["_b"][gi] if r["tested"] else np.nan for r in recs]
    dm = df.compound == "DMSO (pseudo-compound)"
    if dm.any() and df.loc[dm, "tested"].iloc[0]:
        df.loc[dm, "rank"] = 1 + (df.loc[df.in_ranking, "T_obs"] > df.loc[dm, "T_obs"].iloc[0]).sum()
    a = df[df.compound == pc["compound"]].iloc[0]
    cut = int(np.floor(cfg.const("pc1_top_fraction") * a.n_ranked))
    row = dict(gate="PC1", scale="level5", tensor=f"{P.pn}_resp.npy", panel=P.pn, compound=pc["compound"],
               n_lines=int(mask.sum()), n_g1=int(np.nansum(x[mask] == 1)), n_g0=int(np.nansum(x[mask] == 0)),
               n_lineages=len(set(P.lin[mask])), T_obs=float(a.T_obs), p=float(a.p), rank=float(a["rank"]),
               n_ranked=int(a.n_ranked), top_cut=cut, beta_gene=float(a.beta_gene), panel_stat=obs, panel_p=pp, B=B,
               passes=bool(a.p < cfg.const("pc1_p") and a["rank"] <= cut))
    if dm.any():
        r = df[dm].iloc[0]
        row.update(dmso_T=float(r.T_obs), dmso_p=float(r.p), dmso_rank_if_inserted=float(r["rank"]))
    return row, df


def pc1_reproduces(row, prior):
    """Exact reproduction of task D's PC1 row (P1, TP53 (CMP))."""
    ok = (row["rank"] == prior["rank"] and row["n_ranked"] == prior["n_ranked"] and math.isclose(row["p"], prior["p"], rel_tol=1e-12)
          and math.isclose(row["T_obs"], prior["T_obs"], rel_tol=1e-9) and math.isclose(row["beta_gene"], prior["beta_CDKN1A"], rel_tol=1e-9)
          and math.isclose(row["panel_p"], prior["panel_p"], rel_tol=1e-12))
    if "dmso_p" in row and np.isfinite(prior.get("dmso_p", np.nan)):
        ok = ok and math.isclose(row["dmso_p"], prior["dmso_p"], rel_tol=1e-12)
    return bool(ok)


# ======================= K4 (level 3, task E's engine) =======================
def k4(R3, Z5, lines_E, comps_E, genes, pc1_lines, cfg, B, tensor_name):
    pos = {l: i for i, l in enumerate(lines_E)}
    assert set(pc1_lines.columns) == {"cell_iname", "cell_lineage", "tp53_mutant"} and pc1_lines.cell_iname.is_unique
    idx = np.array([pos[l] for l in pc1_lines.cell_iname])
    lin = pc1_lines.cell_lineage.values
    x = pc1_lines.tp53_mutant.astype(float).values
    rng = np.random.default_rng(cfg.run("k4_seed"))
    perms = engine_E.within_lineage_perms(lin, B, rng)
    g = genes.index(cfg.const("pc1")["gene"])
    res = {}
    for scale, T in (("level3", R3), ("level5", Z5)):
        Y = T[idx]
        has = (~np.isnan(Y[:, :, 0])).sum(0) >= 4
        r = engine_E.run_panel(Y[:, has], lin, x, perms)
        df = pd.DataFrame(dict(cmap_name=np.array(comps_E)[has], T_c=r["T"], p_perm=r["p"], beta_gene=r["beta"][:, g]))
        df = df[np.isfinite(df.T_c)].copy()
        df["rank"] = df.T_c.rank(ascending=False, method="min").astype(int)
        res[scale] = (df, len(df))
    out = {}
    for scale, (df, n) in res.items():
        a = df[df.cmap_name == cfg.const("pc1")["compound"]].iloc[0]
        out[scale] = dict(rank=int(a["rank"]), n_ranked=n, p=float(a.p_perm), T_c=float(a.T_c), beta_gene=float(a.beta_gene))
    l3, l5 = out["level3"], out["level5"]
    pc1_pass = l3["p"] < cfg.const("pc1_p") and l3["rank"] / l3["n_ranked"] <= cfg.const("pc1_top_fraction") and l3["beta_gene"] < 0
    same_sign = np.sign(l3["beta_gene"]) == np.sign(l5["beta_gene"])
    return dict(gate="K4", scale="level3", tensor=tensor_name, panel="P1", compound=cfg.const("pc1")["compound"],
                n_lines=len(idx), n_g1=int(x.sum()), n_g0=int((1 - x).sum()), n_lineages=len(set(lin)),
                T_obs=l3["T_c"], p=l3["p"], rank=l3["rank"], n_ranked=l3["n_ranked"], beta_gene=l3["beta_gene"],
                beta_gene_level5=l5["beta_gene"], rank_level5=l5["rank"], p_level5=l5["p"], B=B,
                passes=bool(pc1_pass and same_sign))


def k4_reproduces(row, prior_l3, prior_l5):
    return bool(row["rank"] == prior_l3["rank"] and row["n_ranked"] == prior_l3["n_ranked"]
                and math.isclose(row["p"], prior_l3["p"], rel_tol=1e-12) and math.isclose(row["T_obs"], prior_l3["T"], rel_tol=1e-9)
                and math.isclose(row["beta_gene"], prior_l3["beta"], rel_tol=1e-9)
                and math.isclose(row["beta_gene_level5"], prior_l5["beta"], rel_tol=1e-9))


# ======================= calibration (pseudo-labels) =======================
_W = {}


def _init(panel_args):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    _W["P"] = C.Panel(*panel_args)


def xo_from_codes(P, tpX, codes_full, family, term):
    """Pseudo-label codes (NaN outside the file's lines) to the engine design (covariates first)."""
    fin = np.isfinite(codes_full)
    if family == "factorial":
        Y = np.where(fin, np.isin(codes_full, CODES_FACT["Y"]).astype(float), np.nan)
        X = np.where(fin, np.isin(codes_full, CODES_FACT["X"]).astype(float), np.nan)
        Xo = np.c_[tpX, Y, X]
        iy = tpX.shape[1]
        return Xo, (iy if term == "factY" else iy + 1), [iy, iy + 1]
    T = np.where(fin, (codes_full == 1).astype(float), np.nan)
    Xo = np.c_[tpX, T]
    return Xo, tpX.shape[1], [tpX.shape[1]]


def calib_rows(P, tpX, codes, elig, family, term, B, rule, s_idx, h1_subset=None):
    """D's calibration row for each pseudo-labelling s in s_idx (codes: (L, S) with NaN outside the file's lines)."""
    out = []
    for s in s_idx:
        Xo, k, tcols = xo_from_codes(P, tpX, codes[:, s], family, term)
        mask = C.analysis_set(P.lin, Xo, tcols, elig)
        df, permT, recs = C.run_term(P, Xo, k, mask, B, rule=rule)
        base = dict(term=term, pseudo=int(s), B=B, n_lines=int(mask.sum()), n_g1=int(np.nansum(Xo[mask, k] == 1)),
                    n_g0=int(np.nansum(Xo[mask, k] == 0)), n_lineages=len(set(P.lin[mask])))
        if "p" not in df.columns or not df.tested.any():
            base.update(n_tested=0, n_p_lt_005=np.nan, share_p_lt_005=np.nan, n_q_lt_01=np.nan, min_p=np.nan,
                        panel_stat=np.nan, panel_p=np.nan, n_untested_by_7iA_rule=0,
                        untested_reasons="; ".join(f"{a}: {b}" for a, b in df.reason.value_counts().items()))
            out.append(base); continue
        df, obs, pp, null = C.finish_table(df, permT)
        t = df[df.tested]
        bit = df[(~df.tested) & df.reason.str.contains("< 3 lines|relabellings", regex=True)]
        base.update(n_tested=len(t), n_p_lt_005=int((t.p < 0.05).sum()), share_p_lt_005=float((t.p < 0.05).mean()),
                    n_q_lt_01=int((t.q < 0.1).sum()), min_p=float(t.p.min()), panel_stat=obs, panel_p=pp,
                    n_untested_by_7iA_rule=len(bit),
                    untested_reasons="; ".join(f"{a}: {b}" for a, b in df[~df.tested].reason.value_counts().items()))
        if h1_subset is not None:
            sub = df.compound.isin(h1_subset(mask)).values
            o2, p2, _, n2 = C.panel_stat(df, permT, sub)
            base.update(panel_stat_h1set=o2, panel_p_h1set=p2, n_h1set=n2)
        out.append(base)
    return out


def _calib_task(args):
    tpX, codes, elig, family, term, B, rule, s_idx = args
    return calib_rows(_W["P"], tpX, codes, elig, family, term, B, rule, s_idx)


def run_calibration(panel_args, P, tpX, codes, elig, family, term, B, rule, s_idx, workers, chunk=10):
    """Runs calib_rows over s_idx, in `workers` processes (each single-threaded) when workers > 1."""
    s_idx = list(s_idx)
    if workers <= 1 or len(s_idx) <= chunk:
        return calib_rows(P, tpX, codes, elig, family, term, B, rule, s_idx)
    parts = [s_idx[i:i + chunk] for i in range(0, len(s_idx), chunk)]
    rows = []
    with ProcessPoolExecutor(max_workers=workers, initializer=_init, initargs=(panel_args,)) as ex:
        for r in ex.map(_calib_task, [(tpX, codes, elig, family, term, B, rule, p) for p in parts]):
            rows.extend(r)
    return sorted(rows, key=lambda r: r["pseudo"])


def summarize(rows, gate):
    """Task D's summary fields (src/D/08_summarize.py, as in task D2's 05_calibration_summary.py)."""
    g = pd.DataFrame(rows)
    t = g[g.n_tested > 0]
    r = dict(n_pseudo=len(g), n_pseudo_with_tests=len(t), n_lines_median=float(g.n_lines.median()),
             n_g1_median=float(g.n_g1.median()), n_g0_median=float(g.n_g0.median()),
             n_lineages_median=float(g.n_lineages.median()), n_tested_median=float(g.n_tested.median()))
    if len(t):
        s = t.share_p_lt_005
        r.update(mean_share_p_lt_005=float(s.mean()), median_share_p_lt_005=float(s.median()), p05_share=float(s.quantile(0.05)),
                 p95_share=float(s.quantile(0.95)), mean_n_p_lt_005=float(t.n_p_lt_005.mean()),
                 p95_n_p_lt_005=float(t.n_p_lt_005.quantile(0.95)), mean_n_q_lt_01=float(t.n_q_lt_01.mean()),
                 share_pseudo_any_q_lt_01=float((t.n_q_lt_01 > 0).mean()), p95_n_q_lt_01=float(t.n_q_lt_01.quantile(0.95)),
                 share_pseudo_panel_p_lt_005=float((t.panel_p < 0.05).mean()),
                 panel_stat_p05=float(t.panel_stat.quantile(0.05)), panel_stat_p95=float(t.panel_stat.quantile(0.95)),
                 gate_pass=bool(gate[0] <= s.mean() <= gate[1]))
    else:
        r.update(gate_pass=False)
    return r


COMPARE = ["n_lines", "n_g1", "n_g0", "n_lineages", "n_tested", "n_p_lt_005", "n_q_lt_01"]


def compare_rows(new, prior):
    """Exact comparison of per-pseudo-labelling rows (counts equal; p values equal; panel statistic to rel 1e-9)."""
    a = pd.DataFrame(new).set_index("pseudo")
    b = prior.set_index("pseudo").loc[a.index]
    same = np.ones(len(a), bool)
    for c in COMPARE:
        same &= (a[c].fillna(-1).values == b[c].fillna(-1).values)
    for c in ("min_p", "panel_p"):
        same &= np.isclose(a[c].fillna(-1).values, b[c].fillna(-1).values, rtol=1e-12, atol=0)
    same &= np.isclose(a.panel_stat.fillna(-1).values, b.panel_stat.fillna(-1).values, rtol=1e-9, atol=0)
    return pd.Series(same, index=a.index)


def fresh_codes(lines_sorted, counts, S, seed):
    """Task D's pseudo-label scheme: lines sorted by lineage and name; within each lineage a uniformly random
    arrangement of the lineage's group codes, exactly the given counts. Returns (lines, lineages, codes (n, S))."""
    rng = np.random.default_rng(seed)
    e = pd.DataFrame(lines_sorted, columns=["cell_iname", "lineage"]).sort_values(["lineage", "cell_iname"]).reset_index(drop=True)
    c = counts.set_index("lineage")
    lab = np.zeros((len(e), S), dtype=np.int8)
    for lin, idx in e.groupby("lineage").groups.items():
        r = c.loc[lin]
        cd = np.concatenate([np.full(int(r[f"n_g{g}"]), g) for g in (1, 2, 3, 4)])
        assert len(cd) == len(idx)
        for s in range(S):
            lab[np.asarray(idx), s] = rng.permutation(cd)
    return e.cell_iname.tolist(), e.lineage.tolist(), lab


def codes_full(P, file_lines, codes):
    pos = {l: i for i, l in enumerate(P.lines)}
    full = np.full((len(P.lines), codes.shape[1]), np.nan)
    ix = np.array([pos[l] for l in file_lines])
    full[ix] = codes
    elig = np.zeros(len(P.lines), bool); elig[ix] = True
    return full, elig


def read_pseudo(path):
    ps = pd.read_csv(path, sep="\t")
    C.guard_columns([c for c in ps.columns if c not in ("cell_iname", "lineage")])
    cols = [c for c in ps.columns if c.startswith("s")]
    return ps.cell_iname.tolist(), ps.lineage.tolist(), ps[cols].values.astype(float)


def implied_set(P, tpX, file_lines, codes, family):
    """The analysis set and per-lineage counts implied by a pseudo-label file (from its first pseudo-labelling)."""
    full, elig = codes_full(P, file_lines, codes[:, :1])
    Xo, k, tcols = xo_from_codes(P, tpX, full[:, 0], family, "factY" if family == "factorial" else "annotated")
    m = C.analysis_set(P.lin, Xo, tcols, elig)
    cnt = []
    for l in sorted(set(P.lin[m])):
        v = full[m & (P.lin == l), 0]
        cnt.append(dict(lineage=l, n_lines=len(v), **{f"n_g{g}": int(np.sum(v == g)) for g in (1, 2, 3, 4)}))
    return [P.lines[i] for i in np.flatnonzero(m)], pd.DataFrame(cnt)


def same_counts(a, b):
    """Equal per-lineage counts (missing counts read as 0, numbers compared as integers)."""
    cols = ["n_lines", "n_g1", "n_g2", "n_g3", "n_g4"]
    f = lambda t: t.assign(**{c: pd.to_numeric(t[c], errors="coerce").fillna(0).astype(int) for c in cols}).sort_values("lineage")[
        ["lineage"] + cols].astype(str).reset_index(drop=True)
    return f(a).equals(f(b))


# ======================= label-free E3 precision (task D2's scheme) =======================
def _fit_d2(Sv, lin, tp, G):
    """Task D2's fit (from task E's 07_response_size.py): OLS of S on lineage + TP53 (mutant, unknown) + G, classical
    and HC3 SE of the G columns."""
    levels = sorted(set(lin))
    cols = [(lin == l).astype(float) for l in levels]
    for t in ["mutant", "unknown"]:
        cols.append((tp == t).astype(float))
    X = np.column_stack(cols)
    keep = []
    for j in range(X.shape[1]):
        if np.linalg.matrix_rank(X[:, keep + [j]]) == len(keep) + 1:
            keep.append(j)
    q = G.shape[1]
    X = np.column_stack([X[:, keep], G.astype(float)])
    if np.linalg.matrix_rank(X) < X.shape[1]:
        return np.full(q, np.nan), np.full(q, np.nan), np.full(q, np.nan), np.nan
    XtXi = np.linalg.inv(X.T @ X)
    A = XtXi @ X.T
    beta = A @ Sv
    e = Sv - X @ beta
    n, kk = X.shape
    df = n - kk
    se = np.sqrt(e @ e / df * np.diag(XtXi)[-q:])
    h = np.einsum("ij,ji->i", X, A)
    ok = h < 1 - 1e-8
    se3 = np.array([np.sqrt((np.where(ok, (A[-q + j] / np.where(ok, 1 - h, 1)) ** 2, 0.0)) @ (e ** 2)) for j in range(q)])
    return beta[-q:], se, se3, df


def precision_e3(S_all, panel_lines, set_lines, set_lineage, set_tp, counts, factorial, draws, seed, margin, levels=(15, 20, 25)):
    """Task D2's label-free precision (src/D2/E/07_precision.py): within-lineage random arrangements of the set's
    counts; S ~ lineage + TP53 + term(s); estimates in percent of the mean S of the model lines."""
    pos = {l: i for i, l in enumerate(panel_lines)}
    ix = np.sort(np.array([pos[l] for l in set_lines]))
    order = [panel_lines[i] for i in ix]
    meta = pd.DataFrame(dict(cell_iname=set_lines, lineage=set_lineage, tp=set_tp)).set_index("cell_iname").loc[order]
    lin, tp = meta.lineage.values, meta.tp.values
    Sv = S_all[ix]
    cnt = counts.set_index("lineage")
    assert set(cnt.index) == set(lin) and int(cnt.n_lines.sum()) == len(ix)
    rng = np.random.default_rng(seed)
    est, hw, hw3 = [], [], []
    for _ in range(draws):
        grp = np.zeros(len(ix), int)
        for l, r in cnt.iterrows():
            m = np.nonzero(lin == l)[0]
            ns = [int(r[f"n_g{k}"]) if pd.notna(r[f"n_g{k}"]) else 0 for k in (1, 2, 3, 4)]
            pick = rng.permutation(m)[:sum(ns)]
            assert len(pick) == len(m)
            grp[pick] = np.repeat([1, 2, 3, 4], ns)
        G = np.column_stack([np.isin(grp, [1, 4]), np.isin(grp, [2, 4])]) if factorial else (grp == 1)[:, None]
        b, se, se3, df = _fit_d2(Sv, lin, tp, G)
        f100 = 100.0 / Sv.mean()
        qt = stats.t.ppf(0.95, df)
        est.append(b * f100); hw.append(qt * se * f100); hw3.append(qt * se3 * f100)
    est, hw, hw3 = map(np.array, (est, hw, hw3))
    terms = ["Y_plus", "X_plus"] if factorial else ["g1_vs_g2"]
    rows = []
    for j, term in enumerate(terms):
        e_, h_, h3_ = est[:, j], hw[:, j], hw3[:, j]
        ok = np.isfinite(e_) & np.isfinite(h_) & np.isfinite(h3_)
        e_, h_, h3_ = e_[ok], h_[ok], h3_[ok]
        r = dict(term=term, n_draws=draws, n_draws_estimable=int(ok.sum()), n_lines_in_model=len(ix), n_lineages=len(set(lin)),
                 mean_S_lines_in_model=float(Sv.mean()), median_halfwidth_pct=float(np.median(h_)),
                 median_halfwidth_hc3_pct=float(np.median(h3_)), sd_estimate_pct=float(e_.std(ddof=1)))
        for X in levels:
            r[f"P_ci_inside_pm{X}pct"] = float((np.abs(e_) + h_ < X).mean())
            r[f"P_ci_inside_pm{X}pct_hc3"] = float((np.abs(e_) + h3_ < X).mean())
        rows.append(r)
    return rows
