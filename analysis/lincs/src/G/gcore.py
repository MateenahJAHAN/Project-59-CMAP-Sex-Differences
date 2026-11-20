"""Core helpers for task G.

Adapted from task D's src/D/common.py (origin/task-D a84f566, sha256 1fcea3e2...). The engine itself
(src/G/engine.py) is task D's final engine.py, byte-identical (sha256 2925ff6d...).

Changes from task D's common.py, each marked "G:" below:
  1. Paths and permutation keys are arguments (no fixed task-D layout).
  2. compound_design() takes the minimum group size and the minimum number of distinct within-lineage
     relabellings (DECIDE 7i). With min_group=2 and min_relab=0 it is task D's function.
  3. run_term() takes an optional block vector `lin` (used by S8 and S9) and a `rule` dict; with neither it is
     task D's function.
  4. New: run_term_signed() gives per-permutation reductions of the signed t vector and of the coefficient
     vector of the tested term (X-linked score, one-sided scalar tests, P2 shape agreement).
  5. New: panel_stat() for the panel statistic over any compound subset, and scalar_fit2() with HC3 and
     classical intervals at any level.
  6. log() never prints a statistic: callers pass progress messages only.
"""
import os, time, math
import numpy as np, pandas as pd
from scipy import stats
import engine as E

FORBIDDEN = ["state", "sex", "gender", "yhi", "xhi", "xist", "y_score", "conflict", "confidence", "fact_", "strict_xy", "annot_m"]
TIE_TOL = 1e-5
_LOG = {"path": None}


def set_log(path):
    _LOG["path"] = path


def log(tag, msg, t0=None):
    el = f" +{time.time()-t0:.0f}s" if t0 is not None else ""
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} [{tag}{el}] {msg}"
    print(line, flush=True)
    if _LOG["path"]:
        with open(_LOG["path"], "a") as f:
            f.write(line + "\n")


def guard_columns(cols):
    bad = [c for c in cols if any(k in c.lower() for k in FORBIDDEN)]
    if bad:
        raise RuntimeError(f"BLINDING GUARD: forbidden columns {bad}")


def read_covariates(path):
    hdr = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    guard_columns(hdr)
    df = pd.read_csv(path, sep="\t")
    assert df.cell_iname.is_unique
    return df


class Panel:
    """G: task D's Panel with explicit paths. R: (lines, compounds, 978) float64; keys: (B_full, lines)."""

    def __init__(self, pn, tensor_path, lines_path, comps_path, genes, cov, key_seed, B_full, dmso_path=None):
        self.pn = pn
        self.R = np.load(tensor_path).astype(np.float64)
        self.lines = pd.read_csv(lines_path, sep="\t").cell_iname.tolist()
        self.comps = pd.read_csv(comps_path, sep="\t").cmap_name.tolist()
        assert self.R.shape == (len(self.lines), len(self.comps), 978)
        self.genes = list(genes)
        c = pd.DataFrame({"cell_iname": self.lines}).merge(cov, on="cell_iname", how="left", validate="1:1")
        assert len(c) == len(self.lines) and c.lineage.notna().all()
        self.cov = c
        self.lin = np.array(c.lineage.astype(str).tolist(), dtype="U64")
        self.avail = np.isfinite(self.R[:, :, 0])
        self.keys = E.make_keys(len(self.lines), B_full, key_seed)
        self.dmso = np.load(dmso_path).astype(np.float64) if dmso_path and os.path.exists(dmso_path) else None


def analysis_set(lin, Xo, tested_cols, base_mask):
    """Task D's rule: lines with all model columns known (base_mask & finite Xo), in lineages with >= 2 such lines
    in which at least one tested column varies."""
    ok = base_mask & np.all(np.isfinite(Xo), axis=1)
    keep = np.zeros_like(ok)
    for l in np.unique(lin[ok]):
        idx = ok & (lin == l)
        if idx.sum() >= 2 and any(np.unique(Xo[idx, k]).size > 1 for k in tested_cols):
            keep |= idx
    return keep


def n_relabellings(lin, t):
    """G: number of distinct within-lineage arrangements of the tested values t (product of multinomials)."""
    n = 1
    for l in np.unique(lin):
        v = t[lin == l]
        _, cnt = np.unique(v, return_counts=True)
        m = math.factorial(len(v))
        for c in cnt:
            m //= math.factorial(int(c))
        n *= m
    return n


def compound_design(lin, Xo, k, rows, min_group=2, min_relab=0):
    """Task D's per-compound design (drop lineages with < 2 available lines; untested if a group has fewer than
    min_group lines, the design is singular or a line has leverage 1). G: min_group and min_relab (7i)."""
    l = lin[rows]
    u, cnt = np.unique(l, return_counts=True)
    rows = rows[np.isin(l, u[cnt >= 2])]
    if len(rows) < 4:
        return rows, None, "too few lines"
    t = Xo[rows, k]
    vals = np.unique(t)
    if vals.size == 2 and min((t == vals[0]).sum(), (t == vals[1]).sum()) < min_group:
        return rows, None, f"group with < {min_group} lines"
    if vals.size < 2:
        return rows, None, "tested term constant"
    if min_relab > 0 and n_relabellings(lin[rows], t) < min_relab:
        return rows, None, f"fewer than {min_relab} distinct within-lineage relabellings"
    try:
        d = E.Design(lin[rows], Xo[rows], k)
    except ValueError as ex:
        return rows, None, str(ex)
    if d.df < 1:
        return rows, None, "no residual df"
    return rows, d, ""


def _rule(rule):
    rule = rule or {}
    return int(rule.get("min_group", 2)), int(rule.get("min_relab", 0))


def run_term(P, Xo, k, mask, B, comps=None, Y_override=None, want_perm=True, dtype=np.float32, extra=None,
             lin=None, rule=None):
    """Task D's run_term. G: `lin` replaces the lineage blocks (S8, S9); `rule` sets 7i (see compound_design).
    Returns (table, permT (B x targets) or None, recs)."""
    lin = P.lin if lin is None else lin
    mg, mr = _rule(rule)
    keys = P.keys[:B]
    rows_all = np.flatnonzero(mask)
    targets = []
    cidx = range(len(P.comps)) if comps is None else comps
    for j in cidx:
        Y = P.R[:, j, :] if (Y_override is None or j not in Y_override) else Y_override[j]
        targets.append((P.comps[j], j, Y))
    if extra:
        for nm, Y in extra.items():
            targets.append((nm, -1, Y))
    recs, groups = [], {}
    for name, j, Y in targets:
        av = np.isfinite(Y[:, 0])
        rows = rows_all[av[rows_all]]
        rows, d, why = compound_design(lin, Xo, k, rows, mg, mr)
        if d is None:
            recs.append(dict(compound=name, tested=False, reason=why, n_lines=len(rows)))
            continue
        b, se, t, sig = d.observed(Y[rows])
        tv = Xo[rows, k]
        g1 = int((tv == tv.max()).sum())
        rec = dict(compound=name, tested=True, reason="", n_lines=len(rows), n_g1=g1, n_g0=len(rows) - g1,
                   n_lineages=d.n_fe, df=d.df, T_obs=float(np.mean(t ** 2)))
        rec["_b"], rec["_se"], rec["_t"], rec["_sig"], rec["_d"], rec["_rows"], rec["_j"] = b, se, t, sig, d, rows, j
        recs.append(rec)
        if want_perm:
            key = tuple(rows)
            groups.setdefault(key, []).append((len(recs) - 1, d, d.reduced_resid(Y[rows])))
    permT = None
    tested = [i for i, r in enumerate(recs) if r["tested"]]
    if want_perm and tested:
        permT = np.full((B, len(recs)), np.nan)
        for key, items in groups.items():
            rows = np.array(key)
            d = items[0][1]
            perm = E.induced_perm(keys[:, rows], lin[rows])
            Ecat = np.concatenate([it[2] for it in items], axis=1)
            G = items[0][2].shape[1]
            chunk = int(max(8, min(1000, 4e7 / (2 * d.p * Ecat.shape[1]))))
            out = d.perm_T(Ecat, perm, dtype=dtype, chunk=chunk, n_genes=G)
            for m, it in enumerate(items):
                permT[:, it[0]] = out[:, m]
        for i in tested:
            recs[i]["p"] = (1 + np.sum(permT[:, i] >= recs[i]["T_obs"] * (1 - TIE_TOL))) / (B + 1)
    df = pd.DataFrame([{k2: v for k2, v in r.items() if not k2.startswith("_")} for r in recs])
    return df, permT, recs


def perm_signed(d, Eres, perm, reducer, chunk=200, dtype=np.float64):
    """G: for the tested coefficient k of Design d, the signed HC3 t vector and the coefficient vector under each
    permutation (same algebra as engine.Design.perm_T). reducer(t, beta) -> (b, ...) array per chunk.
    Eres: (n, G) reduced residuals; perm: (B, n)."""
    B = perm.shape[0]
    n, G = Eres.shape
    inv = np.argsort(perm, axis=1)
    E_ = np.ascontiguousarray(Eres, dtype=dtype)
    E2 = np.ascontiguousarray(Eres ** 2, dtype=dtype)
    WX = (d.w[:, None] * d.Xd).T
    Rm = 2.0 * WX - d.K @ d.A
    p, k = d.p, d.k
    M = np.vstack([d.A, Rm, d.w[None, :]]).astype(dtype)
    outs = []
    for s in range(0, B, chunk):
        iv = inv[s:s + chunk]
        b = iv.shape[0]
        Mb = M[:, iv].transpose(1, 0, 2)
        UR = (Mb[:, :2 * p].reshape(b * 2 * p, n) @ E_).reshape(b, 2 * p, G)
        S1 = Mb[:, 2 * p] @ E2
        V = S1 - np.einsum("bjg,bjg->bg", UR[:, :p], UR[:, p:])
        beta = UR[:, k]
        t = np.divide(beta, np.sqrt(np.maximum(V, 0.0)), out=np.zeros_like(V), where=V > 0)
        outs.append(np.asarray(reducer(t, beta)))
    return np.concatenate(outs, axis=0)


def run_term_signed(P, Xo, k, mask, B, reducer_for, comps=None, Y_override=None, lin=None, rule=None):
    """G: like run_term, but returns, per tested target, reducer(t, beta) over B permutations (in recs[i]['perm'])
    and reducer applied to the observed t and beta (recs[i]['obs']). reducer_for(name, rec) returns the reducer
    for one target (or None to skip its permutations)."""
    lin = P.lin if lin is None else lin
    mg, mr = _rule(rule)
    keys = P.keys[:B]
    rows_all = np.flatnonzero(mask)
    recs = []
    cidx = range(len(P.comps)) if comps is None else comps
    for j in cidx:
        Y = P.R[:, j, :] if (Y_override is None or j not in Y_override) else Y_override[j]
        av = np.isfinite(Y[:, 0])
        rows = rows_all[av[rows_all]]
        rows, d, why = compound_design(lin, Xo, k, rows, mg, mr)
        if d is None:
            recs.append(dict(compound=P.comps[j], tested=False, reason=why, n_lines=len(rows)))
            continue
        b, se, t, sig = d.observed(Y[rows])
        rec = dict(compound=P.comps[j], tested=True, reason="", n_lines=len(rows), n_lineages=d.n_fe, df=d.df,
                   _b=b, _t=t, _se=se, _sig=sig, _rows=rows, _j=j, _d=d)
        red = reducer_for(P.comps[j], rec)
        if red is not None:
            rec["obs"] = np.asarray(red(t[None, :], b[None, :]))[0]
            perm = E.induced_perm(keys[:, rows], lin[rows])
            rec["perm"] = perm_signed(d, d.reduced_resid(Y[rows]), perm, red)
        recs.append(rec)
    return recs


def bh(p):
    p = np.asarray(p, dtype=float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    o = np.argsort(pv)
    r = pv[o] * n / np.arange(1, n + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    qq = np.empty(n); qq[o] = np.minimum(r, 1)
    q[ok] = qq
    return q


def finish_table(df, permT, tested_only_mean=True, compounds_only=None):
    """Task D's: adds q and rank among tested panel compounds; returns (df, panel stat, panel p, null)."""
    df = df.copy()
    t = df.tested.values.astype(bool)
    if compounds_only is not None:
        t = t & df.compound.isin(compounds_only).values
    df["in_ranking"] = t
    df["q"] = np.nan
    df.loc[t, "q"] = bh(df.loc[t, "p"].values)
    df["rank"] = np.nan
    df.loc[t, "rank"] = df.loc[t, "T_obs"].rank(ascending=False, method="min").values
    df["n_ranked"] = int(t.sum())
    obs = float(df.loc[t, "T_obs"].mean())
    null = np.nanmean(permT[:, t], axis=1) if permT is not None else None
    pp = float((1 + np.sum(null >= obs * (1 - TIE_TOL))) / (len(null) + 1)) if null is not None else np.nan
    return df, obs, pp, null


def panel_stat(df, permT, subset_mask):
    """G: panel statistic (mean T_c) over a boolean subset of the table rows, with its permutation p."""
    sel = subset_mask & df.tested.values.astype(bool)
    if sel.sum() == 0:
        return np.nan, np.nan, None, 0
    obs = float(df.loc[sel, "T_obs"].mean())
    null = np.nanmean(permT[:, sel], axis=1)
    pp = float((1 + np.sum(null >= obs * (1 - TIE_TOL))) / (len(null) + 1))
    return obs, pp, null, int(sel.sum())


def response_size(P_R, avail):
    """Task D's s_ic = r_ic . m_c,-i / |m_c,-i|^2, m_c,-i the mean of the other available panel lines (label-free)."""
    L, C, G = P_R.shape
    S = np.full((L, C), np.nan)
    for j in range(C):
        av = np.flatnonzero(avail[:, j])
        Y = P_R[av, j, :]
        tot = Y.sum(axis=0)
        n = len(av)
        for a, i in enumerate(av):
            m = (tot - Y[a]) / (n - 1)
            S[i, j] = float(Y[a] @ m / (m @ m))
    return S


def scalar_fit2(lin, Xo, k, rows, s, rule=None, levels=(0.90, 0.95)):
    """G: OLS of a scalar per-line outcome s on the design (task D's compound_design rows), with HC3 and classical
    SE and CIs at the given levels. Returns dict or None."""
    mg, mr = _rule(rule)
    rows = rows[np.isfinite(s[rows])]
    rows, d, why = compound_design(lin, Xo, k, rows, mg, mr)
    if d is None:
        return None
    b, se, t, sig = d.observed(s[rows][:, None])
    se_cl = float(sig[0] * np.sqrt(d.xtxi_kk))
    out = dict(beta=float(b[0]), se_hc3=float(se[0]), se_classical=se_cl, df=d.df, n=len(rows), t_hc3=float(t[0]),
               p_hc3=float(2 * stats.t.sf(abs(t[0]), d.df)), mean_s=float(np.mean(s[rows])), rows=rows)
    for lv in levels:
        tq = stats.t.ppf(0.5 + lv / 2, d.df)
        out[f"hw{int(round(lv*100))}_hc3"] = float(tq * se[0])
        out[f"hw{int(round(lv*100))}_classical"] = float(tq * se_cl)
    return out
