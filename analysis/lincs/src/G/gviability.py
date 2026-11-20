"""Viability step of task G (task draft step 10; REC 4d H2 rule 4; PR H8). Panel P1, factorial Y and X terms.

1. E2 with and without the compound's PRISM 19Q4 AUC as a covariate, on the lines that have an AUC.
2. Viability state effect per compound in the same lines with task C's model (OLS of AUC on lineage fixed
   effects and the state terms, no other covariate; HC3).
3. Across compounds: Spearman correlation of E2 (same lines, unadjusted) with the viability effect, with a
   relabelling null: the (Y, X) labels are permuted jointly within lineage over the analysis set (first N rows of
   the panel's permutation keys), and both effects are recomputed. Descriptive.
"""
import numpy as np, pandas as pd
from scipy import stats
import gcore as C
import engine as E
from ganalyses import express


def _coef(lin, X, y, k):
    """OLS coefficient k of y on lineage fixed effects + X (within-lineage centring)."""
    Xd = E.centre(X, lin)
    yd = E.centre(y[:, None], lin)[:, 0]
    b, *_ = np.linalg.lstsq(Xd, yd, rcond=None)
    return b[k]


def run(ctx, res, B_relabel):
    cfg, P, d = ctx.cfg, ctx.P["P1"], res["design"]
    auc = ctx.auc.get("P1")
    if auc is None:
        return pd.DataFrame(), dict(note="no PRISM data in this run")
    s = ctx.s_level5["P1"]
    tcols = d["tcols"]
    k_in_t = tcols.index(d["k"])
    rows_out, keep = [], []
    for r in res["recs"]:
        if not r["tested"] or r.get("_j", -1) < 0:
            continue
        j = r["_j"]
        a = auc[:, j]
        rows = r["_rows"][np.isfinite(a[r["_rows"]])]
        row = dict(term=d["term"], compound=r["compound"], n_lines_with_auc=len(rows))
        f0 = C.scalar_fit2(d["lin"], d["Xo"], d["k"], r["_rows"], s[:, j], rule=d["rule"])
        e0 = express(f0, cfg)
        if len(rows) < 4:
            row.update(E2_primary=e0[0], note="fewer than 4 lines with an AUC"); rows_out.append(row); continue
        f1 = C.scalar_fit2(d["lin"], d["Xo"], d["k"], rows, s[:, j], rule=d["rule"])
        Xa = np.c_[a, d["Xo"]]
        f2 = C.scalar_fit2(d["lin"], Xa, d["k"] + 1, rows, s[:, j], rule=d["rule"])
        Xv = d["Xo"][:, tcols]
        fv = C.scalar_fit2(d["lin"], Xv, k_in_t, rows, a, rule=d["rule"])
        e1, e2 = express(f1, cfg), express(f2, cfg)
        row.update(E2_primary=e0[0], E2_same_lines=e1[0], E2_same_lines_lo90=e1[1], E2_same_lines_hi90=e1[2],
                   E2_auc_adjusted=e2[0], E2_auc_adjusted_lo90=e2[1], E2_auc_adjusted_hi90=e2[2],
                   E2_change_more_than_half=bool(np.isfinite(e1[0]) and np.isfinite(e2[0]) and abs(e2[0] - e1[0]) > 0.5 * abs(e1[0])),
                   viability_beta=fv["beta"] if fv else np.nan, viability_lo90=fv["beta"] - fv["hw90_hc3"] if fv else np.nan,
                   viability_hi90=fv["beta"] + fv["hw90_hc3"] if fv else np.nan, note="")
        rows_out.append(row)
        if f1 is not None and fv is not None:
            keep.append((r["compound"], f1["rows"], j, e1[0], fv["beta"], f1["mean_s"]))
    tab = pd.DataFrame(rows_out)
    summ = dict(term=d["term"], n_compounds=len(keep))
    if len(keep) >= 5:
        e2v = np.array([k[3] for k in keep]); vb = np.array([k[4] for k in keep])
        rho = float(stats.spearmanr(e2v, vb).statistic)
        m = np.flatnonzero(d["mask"])
        perm = E.induced_perm(P.keys[:B_relabel][:, m], d["lin"][m])
        lab = d["Xo"][:, tcols].copy()
        null = np.empty(B_relabel)
        unit = cfg.decide("7j")["unit"] == "A"
        for b in range(B_relabel):
            labb = lab.copy()
            labb[m] = lab[m][perm[b]]
            Xb = d["Xo"].copy(); Xb[:, tcols] = labb
            ev, vv = [], []
            for c, rows, j, _, _, ms in keep:
                lr = d["lin"][rows]
                e = _coef(lr, Xb[rows], s[rows, j], d["k"]) * (100.0 / ms if unit else 100.0)
                v = _coef(lr, labb[rows], auc[rows, j], k_in_t)
                ev.append(e); vv.append(v)
            null[b] = stats.spearmanr(ev, vv).statistic
        summ.update(spearman_E2_vs_viability=rho, relabel_n=B_relabel, null_p025=float(np.quantile(null, 0.025)),
                    null_p975=float(np.quantile(null, 0.975)), null_mean=float(np.mean(null)),
                    p_two_sided=float((1 + np.sum(np.abs(null) >= abs(rho))) / (B_relabel + 1)))
    return tab, summ
