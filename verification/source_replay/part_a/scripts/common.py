"""Shared constants and the frozen rules (plan sections 3-7)."""
import numpy as np, pandas as pd
from scipy.stats import gaussian_kde

SEED = 20261115
DEEP_MIN = 2500        # plan 7d
NBINS = 20             # random bins of marker-positive deep cells (plan 7d)
NGROUP = 2 + 2 * NBINS # groups per line in the identity accumulator
MIN_UNIT_CELLS = 10    # plan 3
Y_GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y", "NLGN4Y"]
YP_GENES = ["RPS4Y1", "ZFY"]
YQ_GENES = ["DDX3Y", "USP9Y", "UTY", "NLGN4Y", "KDM5D", "EIF1AY"]
X_GAMETOLOGS = ["RPS4X", "DDX3X", "KDM5C", "KDM6A", "EIF1AX", "ZFX", "USP9X", "NLGN4X"]
MARKERS = ["XIST", "RPS4Y1", "ZFY", "DDX3Y", "USP9Y", "UTY", "NLGN4Y", "KDM5D", "EIF1AY"] + X_GAMETOLOGS
assert sorted(MARKERS[1:9]) == sorted(Y_GENES) and sorted(YP_GENES + YQ_GENES) == sorted(Y_GENES)
CCLE_CUT = {"XIST": 1.983388, "y_score": 1.306546, "yp_score": 1.306546, "yq_score": 1.306546}
STATES = {"XY", "XX-like", "Y-lost", "XIST-lost", "discordant", "ambiguous"}
ARM_STATES = {"Y intact", "Yq lost", "Yp lost", "Y lost"}


def l2cpm(counts, depth):
    return np.log2(1e6 * np.asarray(counts, float) / np.asarray(depth, float) + 1)


def pseudobulk_scores(df, by):
    """df: per-cell table with marker counts and depth. Returns XIST, y/yp/yq scores per group."""
    g = df.groupby(by, observed=True)[["depth"] + MARKERS].sum()
    n = df.groupby(by, observed=True).size().rename("n_cells")
    out = pd.DataFrame(index=g.index)
    out["n_cells"] = n
    out["depth"] = g["depth"]
    for s in MARKERS:
        out[s + "_l2cpm"] = l2cpm(g[s], g["depth"])
    out["XIST"] = out["XIST_l2cpm"]
    out["y_score"] = out[[s + "_l2cpm" for s in Y_GENES]].mean(axis=1)
    out["yp_score"] = out[[s + "_l2cpm" for s in YP_GENES]].mean(axis=1)
    out["yq_score"] = out[[s + "_l2cpm" for s in YQ_GENES]].mean(axis=1)
    return out.reset_index()


def kde_cutoff(values):
    """Plan 4: Scott-bandwidth Gaussian KDE on a 1,024-point grid spanning the observed range;
    cut-off at the density minimum between the two tallest modes. Returns (cut, info dict)."""
    v = np.asarray(values, float); v = v[np.isfinite(v)]
    grid = np.linspace(v.min(), v.max(), 1024)
    kde = gaussian_kde(v, bw_method="scott")
    d = kde(grid)
    modes = []
    for i in range(len(d)):
        left = d[i] > d[i - 1] if i > 0 else d[i] > d[i + 1]
        right = d[i] >= d[i + 1] if i < len(d) - 1 else d[i] > d[i - 1]
        if (i == 0 and left) or (i == len(d) - 1 and right) or (0 < i < len(d) - 1 and left and right):
            modes.append(i)
    info = dict(n_units=len(v), bandwidth=float(kde.factor * np.std(v, ddof=1)), n_modes=len(modes),
                grid_min=float(grid[0]), grid_max=float(grid[-1]))
    if len(modes) < 2:
        info["method"] = "fallback_ccle"
        return np.nan, info
    top = sorted(sorted(modes, key=lambda i: -d[i])[:2])
    j = top[0] + int(np.argmin(d[top[0]:top[1] + 1]))
    info.update(method="kde", mode1=float(grid[top[0]]), mode2=float(grid[top[1]]),
                mode1_density=float(d[top[0]]), mode2_density=float(d[top[1]]), cut_density=float(d[j]))
    return float(grid[j]), info


def call_state(ann, xist, y, cx, cy):
    """Plan 5 rule table (from Task A output)."""
    if not (np.isfinite(xist) and np.isfinite(y)):
        return "ambiguous"
    yhi, xhi = y >= cy, xist >= cx
    female, male = ann == "Female", ann == "Male"
    if yhi:
        return "discordant" if female else "XY"
    if xhi:
        return "discordant" if male else "XX-like"
    if female:
        return "XIST-lost"
    if male:
        return "Y-lost"
    return "ambiguous"


def call_arm(yp, yq, cp, cq):
    ph, qh = yp >= cp, yq >= cq
    return {(True, True): "Y intact", (True, False): "Yq lost", (False, True): "Yp lost", (False, False): "Y lost"}[(ph, qh)]


def em_mixture(y, d, b, tol=1e-9, maxit=5000):
    """Plan 7a: y ~ pi Pois(b d) + (1-pi) Pois(lam d), b fixed per cell (its batch background).
    Returns dict(pi, lam, loglik, loglik1, iters)."""
    from scipy.special import gammaln, xlogy
    y = np.asarray(y, float); d = np.asarray(d, float); b = np.asarray(b, float)
    mu0 = b * d
    lg = gammaln(y + 1)
    l0 = xlogy(y, mu0) - mu0 - lg
    lam = 2 * y.sum() / d.sum(); pi = 0.5; ll_old = -np.inf
    for it in range(maxit):
        mu1 = lam * d
        l1 = xlogy(y, mu1) - mu1 - lg
        a0 = np.log(pi) + l0; a1 = np.log1p(-pi) + l1
        mx = np.maximum(a0, a1)
        lse = mx + np.log(np.exp(a0 - mx) + np.exp(a1 - mx))
        ll = lse.sum()
        r0 = np.exp(a0 - lse)
        pi = min(max(r0.mean(), 1e-12), 1 - 1e-12)
        w = 1 - r0
        lam = max((w * y).sum() / max((w * d).sum(), 1e-300), 1e-300)
        if abs(ll - ll_old) <= tol * abs(ll):
            break
        ll_old = ll
    lam1 = y.sum() / d.sum()
    ll1 = (xlogy(y, lam1 * d) - lam1 * d - lg).sum()
    return dict(pi=float(pi), lam=float(lam), loglik=float(ll), loglik_one=float(ll1), iters=it + 1)


def mh_or(e, f, strata):
    """Mantel-Haenszel OR of binary e with binary f over strata."""
    num = den = 0.0
    for s in np.unique(strata):
        m = strata == s
        a = np.sum(e[m] & f[m]); b = np.sum(e[m] & ~f[m]); c = np.sum(~e[m] & f[m]); dd = np.sum(~e[m] & ~f[m])
        n = a + b + c + dd
        if n == 0:
            continue
        num += a * dd / n; den += b * c / n
    return float(num / den) if den > 0 else (np.inf if num > 0 else np.nan)
