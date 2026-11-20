"""Per-compound OLS with lineage fixed effects + one tested term, HC3 t per gene,
T_c = mean over genes of t^2, Freedman-Lane permutation within lineage shared across compounds.
"""
import numpy as np


def design(lineage, x):
    """lineage: array of labels (n,), x: tested term (n,). Returns X (n, k+1), Z (n, k)."""
    levels = sorted(set(lineage))
    Z = np.stack([(np.asarray(lineage) == l).astype(float) for l in levels], 1)
    X = np.column_stack([Z, np.asarray(x, float)])
    return X, Z


def fixed_parts(X, Z):
    A = np.linalg.solve(X.T @ X, X.T)          # (p, n)
    a = A[-1]                                   # tested-term row
    H = X @ A
    h = np.diag(H)
    M = np.eye(len(X)) - H
    H0 = Z @ np.linalg.solve(Z.T @ Z, Z.T)
    return a, h, M, H0


def tstats(Yperm_resid, a, h, M):
    """Yperm_resid: (n, m) reduced-model residuals, already permuted. Full-model coefficient and
    HC3 t of the tested term for each column. Valid because the reduced fit lies in col(X)."""
    beta = a @ Yperm_resid
    e = M @ Yperm_resid
    w = (a / (1.0 - h)) ** 2
    var = w @ (e ** 2)
    return beta, beta / np.sqrt(var)


def within_lineage_perms(lineage, B, rng):
    """B permutations of positions 0..n-1 within lineage. perms[b, j] = source index for position j."""
    lineage = np.asarray(lineage)
    n = len(lineage)
    perms = np.empty((B, n), dtype=np.int64)
    blocks = [np.nonzero(lineage == l)[0] for l in sorted(set(lineage))]
    for b in range(B):
        p = np.arange(n)
        for idx in blocks:
            p[idx] = rng.permutation(idx)
        perms[b] = p
    return perms


def induced(perms, lineage, obs):
    """Restrict full permutations to observed positions `obs` (sorted index array into the full set).
    Returns (B, len(obs)) permutations expressed in local indices 0..len(obs)-1."""
    lineage = np.asarray(lineage)
    obs = np.asarray(obs)
    loc = -np.ones(len(lineage), dtype=np.int64); loc[obs] = np.arange(len(obs))
    out = np.empty((perms.shape[0], len(obs)), dtype=np.int64)
    blocks = [np.nonzero(lineage == l)[0] for l in sorted(set(lineage))]
    for b in range(perms.shape[0]):
        for idx in blocks:
            tgt = idx[loc[idx] >= 0]                    # observed targets, original order
            src = perms[b, idx]                          # sources in target order
            src = src[loc[src] >= 0]                     # observed sources, permuted order
            out[b, loc[tgt]] = loc[src]
    return out


def run_panel(Y, lineage, x, perms, min_per_group=1, log=None):
    """Y: (L, C, G) with NaN for missing cells. lineage, x: (L,). perms: (B, L) within-lineage perms.
    For each compound, lines with an observed cell are used; lineages that lose one of the two
    groups among observed lines are dropped for that compound.
    Returns dict of arrays: T (C,), p (C,), beta (C, G), t (C, G), n, n1, n0, n_lineages."""
    L, C, G = Y.shape
    B = perms.shape[0]
    lineage = np.asarray(lineage); x = np.asarray(x, float)
    res = dict(T=np.full(C, np.nan), p=np.full(C, np.nan), beta=np.full((C, G), np.nan), t=np.full((C, G), np.nan),
               n=np.zeros(C, int), n1=np.zeros(C, int), n0=np.zeros(C, int), n_lineages=np.zeros(C, int),
               Tperm=np.full((C, B), np.nan))
    # group compounds by usable line set
    keys = {}
    for c in range(C):
        obs = ~np.isnan(Y[:, c, 0])
        for l in set(lineage):
            m = obs & (lineage == l)
            if (x[m] == 1).sum() < min_per_group or (x[m] == 0).sum() < min_per_group:
                obs = obs & ~(lineage == l)
        keys.setdefault(tuple(np.nonzero(obs)[0]), []).append(c)
    for obs, cs in keys.items():
        obs = np.array(obs)
        X, Z = design(lineage[obs], x[obs])
        a, h, M, H0 = fixed_parts(X, Z)
        Yo = np.concatenate([Y[obs, c, :] for c in cs], axis=1)        # (n, len(cs)*G)
        R0 = Yo - H0 @ Yo
        beta, t = tstats(R0, a, h, M)
        Tobs = (t.reshape(len(cs), G) ** 2).mean(1)
        P = perms if len(obs) == L else induced(perms, lineage, obs)
        cnt = np.zeros(len(cs)); Tp = np.empty((len(cs), B))
        for b in range(B):
            _, tb = tstats(R0[P[b]], a, h, M)
            Tb = (tb.reshape(len(cs), G) ** 2).mean(1)
            Tp[:, b] = Tb
            cnt += Tb >= Tobs - 1e-12 * np.abs(Tobs)
        for k, c in enumerate(cs):
            res["T"][c] = Tobs[k]; res["p"][c] = (1 + cnt[k]) / (B + 1)
            res["beta"][c] = beta[k * G:(k + 1) * G]; res["t"][c] = t[k * G:(k + 1) * G]
            res["n"][c] = len(obs); res["n1"][c] = int((x[obs] == 1).sum()); res["n0"][c] = int((x[obs] == 0).sum())
            res["n_lineages"][c] = len(set(lineage[obs])); res["Tperm"][c] = Tp[k]
        if log:
            log(f"group of {len(cs)} compounds, n = {len(obs)} lines done")
    return res
