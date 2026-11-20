"""Line-level permutation engine for task D (and task G).

Model per compound c and gene g (lines are the unit):
    y_ig = lineage FE + Z_i' gamma + T_i' beta + e_ig
The tested term is one column k of T. The statistic per gene is the HC3 t of
beta_k; the compound statistic is T_c = mean over genes of t^2.

Inference: Freedman-Lane. Fit the reduced model (without column k), permute
its residuals within lineage, add back the reduced fit, refit the full model.
Because the reduced fit lies in the column space of the full model, adding it
back changes neither beta_k nor the full-model residuals, so the engine works
with the permuted reduced residuals directly (exact identity, not an
approximation). Lineage fixed effects are removed by within-lineage centring
(Frisch-Waugh-Lovell); within-lineage permutations commute with that centring.
Leverages include the lineage part (1 / n_lineage).

One set of permutations of line indices is shared by all compounds: a matrix
of random keys U (B x N lines); within each lineage, lines are ordered by key.
For a compound that lacks some lines, the induced permutation of its available
lines in each lineage is the order of their keys, so compounds with all lines
use exactly the same permutation.

p = (1 + #{b: T_c^(b) >= T_c^obs}) / (B + 1).
"""
import numpy as np
try:
    import numexpr as _ne
except ImportError:  # optional speed-up
    _ne = None

EPS_LEV = 1e-8


def make_keys(n_lines, B, seed):
    rng = np.random.default_rng(seed)
    return rng.random((B, n_lines))


def induced_perm(keys, lin_codes):
    """keys: (B, n) random keys for the available lines (in local order).
    lin_codes: (n,) lineage code per available line.
    Returns perm (B, n) of local indices: permuted y_i = e[perm[b, i]]."""
    B, n = keys.shape
    perm = np.empty((B, n), dtype=np.int64)
    for l in np.unique(lin_codes):
        idx = np.flatnonzero(lin_codes == l)
        order = np.argsort(keys[:, idx], axis=1, kind="stable")
        perm[:, idx] = idx[order]
    return perm


def centre(M, lin_codes):
    M = np.asarray(M, dtype=float)
    out = M.copy()
    for l in np.unique(lin_codes):
        idx = lin_codes == l
        out[idx] -= M[idx].mean(axis=0, keepdims=True)
    return out


class Design:
    """Design for one line subset. X_other: (n, p') non-FE columns (covariates
    and all tested-model terms); k: index of tested column within X_other."""

    def __init__(self, lin_codes, X_other, k):
        self.lin = np.asarray(lin_codes)
        self.n = len(self.lin)
        Xd = centre(X_other, self.lin)
        self.Xd = Xd
        self.k = k
        self.p = Xd.shape[1]
        nl = {l: (self.lin == l).sum() for l in np.unique(self.lin)}
        self.n_fe = len(nl)
        self.df = self.n - self.n_fe - self.p
        XtX = Xd.T @ Xd
        cond = np.linalg.cond(XtX)
        if not np.isfinite(cond) or cond > 1e12:
            raise ValueError(f"singular design (cond={cond:.3g})")
        XtXi = np.linalg.inv(XtX)
        self.A = XtXi @ Xd.T  # (p, n)
        h_fe = np.array([1.0 / nl[l] for l in self.lin])
        h = h_fe + np.einsum("ij,jk,ik->i", Xd, XtXi, Xd)
        if np.any(h > 1 - EPS_LEV):
            raise ValueError("leverage 1 for some line (HC3 undefined)")
        self.h = h
        self.w = self.A[k] ** 2 / (1.0 - h) ** 2  # HC3 weights for term k
        self.K = (Xd * self.w[:, None]).T @ Xd  # (p, p)
        red = [j for j in range(self.p) if j != k]
        Xr = Xd[:, red]
        if Xr.shape[1] > 0:
            self.Pr = Xr @ np.linalg.solve(Xr.T @ Xr, Xr.T)
        else:
            self.Pr = np.zeros((self.n, self.n))
        self.xtxi_kk = XtXi[k, k]

    def reduced_resid(self, Y):
        Yd = centre(Y, self.lin)
        return Yd - self.Pr @ Yd

    def observed(self, Y):
        """Returns beta_k (G,), HC3 se (G,), t (G,), resid SD (G,)."""
        Yd = centre(Y, self.lin)
        beta = self.A @ Yd  # (p, G)
        e = Yd - self.Xd @ beta
        var = (self.w[:, None] * e ** 2).sum(axis=0)
        se = np.sqrt(var)
        sig = np.sqrt((e ** 2).sum(axis=0) / self.df)
        # convention: a gene with beta = 0 and se = 0 (e.g. identical values in
        # every line that carries weight) has t = 0, in observed and permuted data
        t = np.divide(beta[self.k], se, out=np.zeros_like(se), where=se > 0)
        return beta[self.k], se, t, sig

    def perm_T(self, E, perm, dtype=np.float64, chunk=1000, n_genes=None):
        """E: (n, G_total) reduced residuals (several compounds with this same
        design may be concatenated along axis 1). perm: (B, n).
        Returns tstat^2 summed per block of n_genes columns: (B, G_total/n_genes)."""
        # With y* = P e:  u = A y*  (all p coefficients),
        # V = sum_i w_i (y*_i - Xd_i u)^2 = sum_i w_i y*_i^2 - u'(2 X'W y* - K u)
        # and 2 X'W y* - K u = (2 WX - K A) y*, so both u and R = (2WX - KA) y*
        # are single matrix products of the permuted design with E.
        B = perm.shape[0]
        G = E.shape[1]
        if n_genes is None:
            n_genes = G
        nb = G // n_genes
        assert nb * n_genes == G
        inv = np.argsort(perm, axis=1)
        E_ = np.ascontiguousarray(E, dtype=dtype)
        E2 = np.ascontiguousarray(E ** 2, dtype=dtype)
        WX = (self.w[:, None] * self.Xd).T  # (p, n)
        Rm = 2.0 * WX - self.K @ self.A  # (p, n)
        p, k = self.p, self.k
        # row order: u_0..u_{p-1}, R_0..R_{p-1}, w
        M = np.vstack([self.A, Rm, self.w[None, :]]).astype(dtype)  # (2p+1, n)
        m = M.shape[0]
        out = np.empty((B, nb), dtype=np.float64)
        for s in range(0, B, chunk):
            iv = inv[s:s + chunk]
            b = iv.shape[0]
            Mb = M[:, iv].transpose(1, 0, 2)  # (b, m, n)
            UR = (Mb[:, :2 * p].reshape(b * 2 * p, self.n) @ E_).reshape(b, 2 * p, G)
            S1 = Mb[:, 2 * p] @ E2  # (b, G)
            if _ne is not None:
                loc = {"S1": S1, "uk": UR[:, k]}
                expr = "S1"
                for j in range(p):
                    loc[f"u{j}"] = UR[:, j]
                    loc[f"r{j}"] = UR[:, p + j]
                    expr += f" - u{j}*r{j}"
                t2 = _ne.evaluate(f"where(({expr}) > 0, uk*uk/({expr}), 0)", local_dict=loc)
            else:
                V = S1 - np.einsum("bjg,bjg->bg", UR[:, :p], UR[:, p:])
                t2 = np.divide(UR[:, k] ** 2, V, out=np.zeros_like(V), where=V > 0)
            out[s:s + b] = t2.reshape(b, nb, n_genes).mean(axis=2)
        return out


def brute_force_T(lin, X_other, k, Y, perm):
    """Reference implementation by explicit refits (slow; for tests only)."""
    n = len(lin)
    levels = np.unique(lin)
    FE = (lin[:, None] == levels[None, :]).astype(float)
    X = np.hstack([FE, X_other])
    kk = FE.shape[1] + k
    Xr = np.delete(X, kk, axis=1)
    br = np.linalg.lstsq(Xr, Y, rcond=None)[0]
    fit_r = Xr @ br
    er = Y - fit_r
    XtXi = np.linalg.inv(X.T @ X)
    H = X @ XtXi @ X.T
    h = np.diag(H)
    res = []
    for b in range(perm.shape[0]):
        Ys = fit_r + er[perm[b]]
        bet = XtXi @ X.T @ Ys
        e = Ys - X @ bet
        a = (XtXi @ X.T)[kk]
        var = ((a ** 2)[:, None] * e ** 2 / ((1 - h) ** 2)[:, None]).sum(axis=0)
        res.append(np.mean(bet[kk] ** 2 / var))
    return np.array(res)
