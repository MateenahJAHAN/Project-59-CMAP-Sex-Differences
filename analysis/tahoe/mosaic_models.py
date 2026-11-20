"""Estimators for the mosaic tests (a) to (e)."""
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln, digamma, expit
from scipy.stats import nbinom
from statsmodels.stats.contingency_tables import StratifiedTable

D0 = 2000.0
DEPTHS = [1000, 2000, 4000, 8000]


def pois_logpmf(y, lam):
    lam = np.maximum(lam, 1e-300)
    return y * np.log(lam) - lam - gammaln(y + 1)


def nb_logpmf(y, mu, k):
    return (gammaln(y + k) - gammaln(k) - gammaln(y + 1) + k * np.log(k / (k + mu)) + y * np.log(mu / (k + mu)))


# ---------------- (b) NB high component + Poisson low component at the plate background, logistic share
def nb_mix_nll(theta, y, d, b, z):
    a, be, lr, lk = theta
    r, k = np.exp(lr), np.exp(lk)
    eta = a + be * z
    log_pi = -np.logaddexp(0, -eta)
    log_1mpi = -np.logaddexp(0, eta)
    A = log_pi + pois_logpmf(y, d * b)
    mu = d * r
    B = log_1mpi + nb_logpmf(y, mu, k)
    ll = np.logaddexp(A, B)
    w = np.exp(A - ll)
    pi = expit(eta)
    g_a = w - pi
    g_b = (w - pi) * z
    g_lr = (1 - w) * k * (y - mu) / (k + mu)
    g_lk = (1 - w) * k * (digamma(y + k) - digamma(k) + np.log(k / (k + mu)) + (mu - y) / (k + mu))
    return -ll.sum(), -np.array([g_a.sum(), g_b.sum(), g_lr.sum(), g_lk.sum()])


def fit_nb_mix(y, d, b):
    y = np.asarray(y, float); d = np.asarray(d, float); b = np.asarray(b, float)
    z = np.log(d / D0)
    pos = y >= 1
    r0 = y[pos].sum() / d[pos].sum() if pos.any() else 1e-3
    bounds = [(-20, 20), (-20, 20), (np.log(1e-7), np.log(0.5)), (np.log(1e-3), np.log(1e4))]
    best = None
    for a0 in (-2.0, 0.0, 1.0):
        for k0 in (1.0, 10.0):
            x0 = np.array([a0, 0.0, np.log(r0), np.log(k0)])
            res = minimize(nb_mix_nll, x0, args=(y, d, b, z), jac=True, method='L-BFGS-B', bounds=bounds,
                           options=dict(maxiter=2000, ftol=1e-12, gtol=1e-7))
            if best is None or res.fun < best.fun:
                best = res
    a, be, lr, lk = best.x
    out = dict(alpha=a, beta=be, r_high=np.exp(lr), k_high=np.exp(lk), nll=best.fun, converged=bool(best.success),
               n_cells=len(y), at_bound=bool(any(np.isclose(best.x[i], bounds[i][j]) for i in range(4) for j in range(2))))
    for D in DEPTHS:
        out[f'pi_{D}'] = float(expit(a + be * np.log(D / D0)))
    # delta-method CI at 2,000 UMIs (pi = expit(alpha)), Hessian by central differences of the analytic gradient
    try:
        H = np.zeros((4, 4))
        for i in range(4):
            e = np.zeros(4); h = 1e-4 * max(1.0, abs(best.x[i])); e[i] = h
            H[:, i] = (nb_mix_nll(best.x + e, y, d, b, z)[1] - nb_mix_nll(best.x - e, y, d, b, z)[1]) / (2 * h)
        H = (H + H.T) / 2
        cov = np.linalg.inv(H)
        se = float(np.sqrt(cov[0, 0])) if cov[0, 0] > 0 else np.nan
    except np.linalg.LinAlgError:
        se = np.nan
    out['pi_2000_lo'] = float(expit(a - 1.96 * se)) if np.isfinite(se) else np.nan
    out['pi_2000_hi'] = float(expit(a + 1.96 * se)) if np.isfinite(se) else np.nan
    return out


# ---------------- (c) Poisson mixture with the plate background, constant share, EM
def fit_pois_mix(y, d, b, tol=1e-9, max_iter=20000):
    y = np.asarray(y, float); d = np.asarray(d, float); b = np.asarray(b, float)
    pos = y >= 1
    r = y[pos].sum() / d[pos].sum() if pos.any() else 1e-3
    pi = 0.5
    lp_low = pois_logpmf(y, d * b)
    prev = -np.inf
    for it in range(1, max_iter + 1):
        A = np.log(pi) + lp_low
        B = np.log1p(-pi) + pois_logpmf(y, d * r)
        ll = np.logaddexp(A, B)
        tot = ll.sum()
        w = np.exp(A - ll)
        pi = float(np.clip(w.mean(), 1e-12, 1 - 1e-12))
        r = float(((1 - w) * y).sum() / ((1 - w) * d).sum())
        if abs(tot - prev) < tol * max(1.0, abs(tot)):
            break
        prev = tot
    return dict(pi=pi, r_high=r, ll=tot, iterations=it, converged=it < max_iter)


# ---------------- (a) deep-cell zero fraction
def deep_zero(y, d, thr, r_pos, nbpar, rng, B=2000):
    m = d >= thr
    n = int(m.sum())
    if n == 0:
        return dict(n_deep=0, zero_frac=np.nan, ci_lo=np.nan, ci_hi=np.nan, exp_zero_poisson=np.nan, exp_zero_nb=np.nan)
    z = (y[m] == 0)
    p = z.mean()
    # nonparametric bootstrap over cells; for a 0/1 statistic the resampled count of zeros is Binomial(n, p)
    boot = rng.binomial(n, p, size=B) / n
    lo, hi = np.percentile(boot, [2.5, 97.5])
    exp_p = float(np.exp(-d[m] * r_pos).mean())
    exp_nb = float(np.exp(nb_logpmf(0, d[m] * nbpar[0], nbpar[1])).mean()) if nbpar is not None else np.nan
    return dict(n_deep=n, n_zero=int(z.sum()), zero_frac=float(p), ci_lo=float(lo), ci_hi=float(hi),
                exp_zero_poisson=exp_p, exp_zero_nb=exp_nb)


# ---------------- (e) Mantel-Haenszel odds ratio of zero co-occurrence by depth decile
def mh_or(zq, zp, dec):
    tabs = []
    for g in np.unique(dec):
        m = dec == g
        a = int((zq[m] & zp[m]).sum()); b_ = int((zq[m] & ~zp[m]).sum())
        c = int((~zq[m] & zp[m]).sum()); dd = int((~zq[m] & ~zp[m]).sum())
        tabs.append(np.array([[a, b_], [c, dd]], float))
    num = sum(t[0, 0] * t[1, 1] / t.sum() for t in tabs)
    den = sum(t[0, 1] * t[1, 0] / t.sum() for t in tabs)
    if den == 0 or num == 0:
        return dict(or_mh=np.inf if den == 0 else 0.0, or_lo=np.nan, or_hi=np.nan, n_strata=len(tabs),
                    n_both_zero=int(sum(t[0, 0] for t in tabs)))
    st = StratifiedTable(tabs)
    lo, hi = st.oddsratio_pooled_confint()
    orr = st.oddsratio_pooled
    assert abs(orr - num / den) < 1e-8 * max(1, orr)
    return dict(or_mh=float(orr), or_lo=float(lo), or_hi=float(hi), n_strata=len(tabs),
                n_both_zero=int(sum(t[0, 0] for t in tabs)))


def depth_deciles(d):
    return pd.qcut(pd.Series(d), 10, labels=False, duplicates='drop').to_numpy()
