"""Synthetic checks of the mosaic estimators: gradient, recovery of a known share, MH odds ratio."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.special import expit
from mosaic_models import nb_mix_nll, fit_nb_mix, fit_pois_mix, mh_or, depth_deciles, deep_zero

rng = np.random.default_rng(7)
n = 40000
d = np.exp(rng.normal(np.log(1800), 0.6, n)).round().clip(300, 30000)
b = np.full(n, 2e-5)
# truth: share negative 0.35 at 2,000 UMIs, rising slightly with depth
alpha, beta = np.log(0.35 / 0.65), 0.3
pi = expit(alpha + beta * np.log(d / 2000))
neg = rng.random(n) < pi
r, k = 3e-3, 4.0
mu = d * r
y_pos = rng.negative_binomial(k, k / (k + mu))
y = np.where(neg, rng.poisson(d * b), y_pos)
# gradient check
th = np.array([-0.4, 0.2, np.log(2.5e-3), np.log(3.0)])
z = np.log(d / 2000)
f0, g0 = nb_mix_nll(th, y, d, b, z)
gn = np.array([(nb_mix_nll(th + e, y, d, b, z)[0] - nb_mix_nll(th - e, y, d, b, z)[0]) / 2e-6 for e in np.eye(4) * 1e-6])
assert np.allclose(g0, gn, rtol=1e-4, atol=1e-3), (g0, gn)
fit = fit_nb_mix(y, d, b)
print('NB mix fit', {k_: round(v, 4) if isinstance(v, float) else v for k_, v in fit.items()})
assert abs(fit['pi_2000'] - 0.35) < 0.03 and abs(fit['beta'] - 0.3) < 0.15 and abs(fit['r_high'] / r - 1) < 0.05
assert fit['pi_2000_lo'] < 0.35 < fit['pi_2000_hi'] or abs(fit['pi_2000'] - 0.35) < 0.01
# homogeneous line: share near 0
y_h = rng.negative_binomial(k, k / (k + mu))
fh = fit_nb_mix(y_h, d, b)
print('homogeneous NB fit pi_2000', round(fh['pi_2000'], 4))
assert fh['pi_2000'] < 0.02
pm = fit_pois_mix(y, d, b)
print('Poisson mix', pm)
# Poisson mixture on Poisson data recovers the share
yp = np.where(rng.random(n) < 0.3, rng.poisson(d * b), rng.poisson(d * r))
pmp = fit_pois_mix(yp, d, b)
assert abs(pmp['pi'] - 0.3) < 0.02, pmp
# MH odds ratio: independent zeros given depth give OR near 1; a mosaic gives OR > 1
dec = depth_deciles(d)
lp, lq = d * 1.0e-3, d * 0.8e-3
zp = rng.poisson(lp) == 0; zq = rng.poisson(lq) == 0
o1 = mh_or(zq, zp, dec)
negm = rng.random(n) < 0.3
zp2 = np.where(negm, True, rng.poisson(lp) == 0); zq2 = np.where(negm, True, rng.poisson(lq) == 0)
o2 = mh_or(zq2, zp2, dec)
print('MH OR independent', o1, 'mosaic', o2)
assert 0.85 < o1['or_mh'] < 1.15 and o2['or_mh'] > 3
dz = deep_zero(y, d, 3000, r, (r, k), rng)
print('deep zero', dz)
print('ALL MOSAIC SYNTHETIC TESTS PASSED')
