"""Shared constants and functions for the Part A local caller (state-caller-local).

Written from the frozen rules (PartA_frozen_rules.md, commit 7cac9a4) and Task A's rule table only.
"""
import os
import time
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

DMSO = '$P59_DATA/tahoe/dmso'
META = '$P59_DATA/tahoe/index/meta'
# output folder; P59_PARTA_OUT allows a rerun into a fresh folder (added after the 09:21 UTC run, no change in computation)
OUT = os.environ.get('P59_PARTA_OUT', '$P59_DATA/tahoe/partA_local')
INPUTS = '$P59_ANALYSIS/inputs'
VEHICLE = 'DMSO_TF'
SEED = 20261115

Y8 = ['RPS4Y1', 'DDX3Y', 'KDM5D', 'UTY', 'EIF1AY', 'ZFY', 'USP9Y', 'NLGN4Y']
YP = ['RPS4Y1', 'ZFY']
YQ = ['DDX3Y', 'USP9Y', 'UTY', 'NLGN4Y', 'KDM5D', 'EIF1AY']
X8 = ['RPS4X', 'DDX3X', 'KDM5C', 'KDM6A', 'EIF1AX', 'ZFX', 'USP9X', 'NLGN4X']
OTHER = ['XIST', 'TSIX', 'TXLNGY', 'PRKY', 'ACTB', 'GAPDH']
MARKERS = ['XIST', 'TSIX'] + Y8 + ['TXLNGY', 'PRKY'] + X8 + ['ACTB', 'GAPDH']
assert set(YP) | set(YQ) == set(Y8) and not set(YP) & set(YQ)

STATES = ['XY', 'XX-like', 'Y-lost', 'XIST-lost', 'discordant', 'ambiguous']
MIN_CELLS_CALL = 200
MIN_CELLS_MOSAIC = 1000
FOREIGN_UMI = 3


def log(msg, fh=None):
    s = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()) + ' ' + msg
    print(s, flush=True)
    if fh is not None:
        fh.write(s + '\n')
        fh.flush()


def plate_num(p):
    s = str(p)
    assert s.startswith('plate'), s
    return int(s[5:])


def load_cells(columns=None):
    """Load cells_markers.parquet and stop unless every row is a vehicle (DMSO_TF) cell."""
    cells = pd.read_parquet(DMSO + '/cells_markers.parquet', columns=columns)
    if 'drug' in cells.columns:
        bad = (cells['drug'] != VEHICLE).sum()
        if bad:
            raise SystemExit(f'STOP: {bad} rows of cells_markers have drug != {VEHICLE}')
    return cells


def log2cpm1(counts, total):
    return np.log2(1e6 * np.asarray(counts, float) / np.asarray(total, float) + 1.0)


def scores_from_sums(df):
    """df has columns umi and the marker sums; returns XIST, y_score, yp_score, yq_score."""
    out = pd.DataFrame(index=df.index)
    T = df['umi'].astype(float)
    out['XIST'] = log2cpm1(df['XIST'], T)
    out['y_score'] = np.mean([log2cpm1(df[g], T) for g in Y8], axis=0)
    out['yp_score'] = np.mean([log2cpm1(df[g], T) for g in YP], axis=0)
    out['yq_score'] = np.mean([log2cpm1(df[g], T) for g in YQ], axis=0)
    return out


def kde_cutoff(values, n_grid=1024):
    """Scott-bandwidth Gaussian KDE on a 1,024-point grid from min to max of the values.
    Cut-off = grid point of lowest density between the two tallest modes."""
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    kde = gaussian_kde(v, bw_method='scott')
    grid = np.linspace(v.min(), v.max(), n_grid)
    d = kde(grid)
    is_max = np.zeros(n_grid, bool)
    is_max[1:-1] = (d[1:-1] > d[:-2]) & (d[1:-1] >= d[2:])
    is_max[0] = d[0] > d[1]
    is_max[-1] = d[-1] > d[-2]
    idx = np.where(is_max)[0]
    res = dict(n_points=len(v), bandwidth_sd=float(np.sqrt(kde.covariance[0, 0])),
               scott_factor=float(kde.factor), grid_min=float(grid[0]), grid_max=float(grid[-1]),
               n_modes=int(len(idx)))
    if len(idx) < 2:
        res.update(cutoff=np.nan, mode_low=np.nan, mode_high=np.nan, dens_mode_low=np.nan,
                   dens_mode_high=np.nan, dens_cutoff=np.nan)
        return res
    top2 = np.sort(idx[np.argsort(d[idx])[-2:]])
    lo, hi = int(top2[0]), int(top2[1])
    ci = lo + int(np.argmin(d[lo:hi + 1]))
    res.update(cutoff=float(grid[ci]), mode_low=float(grid[lo]), mode_high=float(grid[hi]),
               dens_mode_low=float(d[lo]), dens_mode_high=float(d[hi]), dens_cutoff=float(d[ci]),
               all_modes=';'.join(f'{grid[i]:.4f}:{d[i]:.4f}' for i in idx))
    return res


def rule_state(sex, x_hi, y_hi):
    """Task A rule table, evaluated in a fixed order (see PLAN.md)."""
    if sex == 'Female' and y_hi:
        return 'discordant'
    if sex == 'Male' and x_hi and not y_hi:
        return 'discordant'
    if y_hi:
        return 'XY'
    if x_hi:
        return 'XX-like'
    if sex == 'Male':
        return 'Y-lost'
    if sex == 'Female':
        return 'XIST-lost'
    return 'ambiguous'


def arm_state(yp_hi, yq_hi):
    if yp_hi and yq_hi:
        return 'Y intact'
    if yp_hi and not yq_hi:
        return 'Yq lost'
    if yq_hi and not yp_hi:
        return 'Yp lost'
    return 'Y lost'


def norm_name(s):
    return ''.join(ch for ch in str(s).lower() if ch.isalnum())
