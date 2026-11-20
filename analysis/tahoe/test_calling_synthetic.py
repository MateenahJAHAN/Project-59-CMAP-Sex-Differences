"""Synthetic tests of the calling functions (no real data)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from common import kde_cutoff, rule_state, arm_state, Y8, scores_from_sums
from calling import aggregate, call, filtered

rng = np.random.default_rng(1)
# 1. kde_cutoff on a clear bimodal sample
v = np.r_[rng.normal(1, 0.3, 300), rng.normal(6, 0.5, 400)]
c = kde_cutoff(v)
assert 2 < c['cutoff'] < 5, c
assert abs(c['mode_low'] - 1) < 0.3 and abs(c['mode_high'] - 6) < 0.3, c
# trimodal: two tallest modes kept
v3 = np.r_[rng.normal(0, 0.2, 500), rng.normal(3, 0.2, 50), rng.normal(6, 0.2, 500)]
c3 = kde_cutoff(v3)
assert c3['n_modes'] >= 3 and abs(c3['mode_low']) < 0.3 and abs(c3['mode_high'] - 6) < 0.3, c3
# unimodal: no cut-off
c1 = kde_cutoff(rng.normal(0, 1, 500))
assert not np.isfinite(c1['cutoff']), c1
# 2. rule table truth table
exp = {('Male', 0, 0): 'Y-lost', ('Male', 0, 1): 'XY', ('Male', 1, 0): 'discordant', ('Male', 1, 1): 'XY',
       ('Female', 0, 0): 'XIST-lost', ('Female', 0, 1): 'discordant', ('Female', 1, 0): 'XX-like',
       ('Female', 1, 1): 'discordant', ('Unknown', 0, 0): 'ambiguous', ('Unknown', 0, 1): 'XY',
       ('Unknown', 1, 0): 'XX-like', ('Unknown', 1, 1): 'XY'}
for (s, x, y), e in exp.items():
    assert rule_state(s, bool(x), bool(y)) == e, (s, x, y, rule_state(s, bool(x), bool(y)), e)
assert arm_state(True, True) == 'Y intact' and arm_state(True, False) == 'Yq lost'
assert arm_state(False, True) == 'Yp lost' and arm_state(False, False) == 'Y lost'
# 3. full call on synthetic cells: 20 XY, 10 Y-lost, 15 XX-like, 5 XIST-lost; foreign cells planted
rows = []
truth = {}
for i in range(50):
    line = f'L{i:02d}'
    kind = 'XY' if i < 20 else 'Ylost' if i < 30 else 'XX' if i < 45 else 'XISTlost'
    sex = 'Male' if kind in ('XY', 'Ylost') else 'Female'
    truth[line] = (sex, {'XY': 'XY', 'Ylost': 'Y-lost', 'XX': 'XX-like', 'XISTlost': 'XIST-lost'}[kind])
    for p in range(1, 15):
        n = 300
        depth = rng.integers(800, 6000, n)
        yrate = 2e-3 if kind == 'XY' else 2e-5
        xrate = 3e-3 if kind == 'XX' else 1e-5
        foreign = rng.random(n) < 0.02
        yr = np.where(foreign & (kind != 'XY'), 2e-3, yrate)
        xr = np.where(foreign & (kind != 'XX'), 3e-3, xrate)
        df = pd.DataFrame({'cell_line_id': line, 'plate': f'plate{p}', 'total_umi': depth,
                           'XIST': rng.poisson(depth * xr)})
        for g in Y8:
            df[g] = rng.poisson(depth * yr / 8)
        rows.append(df)
cells = pd.concat(rows, ignore_index=True)
cells['y8'] = cells[Y8].sum(axis=1)
cells['stratum'] = 2 * (cells['y8'] >= 3) + (cells['XIST'] >= 3)
ann = pd.DataFrame({'cell_line_id': list(truth), 'cell_name': list(truth), 'organ': 'x',
                    'annotation': [truth[l][0] for l in truth], 'annotation_source': 'synthetic'})
agg = aggregate(cells)
assert agg['n_cells'].sum() == len(cells)
res = call(ann, 'lineplate', agg)
calls = res['calls'].set_index('cell_line_id')
bad = [(l, calls.loc[l, 'state'], truth[l][1]) for l in truth if calls.loc[l, 'state'] != truth[l][1]]
assert not bad, bad
print(res['hist'][['iteration', 'cut_XIST', 'cut_y', 'n_Yneg', 'n_XISTneg', 'cells_kept', 'converged']])
print(calls['confidence'].value_counts().to_dict())
# filter removed exactly the planted foreign-marker cells of negative lines
kept = int(res['filtered']['n_cells'].sum())
yneg = {l for l in truth if truth[l][1] != 'XY'}
xneg = {l for l in truth if truth[l][1] != 'XX-like'}
drop = ((cells.cell_line_id.isin(yneg) & (cells.y8 >= 3)) | (cells.cell_line_id.isin(xneg) & (cells.XIST >= 3))).sum()
assert kept == len(cells) - drop, (kept, len(cells), drop)
# scores: compare against a direct computation for one line x plate
sub = cells[(cells.cell_line_id == 'L00') & (cells.plate == 'plate3') & ~((cells.XIST >= 3))]
T = sub.total_umi.sum()
xs = np.log2(1e6 * sub.XIST.sum() / T + 1)
ys = np.mean([np.log2(1e6 * sub[g].sum() / T + 1) for g in Y8])
lp = res['lineplate'].set_index(['cell_line_id', 'plate'])
assert abs(lp.loc[('L00', 'plate3'), 'XIST'] - xs) < 1e-9 and abs(lp.loc[('L00', 'plate3'), 'y_score'] - ys) < 1e-9
res2 = call(ann, 'line', agg)
assert (res2['calls'].set_index('cell_line_id')['state'] == calls['state']).all()
print('ALL SYNTHETIC TESTS PASSED')
