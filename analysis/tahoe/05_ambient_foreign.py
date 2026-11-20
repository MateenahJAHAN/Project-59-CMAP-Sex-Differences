"""Step 5: ambient Y and XIST per plate (and per well) and foreign-marker shares per line."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.stats import poisson, pearsonr, spearmanr
from common import OUT, Y8, log, load_cells, plate_num

T0 = time.time()
LOG = open(OUT + '/logs/05_ambient_foreign.log', 'a')
log('start 05_ambient_foreign', LOG)
cells = load_cells(['drug', 'cell_line_id', 'cell_name', 'plate', 'sample', 'total_umi', 'XIST', 'y8'] + Y8)
S = pd.read_csv(OUT + '/partA_local_states.tsv', sep='\t')
st = S.set_index('cellosaurus_id')
yneg = set(st.index[~st.y_high])
xneg = set(st.index[~st.xist_high])
y_lost_both_arms = set(st.index[(~st.y_high) & (st.arm_state == 'Y lost')])
partial = sorted(yneg - y_lost_both_arms)
log(f'Y-negative lines {len(yneg)}; with both arms low {len(y_lost_both_arms)}; partial-Y among them: '
    + ', '.join(st.loc[partial, 'cell_name']), LOG)
cells['plate_n'] = cells['plate'].map(plate_num)

rows = []
def rate(df, col, unit, setname, filt):
    d = df if not filt else df[df[col] < 3]
    g = d.groupby(unit, observed=True).agg(marker_umi=(col, 'sum'), umi=('total_umi', 'sum'), n_cells=('total_umi', 'size'),
                                           n_lines=('cell_line_id', 'nunique'))
    g['per_1000'] = 1000 * g['marker_umi'] / g['umi']
    g = g.reset_index().rename(columns={unit: 'unit'})
    g['unit_type'] = unit
    g['marker'] = 'Y (8 genes)' if col == 'y8' else 'XIST'
    g['line_set'] = setname
    g['cells'] = 'marker < 3 (foreign cells removed)' if filt else 'all vehicle cells'
    rows.append(g)

for col, sets in [('y8', [('Y-negative lines (y_score below cut-off)', yneg),
                          ('Y-negative lines with both arms low (partial-Y lines removed)', y_lost_both_arms)]),
                  ('XIST', [('XIST-negative lines', xneg)])]:
    for setname, lines in sets:
        sub = cells[cells.cell_line_id.isin(lines)]
        for unit in ['plate', 'sample']:
            for filt in [False, True]:
                rate(sub, col, unit, setname, filt)
amb = pd.concat(rows, ignore_index=True)
amb = amb[['marker', 'line_set', 'cells', 'unit_type', 'unit', 'n_lines', 'n_cells', 'umi', 'marker_umi', 'per_1000']]
amb.to_csv(OUT + '/partA_local_ambient.tsv', sep='\t', index=False, float_format='%.6g')

# summary of plate rates
pl = amb[amb.unit_type == 'plate']
summ = pl.groupby(['marker', 'line_set', 'cells'])['per_1000'].agg(['min', 'median', 'max']).reset_index()
summ.to_csv(OUT + '/partA_local_ambient_summary.tsv', sep='\t', index=False, float_format='%.6g')
for r in summ.itertuples():
    log(f'ambient {r.marker} | {r.line_set} | {r.cells}: plate per 1,000 UMIs min {r.min:.4f} median {r.median:.4f} max {r.max:.4f}', LOG)

# measured plate background for the mixtures (marker < 3 cells): both Y versions saved
bg = pl[pl.cells.str.startswith('marker < 3')].copy()
bg['per_umi'] = bg['per_1000'] / 1000
bg.to_csv(OUT + '/partA_local_plate_background.tsv', sep='\t', index=False, float_format='%.6g')

# ---------------- foreign-marker shares per line
bmap_y = bg[bg.line_set.str.startswith('Y-negative lines with both')].set_index('unit')['per_umi']
bmap_y_all = bg[bg.line_set.str.startswith('Y-negative lines (y_score')].set_index('unit')['per_umi']
bmap_x = bg[bg.marker == 'XIST'].set_index('unit')['per_umi']
fr = []
for l, d in cells.groupby('cell_line_id', observed=True):
    n = len(d)
    rec = dict(cell_name=st.loc[l, 'cell_name'], cellosaurus_id=l, state=st.loc[l, 'state'], arm_state=st.loc[l, 'arm_state'],
               n_vehicle_cells_before_filter=n, log10_cells=np.log10(n),
               y_negative=l in yneg, xist_negative=l in xneg, partial_y=l in partial)
    yf = d['y8'] >= 3
    xf = d['XIST'] >= 3
    if l in yneg:
        lam = d['total_umi'].to_numpy() * d['plate'].map(bmap_y).to_numpy()
        lam_all = d['total_umi'].to_numpy() * d['plate'].map(bmap_y_all).to_numpy()
        rec.update(n_foreign_Y=int(yf.sum()), share_foreign_Y=yf.mean(),
                   expected_share_Y=float(poisson.sf(2, lam).mean()),
                   expected_share_Y_bg_all_yneg=float(poisson.sf(2, lam_all).mean()))
        rec['obs_over_exp_Y'] = rec['share_foreign_Y'] / rec['expected_share_Y']
    if l in xneg:
        lam = d['total_umi'].to_numpy() * d['plate'].map(bmap_x).to_numpy()
        rec.update(n_foreign_XIST=int(xf.sum()), share_foreign_XIST=xf.mean(),
                   expected_share_XIST=float(poisson.sf(2, lam).mean()))
        rec['obs_over_exp_XIST'] = rec['share_foreign_XIST'] / rec['expected_share_XIST']
    anyf = ((yf & (l in yneg)) | (xf & (l in xneg)))
    rec['n_removed_by_filter'] = int(anyf.sum())
    rec['share_removed_by_filter'] = anyf.mean()
    fr.append(rec)
fr = pd.DataFrame(fr).sort_values('n_vehicle_cells_before_filter').reset_index(drop=True)
chk = S.set_index('cellosaurus_id')
assert (fr.set_index('cellosaurus_id')['n_removed_by_filter'] ==
        (chk['n_vehicle_cells_before_filter'] - chk['n_vehicle_cells']).loc[fr.cellosaurus_id]).all()
fr.to_csv(OUT + '/partA_local_foreign_cells.tsv', sep='\t', index=False, float_format='%.6g')

cor = []
def corr(name, x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[ok], np.asarray(y)[ok]
    pr = pearsonr(x, y); sr = spearmanr(x, y)
    cor.append(dict(comparison=name, n_lines=len(x), pearson_r=pr.statistic, pearson_p=pr.pvalue,
                    spearman_rho=sr.statistic, spearman_p=sr.pvalue))

fy = fr[fr.y_negative]
corr('share y8>=3 vs log10 cells, all Y-negative lines', fy.log10_cells, fy.share_foreign_Y)
fyb = fy[~fy.partial_y]
corr('share y8>=3 vs log10 cells, Y-negative lines with both arms low', fyb.log10_cells, fyb.share_foreign_Y)
fx = fr[fr.xist_negative]
corr('share XIST>=3 vs log10 cells, XIST-negative lines', fx.log10_cells, fx.share_foreign_XIST)
corr('share removed by the filter vs log10 cells, all 50 lines', fr.log10_cells, fr.share_removed_by_filter)
fnp = fr[~fr.partial_y]
corr('share removed by the filter vs log10 cells, partial-Y lines removed', fnp.log10_cells, fnp.share_removed_by_filter)
corr('log10 share y8>=3 vs log10 cells, Y-negative lines with both arms low', fyb.log10_cells,
     np.log10(fyb.share_foreign_Y.where(fyb.share_foreign_Y > 0)))
cor = pd.DataFrame(cor)
cor.to_csv(OUT + '/partA_local_foreign_correlation.tsv', sep='\t', index=False, float_format='%.6g')
for r in cor.itertuples():
    log(f'{r.comparison}: n {r.n_lines}, Pearson r {r.pearson_r:.3f} (p {r.pearson_p:.2g}), Spearman {r.spearman_rho:.3f}', LOG)
log('obs/exp Y foreign share (both-arms-low lines): median %.1f, range %.1f-%.1f' % (
    fyb.obs_over_exp_Y.median(), fyb.obs_over_exp_Y.min(), fyb.obs_over_exp_Y.max()), LOG)
log('obs/exp XIST foreign share: median %.1f, range %.1f-%.1f' % (
    fx.obs_over_exp_XIST.median(), fx.obs_over_exp_XIST.min(), fx.obs_over_exp_XIST.max()), LOG)
log(f'done 05_ambient_foreign in {time.time()-T0:.1f} s', LOG)
