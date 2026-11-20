"""Line-level cut-off sensitivity (frozen rule), plus two exploratory checks of implementation choices:
(1) calls without the foreign-cell filter, (2) cut-offs on a KDE grid extended by 3 bandwidths."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from common import OUT, Y8, log, load_cells, kde_cutoff, rule_state
from calling import aggregate, call, filtered, scores_at

T0 = time.time()
LOG = open(OUT + '/logs/03_sensitivity.log', 'a')
log('start 03_sensitivity', LOG)
cells = load_cells(['drug', 'cell_line_id', 'plate', 'stratum', 'total_umi', 'XIST'] + Y8)
ann = pd.read_csv(OUT + '/partA_local_annotation.tsv', sep='\t')
prim = pd.read_csv(OUT + '/partA_local_states.tsv', sep='\t')
agg = aggregate(cells)
assert int(agg['n_cells'].sum()) == len(cells)

# ---- line-level cut-offs (frozen sensitivity)
res = call(ann, 'line', agg)
c = res['calls']
ll = pd.DataFrame({'cell_name': c['cell_name'], 'cellosaurus_id': c['cell_line_id'], 'annotation': c['annotation'],
                   'n_vehicle_cells': c['n_cells'].astype(int), 'XIST': c['XIST'].round(4),
                   'y_score': c['y_score'].round(4), 'yp_score': c['yp_score'].round(4),
                   'yq_score': c['yq_score'].round(4), 'state': c['state'], 'arm_state': c['arm_state'],
                   'confidence': c['confidence'], 'plates_agree': c['plates_agree'], 'grid_stable': c['grid_stable'],
                   'n_plates_agree': c['n_plates_agree'], 'grid_states': c['grid_states'],
                   'plate_states': c['plate_states']})
ll = ll.sort_values('cell_name', key=lambda s: s.str.lower()).reset_index(drop=True)
ll.to_csv(OUT + '/partA_local_states_linelevel.tsv', sep='\t', index=False)
res['cuts'].to_csv(OUT + '/partA_local_cutoffs_linelevel.tsv', sep='\t', index=False)
res['hist'].to_csv(OUT + '/partA_local_filter_iterations_linelevel.tsv', sep='\t', index=False)
log('line-level cut-offs: ' + '; '.join(f"{r.score} {r.cutoff:.4f}" for r in res['cuts'].itertuples()), LOG)

# ---- exploratory: no foreign-cell filter, line x plate cut-offs
a = filtered(agg, {}, {})
lps, lns = scores_at(a)
nf_cut = {s: kde_cutoff(lps[s])['cutoff'] for s in ['XIST', 'y_score']}
lns = lns.reset_index().merge(ann, on='cell_line_id', validate='1:1')
lns['state_nofilter'] = [rule_state(s, x >= nf_cut['XIST'], y >= nf_cut['y_score'])
                         for s, x, y in zip(lns.annotation, lns.XIST, lns.y_score)]
nf = lns[['cell_line_id', 'n_cells', 'XIST', 'y_score', 'yp_score', 'yq_score', 'state_nofilter']].rename(
    columns={'cell_line_id': 'cellosaurus_id', 'n_cells': 'n_cells_nofilter', 'XIST': 'XIST_nofilter',
             'y_score': 'y_score_nofilter', 'yp_score': 'yp_score_nofilter', 'yq_score': 'yq_score_nofilter'})

# ---- comparison table
cmp = prim[['cell_name', 'cellosaurus_id', 'annotation', 'n_vehicle_cells', 'XIST', 'y_score', 'state', 'confidence',
            'arm_state']].merge(
    ll[['cellosaurus_id', 'n_vehicle_cells', 'XIST', 'y_score', 'state', 'confidence', 'arm_state']],
    on='cellosaurus_id', suffixes=('_primary', '_linelevel'), validate='1:1').merge(nf, on='cellosaurus_id', validate='1:1')
assert len(cmp) == 50
cmp['state_differs_linelevel'] = cmp['state_primary'] != cmp['state_linelevel']
cmp['state_differs_nofilter'] = cmp['state_primary'] != cmp['state_nofilter']
cmp['cut_XIST_nofilter'] = nf_cut['XIST']
cmp['cut_y_nofilter'] = nf_cut['y_score']
cmp.to_csv(OUT + '/partA_local_sensitivity_compare.tsv', sep='\t', index=False, float_format='%.4f')
log(f"line-level: {int(cmp.state_differs_linelevel.sum())} state changes: "
    + '; '.join(f'{r.cell_name} {r.state_primary}->{r.state_linelevel}' for r in cmp[cmp.state_differs_linelevel].itertuples()), LOG)
log(f"no filter (exploratory; cut XIST {nf_cut['XIST']:.4f}, y {nf_cut['y_score']:.4f}): "
    f"{int(cmp.state_differs_nofilter.sum())} state changes: "
    + '; '.join(f'{r.cell_name} {r.state_primary}->{r.state_nofilter}' for r in cmp[cmp.state_differs_nofilter].itertuples()), LOG)

# ---- exploratory: KDE grid extended by 3 bandwidths (check of the min-max grid choice)
fl = pd.read_csv(OUT + '/partA_local_lineplate_scores.tsv', sep='\t')
rows = []
for lev, df in [('lineplate', fl), ('line', None)]:
    for s in ['XIST', 'y_score', 'yp_score', 'yq_score']:
        if lev == 'lineplate':
            v = df[s].to_numpy()
            base = kde_cutoff(v)['cutoff']
        else:
            v = c[s].to_numpy()
            base = res['cuts'].set_index('score').loc[s, 'cutoff']
            assert abs(kde_cutoff(v)['cutoff'] - base) < 1e-12
        kde = gaussian_kde(v, bw_method='scott')
        bw = np.sqrt(kde.covariance[0, 0])
        g = np.linspace(v.min() - 3 * bw, v.max() + 3 * bw, 1024)
        d = kde(g)
        im = np.where((d[1:-1] > d[:-2]) & (d[1:-1] >= d[2:]))[0] + 1
        t2 = np.sort(im[np.argsort(d[im])[-2:]])
        cut = g[t2[0] + np.argmin(d[t2[0]:t2[1] + 1])]
        rows.append(dict(level=lev, score=s, cutoff_minmax_grid=base, cutoff_extended_grid=float(cut),
                         difference=float(cut - base), grid_step_minmax=float((v.max() - v.min()) / 1023)))
pd.DataFrame(rows).to_csv(OUT + '/partA_local_kde_grid_check.tsv', sep='\t', index=False, float_format='%.5f')
log('grid check: ' + '; '.join(f"{r['level']} {r['score']} diff {r['difference']:.4f}" for r in rows), LOG)
log(f'done 03_sensitivity in {time.time()-T0:.1f} s', LOG)
