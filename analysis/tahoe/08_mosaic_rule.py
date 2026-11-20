"""Mosaic rule (fixed in the frozen rules): (b) pi(2,000) in [0.15, 0.85]; (d) |odd - even| <= 0.10;
(e) OR > 2 x control (Y only); (f) within 0.01 of the split-half and at least 0.10 above any other line."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from common import OUT, log

LOG = open(OUT + '/logs/08_mosaic_rule.log', 'a')
log('start 08_mosaic_rule', LOG)
R = pd.read_csv(OUT + '/partA_local_mosaic_abcde.tsv', sep='\t')
I = pd.read_csv(OUT + '/partA_local_identity.tsv', sep='\t')
S = pd.read_csv(OUT + '/partA_local_states.tsv', sep='\t')
M = R.merge(I[['cellosaurus_id', 'marker', 'n_deep_negative', 'n_deep_positive', 'n_pool', 'r_neg_pos', 'r_split_half',
               'r_neg_best_other', 'best_other_line', 'f_gap_to_split_half', 'f_margin_over_best_other', 'pass_f_split',
               'pass_f_other', 'pass_f', 'nc_neg_own_share', 'nc_pos_own_share']],
            on=['cellosaurus_id', 'marker'], how='left', validate='1:1')
assert len(M) == len(R)
verdict, failed = [], []
for r in M.itertuples():
    fails = []
    if not r.pass_b:
        fails.append('b')
    if not r.pass_d:
        fails.append('d')
    if r.marker == 'Y' and not (r.pass_e == 1.0):
        fails.append('e')
    has_f = isinstance(r.pass_f, (bool, np.bool_)) or (r.pass_f in (0, 1, True, False) and not pd.isna(r.pass_f))
    if has_f and not bool(r.pass_f):
        fails.append('f')
    if fails:
        v = 'not mosaic (fails ' + ', '.join(fails) + ')'
    elif not has_f:
        v = 'candidate: passes b, d' + (', e' if r.marker == 'Y' else '') + '; (f) pending (needs full rows)'
    else:
        v = 'mosaic (passes every part)'
    verdict.append(v); failed.append(';'.join(fails))
M['mosaic_verdict'] = verdict
M['failed_parts'] = failed
M['negative_control'] = M.cell_name.isin(['NCI-H460', 'LoVo']) & (M.marker == 'Y')
keep = ['cell_name', 'cellosaurus_id', 'marker', 'scope', 'n_cells', 'negative_control',
        'a3000_n_deep', 'a3000_zero_frac', 'a3000_ci_lo', 'a3000_ci_hi', 'a3000_exp_zero_poisson', 'a3000_exp_zero_nb',
        'a8000_n_deep', 'a8000_zero_frac', 'a8000_ci_lo', 'a8000_ci_hi', 'a8000_exp_zero_poisson', 'a8000_exp_zero_nb',
        'b_pi_1000', 'b_pi_2000', 'b_pi_2000_lo', 'b_pi_2000_hi', 'b_pi_4000', 'b_pi_8000', 'b_r_high', 'b_k_high',
        'c_pi', 'd_odd_pi_2000', 'd_even_pi_2000', 'd_diff_pi_2000', 'e_or_y', 'e_or_y_lo', 'e_or_y_hi', 'e_or_ctrl',
        'e_or_ctrl_lo', 'e_or_ctrl_hi', 'e_or_ctrl5050', 'e_ratio', 'r_neg_pos', 'r_split_half', 'r_neg_best_other',
        'best_other_line', 'f_gap_to_split_half', 'f_margin_over_best_other', 'nc_neg_own_share', 'nc_pos_own_share',
        'share_range_lo', 'share_range_hi', 'share_range_all_lo', 'share_range_all_hi',
        'pass_b', 'pass_d', 'pass_e', 'pass_f', 'mosaic_verdict', 'failed_parts']
M[keep].to_csv(OUT + '/partA_local_mosaic.tsv', sep='\t', index=False, float_format='%.6g')
for r in M.itertuples():
    log(f'{r.cell_name} {r.marker} [{r.scope}]: share range {r.share_range_lo:.3f}-{r.share_range_hi:.3f}; {r.mosaic_verdict}', LOG)
nc = M[M.negative_control]
nc_ok = not nc.mosaic_verdict.str.startswith('mosaic').any() and not nc.mosaic_verdict.str.startswith('candidate').any()
log(f'negative controls (NCI-H460, LoVo) fail the rule: {nc_ok}', LOG)
pend = M[M.mosaic_verdict.str.startswith('candidate')]
pd.DataFrame([dict(cell_name=r.cell_name, cellosaurus_id=r.cellosaurus_id, marker=r.marker,
                   needed='full vehicle-cell rows (genes, expressions) for test (f); not among the six extracted lines')
              for r in pend.itertuples()], columns=['cell_name', 'cellosaurus_id', 'marker', 'needed']).to_csv(
    OUT + '/partA_local_pending.tsv', sep='\t', index=False)
log(f'pending (f): {len(pend)} line x marker pairs', LOG)
# lines with >= 1,000 cells and no marker tested (both markers low)
st = S.set_index('cellosaurus_id')
tested = set(M.cellosaurus_id)
untested = st[(st.n_vehicle_cells >= 1000) & ~st.index.isin(tested)]
pd.DataFrame({'cell_name': untested.cell_name, 'cellosaurus_id': untested.index, 'state': untested.state,
              'n_vehicle_cells': untested.n_vehicle_cells,
              'reason': 'no marker above its cut-off (Y and XIST both low): no marker-positive population to test'}).to_csv(
    OUT + '/partA_local_mosaic_not_tested.tsv', sep='\t', index=False)
log(f'lines with >= 1,000 cells and no marker-high pair: {len(untested)}; lines below 1,000 cells: '
    + ', '.join(st[st.n_vehicle_cells < 1000].cell_name), LOG)
