"""Pseudobulk aggregation, foreign-cell filter loop, cut-offs, states and confidence."""
import itertools
import numpy as np
import pandas as pd
from common import (Y8, YP, YQ, scores_from_sums, kde_cutoff, rule_state, arm_state,
                    MIN_CELLS_CALL, plate_num)

SUM_COLS = ['n_cells', 'umi', 'XIST'] + Y8


def aggregate(cells):
    """Sums per (cell_line_id, plate, stratum). Cells already checked to be vehicle cells."""
    c = cells[['cell_line_id', 'plate', 'stratum', 'total_umi', 'XIST'] + Y8].copy()
    c['n_cells'] = 1
    c = c.rename(columns={'total_umi': 'umi'})
    for k in ['n_cells', 'umi', 'XIST'] + Y8:
        c[k] = c[k].astype(np.int64)
    agg = c.groupby(['cell_line_id', 'plate', 'stratum'], observed=True)[SUM_COLS].sum()
    return agg


def filtered(agg, yneg, xneg):
    """Drop strata of foreign cells: y8 >= 3 (strata 2, 3) in Y-negative lines and
    XIST >= 3 (strata 1, 3) in XIST-negative lines."""
    a = agg.reset_index()
    st = a['stratum'].to_numpy()
    yn = a['cell_line_id'].map(lambda l: yneg.get(l, False)).to_numpy(bool)
    xn = a['cell_line_id'].map(lambda l: xneg.get(l, False)).to_numpy(bool)
    drop = (yn & (st >= 2)) | (xn & ((st == 1) | (st == 3)))
    return a.loc[~drop]


def scores_at(a):
    lp = a.groupby(['cell_line_id', 'plate'], observed=True)[SUM_COLS].sum()
    ln = a.groupby(['cell_line_id'], observed=True)[SUM_COLS].sum()
    ren = {g: g + '_umi' for g in ['XIST'] + Y8}
    lps = pd.concat([lp.rename(columns=ren), scores_from_sums(lp)], axis=1)
    lns = pd.concat([ln.rename(columns=ren), scores_from_sums(ln)], axis=1)
    assert lps.columns.is_unique and lns.columns.is_unique
    return lps, lns


def call(ann, level, agg, max_iter=10):
    """level 'lineplate' (KDE over line x plate points, primary) or 'line' (KDE over line points)."""
    assert level in ('lineplate', 'line')
    lines = sorted(agg.index.get_level_values('cell_line_id').unique())
    status = None
    hist = []
    for it in range(1, max_iter + 1):
        yneg = {} if status is None else status[0]
        xneg = {} if status is None else status[1]
        a = filtered(agg, yneg, xneg)
        lps, lns = scores_at(a)
        pts = lps if level == 'lineplate' else lns
        cx = kde_cutoff(pts['XIST'])
        cy = kde_cutoff(pts['y_score'])
        if not (np.isfinite(cx['cutoff']) and np.isfinite(cy['cutoff'])):
            raise SystemExit(f'STOP: fewer than two modes at level {level}, iteration {it}')
        new_y = {l: bool(lns.loc[l, 'y_score'] < cy['cutoff']) for l in lines}
        new_x = {l: bool(lns.loc[l, 'XIST'] < cx['cutoff']) for l in lines}
        hist.append(dict(level=level, iteration=it,
                         filter_from='none (all cells)' if status is None else f'status of iteration {it-1}',
                         cut_XIST=cx['cutoff'], cut_y=cy['cutoff'],
                         n_Yneg=sum(new_y.values()), n_XISTneg=sum(new_x.values()),
                         cells_kept=int(a['n_cells'].sum()),
                         Yneg_lines=';'.join(sorted(l for l in lines if new_y[l])),
                         XISTneg_lines=';'.join(sorted(l for l in lines if new_x[l]))))
        if status is not None and new_y == status[0] and new_x == status[1]:
            converged = True
            break
        status = (new_y, new_x)
    else:
        converged = False
    hist = pd.DataFrame(hist)
    hist['converged'] = converged
    if not converged:
        raise SystemExit(f'STOP: filter loop did not converge at level {level}')
    # final: filter with the converged status
    yneg, xneg = status
    a = filtered(agg, yneg, xneg)
    lps, lns = scores_at(a)
    pts = lps if level == 'lineplate' else lns
    cuts = {s: kde_cutoff(pts[s]) for s in ['XIST', 'y_score', 'yp_score', 'yq_score']}
    cX, cY = cuts['XIST']['cutoff'], cuts['y_score']['cutoff']
    cP, cQ = cuts['yp_score']['cutoff'], cuts['yq_score']['cutoff']
    assert abs(cX - hist.iloc[-1]['cut_XIST']) < 1e-12 and abs(cY - hist.iloc[-1]['cut_y']) < 1e-12
    # plate states
    lps = lps.reset_index()
    lps = lps.merge(ann[['cell_line_id', 'annotation']], on='cell_line_id', how='left', validate='m:1')
    assert lps['annotation'].notna().all()
    lps['xist_high'] = lps['XIST'] >= cX
    lps['y_high'] = lps['y_score'] >= cY
    lps['yp_high'] = lps['yp_score'] >= cP
    lps['yq_high'] = lps['yq_score'] >= cQ
    lps['plate_state'] = [rule_state(s, x, y) for s, x, y in zip(lps.annotation, lps.xist_high, lps.y_high)]
    lps['plate_arm_state'] = [arm_state(p, q) for p, q in zip(lps.yp_high, lps.yq_high)]
    # line states
    lns = lns.reset_index().merge(ann, on='cell_line_id', how='left', validate='1:1')
    assert len(lns) == len(lines)
    out = []
    for _, r in lns.iterrows():
        l = r['cell_line_id']
        sex = r['annotation']
        x_hi = bool(r['XIST'] >= cX)
        y_hi = bool(r['y_score'] >= cY)
        yp_hi = bool(r['yp_score'] >= cP)
        yq_hi = bool(r['yq_score'] >= cQ)
        n_after = int(r['n_cells'])
        if n_after < MIN_CELLS_CALL:
            st = f'not called: fewer than {MIN_CELLS_CALL} cells'
        else:
            st = rule_state(sex, x_hi, y_hi)
        pl = lps[lps.cell_line_id == l]
        n_pl = len(pl)
        n_agree = int((pl['plate_state'] == st).sum())
        grid = sorted({rule_state(sex, r['XIST'] >= cX + dx, r['y_score'] >= cY + dy)
                       for dx, dy in itertools.product([-1, 0, 1], [-1, 0, 1])})
        grid_stable = (len(grid) == 1 and grid[0] == st)
        frac = n_agree / n_pl
        if st.startswith('not called'):
            conf = 'not called'
        elif frac == 1.0 and grid_stable:
            conf = 'high'
        elif frac >= 0.8:
            conf = 'medium'
        else:
            conf = 'low'
        vc = pl['plate_state'].value_counts()
        modal = ';'.join(sorted(vc[vc == vc.max()].index))
        out.append(dict(cellosaurus_id=l, state=st, arm_state=arm_state(yp_hi, yq_hi), confidence=conf,
                        plates_agree=round(frac, 4), grid_stable=grid_stable,
                        xist_high=x_hi, y_high=y_hi, yp_high=yp_hi, yq_high=yq_hi,
                        n_plates=n_pl, n_plates_agree=n_agree, modal_plate_state=modal,
                        plate_states=';'.join(f'{plate_num(p)}:{s}' for p, s in
                                              sorted(zip(pl.plate, pl.plate_state), key=lambda t: plate_num(t[0]))),
                        grid_states=';'.join(grid),
                        arm_plates_agree=round(float((pl['plate_arm_state'] == arm_state(yp_hi, yq_hi)).mean()), 4),
                        y_negative_filtered=yneg[l], xist_negative_filtered=xneg[l]))
    calls = lns.merge(pd.DataFrame(out), left_on='cell_line_id', right_on='cellosaurus_id', validate='1:1')
    cut_rows = []
    for s, c in cuts.items():
        d = dict(level=level, score=s)
        d.update(c)
        cut_rows.append(d)
    return dict(calls=calls, lineplate=lps, cuts=pd.DataFrame(cut_rows), hist=hist,
                yneg=yneg, xneg=xneg, filtered=a)
