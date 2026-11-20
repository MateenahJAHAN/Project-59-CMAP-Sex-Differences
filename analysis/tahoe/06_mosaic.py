"""Step 6, tests (a) to (e), for every marker-high line with at least 1,000 vehicle cells after the filter.
Y for Y-high lines, XIST for XIST-high lines; COLO 205 (Y) is added as a labelled exploratory case."""
import os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from scipy.special import expit
from common import OUT, SEED, Y8, YP, YQ, X8, MIN_CELLS_MOSAIC, log, load_cells, plate_num
from mosaic_models import (fit_nb_mix, fit_pois_mix, deep_zero, mh_or, depth_deciles, nb_logpmf, DEPTHS)

COLS = ['drug', 'cell_line_id', 'cell_name', 'plate', 'total_umi', 'XIST', 'y8'] + Y8 + X8
BINS = [0, 500, 750, 1000, 1500, 2000, 3000, 4000, 6000, 8000, 12000, np.inf]


def run_one(task):
    line, name, marker, scope, idx, d, y, yp, yq, x8, plate, bg, bg_all = task
    rng = np.random.default_rng(SEED + 1000 * idx + (0 if marker == 'Y' else 1))
    out = dict(cell_name=name, cellosaurus_id=line, marker=marker, scope=scope, n_cells=len(d))
    pos = y >= 1
    r_pos = y[pos].sum() / d[pos].sum()
    out['r_pos'] = r_pos
    # (b) full data
    fb = fit_nb_mix(y, d, bg)
    for k, v in fb.items():
        out['b_' + k] = v
    # (a)
    for thr in (3000, 8000):
        a = deep_zero(y, d, thr, r_pos, (fb['r_high'], fb['k_high']), rng)
        for k, v in a.items():
            out[f'a{thr}_' + k] = v
    # (c)
    fc = fit_pois_mix(y, d, bg)
    out.update(c_pi=fc['pi'], c_r_high=fc['r_high'], c_iterations=fc['iterations'], c_converged=fc['converged'])
    # (d) odd against even plates
    pn = np.array([plate_num(p) for p in plate])
    for half, m in [('odd', pn % 2 == 1), ('even', pn % 2 == 0)]:
        fh = fit_nb_mix(y[m], d[m], bg[m])
        out[f'd_{half}_n_cells'] = int(m.sum())
        out[f'd_{half}_pi_2000'] = fh['pi_2000']
        out[f'd_{half}_converged'] = fh['converged']
        out[f'd_{half}_c_pi'] = fit_pois_mix(y[m], d[m], bg[m])['pi']
        out[f'd_{half}_a3000_zero_frac'] = float((y[m & (d >= 3000)] == 0).mean()) if (m & (d >= 3000)).any() else np.nan
    out['d_diff_pi_2000'] = out['d_odd_pi_2000'] - out['d_even_pi_2000']
    # (e) Y only
    if marker == 'Y':
        dec = depth_deciles(d)
        e = mh_or(yq == 0, yp == 0, dec)
        out.update(e_or_y=e['or_mh'], e_or_y_lo=e['or_lo'], e_or_y_hi=e['or_hi'], e_n_strata=e['n_strata'],
                   e_n_both_zero=e['n_both_zero'], e_mean_yp=yp.mean(), e_mean_yq=yq.mean(), e_mean_x8=x8.mean(),
                   e_zero_yp=float((yp == 0).mean()), e_zero_yq=float((yq == 0).mean()))
        pA, pB = yp.mean() / x8.mean(), yq.mean() / x8.mean()
        scaled = pA + pB > 1
        if scaled:
            s = pA + pB; pA, pB = pA / s, pB / s
        mult = rng.multinomial(x8.astype(np.int64), [pA, pB, max(0.0, 1 - pA - pB)])
        hA, hB = mult[:, 0], mult[:, 1]
        c1 = mh_or(hB == 0, hA == 0, dec)
        out.update(e_ctrl_pA=pA, e_ctrl_pB=pB, e_ctrl_scaled=scaled, e_or_ctrl=c1['or_mh'], e_or_ctrl_lo=c1['or_lo'],
                   e_or_ctrl_hi=c1['or_hi'], e_ctrl_zero_A=float((hA == 0).mean()), e_ctrl_zero_B=float((hB == 0).mean()))
        p50 = min(1.0, y.mean() / x8.mean())
        t = rng.binomial(x8.astype(np.int64), p50)
        h1 = rng.binomial(t, 0.5); h2 = t - h1
        c2 = mh_or(h2 == 0, h1 == 0, dec)
        out.update(e_or_ctrl5050=c2['or_mh'], e_or_ctrl5050_lo=c2['or_lo'], e_or_ctrl5050_hi=c2['or_hi'])
        out['e_ratio'] = out['e_or_y'] / out['e_or_ctrl']
        # background sensitivity: literal background from all Y-negative lines
        fb2 = fit_nb_mix(y, d, bg_all)
        out['b_pi_2000_bg_all_yneg'] = fb2['pi_2000']
        out['c_pi_bg_all_yneg'] = fit_pois_mix(y, d, bg_all)['pi']
    # depth curve
    curve = []
    pi_d = expit(fb['alpha'] + fb['beta'] * np.log(d / 2000.0))
    p0_low = np.exp(-d * bg)
    p0_nb = np.exp(nb_logpmf(0, d * fb['r_high'], fb['k_high']))
    model_zero = pi_d * p0_low + (1 - pi_d) * p0_nb
    bins = pd.cut(pd.Series(d), BINS, right=False, labels=False).to_numpy()
    for bi in range(len(BINS) - 1):
        m = bins == bi
        if not m.any():
            continue
        curve.append(dict(cell_name=name, cellosaurus_id=line, marker=marker, depth_lo=BINS[bi], depth_hi=BINS[bi + 1],
                          n_cells=int(m.sum()), median_depth=float(np.median(d[m])),
                          zero_frac=float((y[m] == 0).mean()), model_zero_frac=float(model_zero[m].mean()),
                          homogeneous_poisson_zero_frac=float(np.exp(-d[m] * r_pos).mean()),
                          homogeneous_nb_zero_frac=float(p0_nb[m].mean()),
                          b_pi_at_median=float(expit(fb['alpha'] + fb['beta'] * np.log(np.median(d[m]) / 2000.0))),
                          mean_marker_per_1000=float(1000 * y[m].sum() / d[m].sum())))
    return out, curve


def main():
    T0 = time.time()
    LOG = open(OUT + '/logs/06_mosaic.log', 'a')
    log('start 06_mosaic', LOG)
    cells = load_cells(COLS)
    S = pd.read_csv(OUT + '/partA_local_states.tsv', sep='\t').set_index('cellosaurus_id')
    bgt = pd.read_csv(OUT + '/partA_local_plate_background.tsv', sep='\t')
    bg_y = bgt[bgt.line_set.str.startswith('Y-negative lines with both')].set_index('unit')['per_umi']
    bg_y_all = bgt[bgt.line_set.str.startswith('Y-negative lines (y_score')].set_index('unit')['per_umi']
    bg_x = bgt[bgt.marker == 'XIST'].set_index('unit')['per_umi']
    assert len(bg_y) == len(bg_y_all) == len(bg_x) == 14
    lines_order = sorted(S.index)
    tasks = []
    for line in lines_order:
        r = S.loc[line]
        for marker in ('Y', 'XIST'):
            high = bool(r.y_high) if marker == 'Y' else bool(r.xist_high)
            scope = None
            if high and r.n_vehicle_cells >= MIN_CELLS_MOSAIC:
                scope = 'in scope'
            if marker == 'Y' and r.cell_name == 'COLO 205':
                scope = 'exploratory (Y below the cut-off; XY with line-level cut-offs)'
            if (marker == 'Y' and r.cell_name == 'A549') or (marker == 'XIST' and r.cell_name == 'SW48'):
                assert scope == 'in scope', (r.cell_name, marker)
            if scope is None:
                continue
            c = cells[cells.cell_line_id == line]
            # foreign-cell filter on the OTHER marker only (the tested marker is never filtered)
            keep = np.ones(len(c), bool)
            if not r.y_high and marker != 'Y':
                keep &= (c['y8'] < 3).to_numpy()
            if not r.xist_high and marker != 'XIST':
                keep &= (c['XIST'] < 3).to_numpy()
            c = c[keep]
            if scope == 'in scope':
                assert len(c) == r.n_vehicle_cells, (r.cell_name, len(c), r.n_vehicle_cells)
            d = c['total_umi'].to_numpy(float)
            y = (c['y8'] if marker == 'Y' else c['XIST']).to_numpy(float)
            yp = c[YP].sum(axis=1).to_numpy(float)
            yq = c[YQ].sum(axis=1).to_numpy(float)
            x8 = c[X8].sum(axis=1).to_numpy(float)
            plate = c['plate'].to_numpy()
            bg = c['plate'].map(bg_y if marker == 'Y' else bg_x).to_numpy(float)
            bga = c['plate'].map(bg_y_all if marker == 'Y' else bg_x).to_numpy(float)
            assert np.isfinite(bg).all() and (bg > 0).all()
            tasks.append((line, r.cell_name, marker, scope, lines_order.index(line), d, y, yp, yq, x8, plate, bg, bga))
    log(f'{len(tasks)} line x marker tests: ' + ', '.join(f'{t[1]} {t[2]}' for t in tasks), LOG)
    results, curves = [], []
    with ProcessPoolExecutor(max_workers=5) as ex:
        for out, curve in ex.map(run_one, tasks):
            results.append(out); curves.extend(curve)
            log(f"{out['cell_name']} {out['marker']}: b pi2000 {out['b_pi_2000']:.3f} c {out['c_pi']:.3f} "
                f"a3000 {out['a3000_zero_frac']:.3f} a8000 {out['a8000_zero_frac']:.3f} d diff {out['d_diff_pi_2000']:.3f}"
                + (f" e OR {out['e_or_y']:.2f} ctrl {out['e_or_ctrl']:.2f}" if out['marker'] == 'Y' else ''), LOG)
    R = pd.DataFrame(results)
    sh = R[['a3000_zero_frac', 'a8000_zero_frac', 'b_pi_2000', 'c_pi']]
    R['share_range_lo'] = sh.min(axis=1); R['share_range_hi'] = sh.max(axis=1)
    sh2 = R[['a3000_zero_frac', 'a8000_zero_frac', 'b_pi_1000', 'b_pi_2000', 'b_pi_4000', 'b_pi_8000', 'c_pi']]
    R['share_range_all_lo'] = sh2.min(axis=1); R['share_range_all_hi'] = sh2.max(axis=1)
    R['pass_b'] = R['b_pi_2000'].between(0.15, 0.85)
    R['pass_d'] = R['d_diff_pi_2000'].abs() <= 0.10
    R['pass_e'] = np.where(R['marker'] == 'Y', R.get('e_or_y', np.nan) > 2 * R.get('e_or_ctrl', np.nan), np.nan)
    R.to_csv(OUT + '/partA_local_mosaic_abcde.tsv', sep='\t', index=False, float_format='%.6g')
    pd.DataFrame(curves).to_csv(OUT + '/partA_local_depth_curve.tsv', sep='\t', index=False, float_format='%.6g')
    log(f'done 06_mosaic in {time.time()-T0:.1f} s', LOG)


if __name__ == '__main__':
    main()
