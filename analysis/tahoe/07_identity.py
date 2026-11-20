"""Step 6 (f): identity of deep marker-negative cells, six lines with full rows.
Pooled deep (>= 2,500 UMIs) marker-negative cells against depth-matched pools of the line's marker-positive cells,
beside a split-half of the positive cells and the best other line (other lines: their filtered vehicle pseudobulk
thinned to the same UMI total); nearest line centroid of single deep cells (multinomial log-likelihood)."""
import os, sys, time, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pyarrow.compute as pc
import scipy.sparse as sp
from common import OUT, DMSO, SEED, VEHICLE, Y8, MARKERS, log, load_cells, plate_num

DEEP = 2500
N_REP = 10
CPM_MIN = 10.0
T0 = time.time()
LOG = open(OUT + '/logs/07_identity.log', 'a')
log('start 07_identity', LOG)
S = pd.read_csv(OUT + '/partA_local_states.tsv', sep='\t').set_index('cellosaurus_id')
genes = pd.read_csv(OUT + '/gene_chrom.tsv', sep='\t')
gset = genes[genes.chrom_class.isin(['autosome', 'X']) & (genes.gene_symbol != 'XIST')]
assert 'XIST' not in set(gset.gene_symbol) and not set(Y8) & set(gset.gene_symbol)

# ---------------- filtered vehicle pseudobulk per line (all genes), for the gene filter and the centroids
pb = pq.read_table(DMSO + '/pseudobulk_counts.parquet', columns=['cell_line_id', 'stratum', 'token_id', 'count']).to_pandas()
pb['cell_line_id'] = pb['cell_line_id'].astype(str)
yneg = set(S.index[~S.y_high]); xneg = set(S.index[~S.xist_high])
drop = (pb.cell_line_id.isin(yneg) & (pb.stratum >= 2)) | (pb.cell_line_id.isin(xneg) & pb.stratum.isin([1, 3]))
pbf = pb.loc[~drop]
cent = pbf.groupby(['cell_line_id', 'token_id'], observed=True)['count'].sum()
tot_line = pbf.groupby('cell_line_id', observed=True)['count'].sum()
assert (tot_line.loc[S.index] == S['pseudobulk_umi']).all(), 'filtered pseudobulk total differs from the state table'
del pb, pbf
allsum = cent.groupby('token_id').sum()
cpm_all = 1e6 * allsum / allsum.sum()
G = np.array(sorted(set(gset.token_id) & set(cpm_all.index[cpm_all >= CPM_MIN])))
log(f'genes: autosomal+X without XIST {len(gset)}; with CPM >= {CPM_MIN} in the all-line vehicle pseudobulk: {len(G)}', LOG)
gidx = np.full(62713, -1, np.int64); gidx[G] = np.arange(len(G))
lines_all = sorted(S.index)
Cmat = np.zeros((len(lines_all), len(G)))
cdf = cent.reset_index()
cdf = cdf[cdf.token_id.isin(set(G))]
Cmat[cdf.cell_line_id.map({l: i for i, l in enumerate(lines_all)}).to_numpy(), gidx[cdf.token_id.to_numpy()]] = cdf['count'].to_numpy()
logP = np.log((Cmat + 1) / (Cmat.sum(axis=1, keepdims=True) + len(G)))
Pcent = Cmat / Cmat.sum(axis=1, keepdims=True)

# ---------------- full rows of the six lines
cells = load_cells(list(dict.fromkeys(['drug', 'cell_line_id', 'cell_name', 'plate', 'BARCODE_SUB_LIB_ID', 'total_umi', 'y8'] + MARKERS)))
assert cells.columns.is_unique
six = {'CVCL_0023': 'Y', 'CVCL_1724': 'XIST', 'CVCL_0459': 'Y', 'CVCL_0399': 'Y', 'CVCL_0504': 'XIST', 'CVCL_0320': 'XIST'}
for l, m in six.items():
    assert (S.loc[l, 'y_high'] if m == 'Y' else S.loc[l, 'xist_high']), (l, m)
c6 = cells[cells.cell_line_id.isin(set(six))].set_index('BARCODE_SUB_LIB_ID')
tok_mk = dict(zip(genes.gene_symbol, genes.token_id))
files = sorted(glob.glob(DMSO + '/six_lines/plate*.parquet'), key=lambda f: int(os.path.basename(f)[5:-8]))
assert len(files) == 14
store = {l: dict(bc=[], depth=[], mk=[], mats=[]) for l in six}
qc = []
seen = 0
for f in files:
    t = pq.read_table(f, columns=['genes', 'expressions', 'drug', 'BARCODE_SUB_LIB_ID', 'cell_line_id', 'plate'])
    drug = t.column('drug').to_numpy(zero_copy_only=False)
    if not (drug == VEHICLE).all():
        raise SystemExit(f'STOP: {f} has rows with drug != {VEHICLE}')
    bc = t.column('BARCODE_SUB_LIB_ID').to_numpy(zero_copy_only=False)
    lid = t.column('cell_line_id').to_numpy(zero_copy_only=False)
    gl = t.column('genes').combine_chunks(); el = t.column('expressions').combine_chunks()
    ln_g = pc.list_value_length(gl).to_numpy(zero_copy_only=False); ln_e = pc.list_value_length(el).to_numpy(zero_copy_only=False)
    assert (ln_g == ln_e).all()
    tok = pc.list_flatten(gl).to_numpy(); val = pc.list_flatten(el).to_numpy().astype(np.float64)
    row = np.repeat(np.arange(len(bc)), ln_g)
    del t, gl, el
    nong = tok < 3
    zero = (~nong) & (val == 0)
    keepm = (~nong) & (~zero)
    assert (val[keepm] > 0).all() and (val[keepm] == np.round(val[keepm])).all(), 'non-integer or negative counts'
    assert tok.max() <= 62712
    special = pd.Series(list(zip(tok[nong], val[nong]))).value_counts().to_dict()
    qc.append(dict(file=os.path.basename(f), rows=len(bc), nongene_entries=int(nong.sum()), nongene_kinds=str(special),
                   zero_entries=int(zero.sum()), cells_with_zero_entry=int(len(np.unique(row[zero])))))
    tok, val, row = tok[keepm], val[keepm], row[keepm]
    tot = np.bincount(row, weights=val, minlength=len(bc))
    ref = c6.loc[bc]
    assert (ref['cell_line_id'].to_numpy() == lid).all()
    assert (tot == ref['total_umi'].to_numpy()).all(), f'total UMI differs in {f}'
    for g in MARKERS:
        mm = tok == tok_mk[g]
        cnt = np.bincount(row[mm], weights=val[mm], minlength=len(bc))
        assert (cnt == ref[g].to_numpy()).all(), f'{g} differs in {f}'
    seen += len(bc)
    # keep deep cells after the primary foreign-cell filter (other marker only)
    for l, marker in six.items():
        r = S.loc[l]
        sel = (lid == l) & (tot >= DEEP)
        if not r.y_high and marker != 'Y':
            sel &= ref['y8'].to_numpy() < 3
        if not r.xist_high and marker != 'XIST':
            sel &= ref['XIST'].to_numpy() < 3
        rows_sel = np.where(sel)[0]
        if len(rows_sel) == 0:
            continue
        remap = np.full(len(bc), -1); remap[rows_sel] = np.arange(len(rows_sel))
        mm = (remap[row] >= 0) & (gidx[tok] >= 0)
        M = sp.csr_matrix((val[mm], (remap[row[mm]], gidx[tok[mm]])), shape=(len(rows_sel), len(G)))
        store[l]['bc'].append(bc[rows_sel]); store[l]['depth'].append(tot[rows_sel])
        store[l]['mk'].append((ref['y8'] if marker == 'Y' else ref['XIST']).to_numpy()[rows_sel])
        store[l]['mats'].append(M)
    log(f'{os.path.basename(f)}: {len(bc)} rows checked', LOG)
assert seen == len(c6) == 381695, (seen, len(c6))
pd.DataFrame(qc).to_csv(OUT + '/partA_local_six_lines_parse_qc.tsv', sep='\t', index=False)
log(f'six-line full rows: {seen} cells; total UMI and 22 marker counts equal cells_markers for every cell', LOG)


def corr_log(a, b):
    la = np.log2(1e6 * a / a.sum() + 1); lb = np.log2(1e6 * b / b.sum() + 1)
    return float(np.corrcoef(la, lb)[0, 1])


res, reps, nc_rows = [], [], []
for li, (l, marker) in enumerate(six.items()):
    st = store[l]
    depth = np.concatenate(st['depth']); mk = np.concatenate(st['mk']); M = sp.vstack(st['mats']).tocsr()
    name = S.loc[l, 'cell_name']
    neg = np.where(mk == 0)[0]; pos = np.where(mk >= 1)[0]
    n = int(min(len(neg), len(pos) // 2))
    rng = np.random.default_rng(SEED + 500 + li)
    others = [o for o in lines_all if o != l]
    oi = [lines_all.index(o) for o in others]
    rec_reps = []
    for rep in range(N_REP):
        nsel = rng.choice(neg, n, replace=False) if len(neg) > n else neg.copy()
        edges = np.quantile(depth[nsel], np.linspace(0, 1, 11))
        edges[0], edges[-1] = -np.inf, np.inf
        bn = np.clip(np.searchsorted(edges, depth[nsel], side='right') - 1, 0, 9)
        bp = np.clip(np.searchsorted(edges, depth[pos], side='right') - 1, 0, 9)
        A, B = [], []
        repl = False
        for k in range(10):
            need = int((bn == k).sum())
            if need == 0:
                continue
            cand = pos[bp == k]
            if len(cand) >= 2 * need:
                pick = rng.choice(cand, 2 * need, replace=False)
            else:
                repl = True
                pick = rng.choice(cand if len(cand) else pos, 2 * need, replace=True)
            A.extend(pick[:need]); B.extend(pick[need:])
        A, B = np.array(A), np.array(B)
        pn = np.asarray(M[nsel].sum(axis=0)).ravel(); pa = np.asarray(M[A].sum(axis=0)).ravel()
        pbb = np.asarray(M[B].sum(axis=0)).ravel()
        r_na, r_nb, r_ab = corr_log(pn, pa), corr_log(pn, pbb), corr_log(pa, pbb)
        Ntot = int(pn.sum())
        r_oth = []
        for o_i in oi:
            th = rng.multinomial(Ntot, Pcent[o_i]).astype(float)
            r_oth.append(corr_log(pn, th))
        r_oth = np.array(r_oth)
        best = int(np.argmax(r_oth))
        # positive pool against the best other line, as a reference
        th_pa = rng.multinomial(int(pa.sum()), Pcent[oi[best]]).astype(float)
        rec = dict(cell_name=name, cellosaurus_id=l, marker=marker, replicate=rep, n_pool=len(nsel),
                   neg_pool_umi=Ntot, posA_pool_umi=int(pa.sum()), posB_pool_umi=int(pbb.sum()),
                   neg_median_depth=float(np.median(depth[nsel])), posA_median_depth=float(np.median(depth[A])),
                   matched_with_replacement=repl, r_neg_posA=r_na, r_neg_posB=r_nb, r_neg_pos=(r_na + r_nb) / 2,
                   r_split_half=r_ab, best_other_line=S.loc[others[best], 'cell_name'], r_neg_best_other=float(r_oth[best]),
                   r_posA_best_other=corr_log(pa, th_pa), r_neg_median_other=float(np.median(r_oth)))
        rec_reps.append(rec)
    R = pd.DataFrame(rec_reps)
    reps.append(R)
    out = dict(cell_name=name, cellosaurus_id=l, marker=marker, n_deep=len(depth), n_deep_negative=len(neg),
               n_deep_positive=len(pos), n_pool=n, n_genes=len(G), n_replicates=N_REP,
               any_matching_with_replacement=bool(R.matched_with_replacement.any()))
    for c in ['r_neg_pos', 'r_split_half', 'r_neg_best_other', 'r_posA_best_other', 'r_neg_median_other']:
        out[c] = R[c].mean(); out[c + '_sd'] = R[c].std()
    out['best_other_line'] = R.best_other_line.mode().iloc[0]
    out['best_other_line_counts'] = ';'.join(f'{k}:{v}' for k, v in R.best_other_line.value_counts().items())
    out['f_gap_to_split_half'] = out['r_neg_pos'] - out['r_split_half']
    out['f_margin_over_best_other'] = out['r_neg_pos'] - out['r_neg_best_other']
    out['pass_f_split'] = out['f_gap_to_split_half'] >= -0.01
    out['pass_f_other'] = out['f_margin_over_best_other'] >= 0.10
    out['pass_f'] = out['pass_f_split'] and out['pass_f_other']
    # nearest centroid of single deep cells
    ll_neg = M[neg] @ logP.T
    own = lines_all.index(l)
    a_neg = np.asarray(ll_neg).argmax(axis=1)
    psel = rng.choice(pos, min(len(pos), 2000), replace=False)
    a_pos = np.asarray(M[psel] @ logP.T).argmax(axis=1)
    out['nc_neg_n'] = len(neg); out['nc_neg_own_share'] = float((a_neg == own).mean())
    out['nc_pos_n'] = len(psel); out['nc_pos_own_share'] = float((a_pos == own).mean())
    vc = pd.Series([S.loc[lines_all[i], 'cell_name'] for i in a_neg[a_neg != own]]).value_counts()
    out['nc_neg_other_assignments'] = ';'.join(f'{k}:{v}' for k, v in vc.head(5).items())
    for kind, arr in [('negative', a_neg), ('positive', a_pos)]:
        vc2 = pd.Series([S.loc[lines_all[i], 'cell_name'] for i in arr]).value_counts()
        for k, v in vc2.items():
            nc_rows.append(dict(cell_name=name, marker=marker, cells=kind, assigned_line=k, n=int(v), share=v / len(arr)))
    res.append(out)
    log(f"{name} {marker}: deep {len(depth)} (neg {len(neg)}, pos {len(pos)}), pool n {n}; r(neg,pos) {out['r_neg_pos']:.4f} "
        f"split {out['r_split_half']:.4f} best other {out['best_other_line']} {out['r_neg_best_other']:.4f}; "
        f"pass_f {out['pass_f']}; nearest centroid own: neg {out['nc_neg_own_share']:.3f}, pos {out['nc_pos_own_share']:.3f}", LOG)
pd.DataFrame(res).to_csv(OUT + '/partA_local_identity.tsv', sep='\t', index=False, float_format='%.6g')
pd.concat(reps).to_csv(OUT + '/partA_local_identity_replicates.tsv', sep='\t', index=False, float_format='%.6g')
pd.DataFrame(nc_rows).to_csv(OUT + '/partA_local_nearest_centroid.tsv', sep='\t', index=False, float_format='%.6g')
log(f'done 07_identity in {time.time()-T0:.1f} s', LOG)
