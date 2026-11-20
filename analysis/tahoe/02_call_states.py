"""Step 3 of the frozen rules: input checks, line x plate pseudobulks, cut-offs, foreign-cell filter,
states, arm states and confidence. Writes partA_local_states.tsv and its sha256 BEFORE any comparison.
"""
import os, sys, json, hashlib, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from common import (DMSO, META, OUT, INPUTS, VEHICLE, Y8, MARKERS, MIN_CELLS_CALL, log, load_cells)
from calling import aggregate, call

T0 = time.time()
os.makedirs(OUT + '/logs', exist_ok=True)
LOG = open(OUT + '/logs/02_call_states.log', 'a')
log('start 02_call_states', LOG)
checks = []


def check(name, ok, detail=''):
    checks.append(dict(check=name, ok=bool(ok), detail=str(detail)))
    log(f'check {name}: {"OK" if ok else "FAIL"} {detail}', LOG)
    if not ok:
        pd.DataFrame(checks).to_csv(OUT + '/partA_local_input_checks.tsv', sep='\t', index=False)
        raise SystemExit(f'STOP: check failed: {name} {detail}')


# ---------------- cells and vehicle-only checks
cells = load_cells()  # stops unless every row has drug DMSO_TF
check('cells_markers rows', len(cells) == 2205786, len(cells))
check('cells_markers drug == DMSO_TF for every row', (cells['drug'] == VEHICLE).all(),
      cells['drug'].value_counts().to_dict())
sm = pd.read_parquet(META + '/sample_metadata.parquet', columns=['sample', 'plate', 'drug'])
dmso_samples = sm[sm['drug'] == VEHICLE]
check('sample_metadata DMSO_TF wells', len(dmso_samples) == 28, len(dmso_samples))
smap = dict(zip(dmso_samples['sample'], dmso_samples['plate']))
check('every cell sample is a DMSO_TF well', cells['sample'].isin(set(smap)).all(),
      cells.loc[~cells['sample'].isin(set(smap)), 'sample'].unique()[:5])
check('cell plate equals the well plate in sample_metadata',
      (cells['sample'].map(smap) == cells['plate']).all())
check('key unique', cells['key'].is_unique)
check('BARCODE_SUB_LIB_ID unique', cells['BARCODE_SUB_LIB_ID'].is_unique)
check('total_umi == obs_tscp_count', (cells['total_umi'] == cells['obs_tscp_count']).all())
check('n_genes == obs_gene_count', (cells['n_genes'] == cells['obs_gene_count']).all())
check('y8 == sum of the 8 Y genes', (cells[Y8].sum(axis=1) == cells['y8']).all())
st = 2 * (cells['y8'] >= 3).astype(int) + (cells['XIST'] >= 3).astype(int)
check('stratum == 2*(y8>=3)+(XIST>=3)', (st == cells['stratum']).all())
check('marker counts non-negative', (cells[MARKERS] >= 0).all().all())
check('total_umi positive', (cells['total_umi'] > 0).all())
n_lines = cells['cell_line_id'].nunique()
n_plates = cells['plate'].nunique()
check('50 lines, 14 plates, 28 wells', (n_lines, n_plates, cells['sample'].nunique()) == (50, 14, 28),
      (n_lines, n_plates, cells['sample'].nunique()))
lp_n = cells.groupby(['cell_line_id', 'plate'], observed=True).size()
check('every line on every plate (700 pairs)', len(lp_n) == 700, len(lp_n))
check('cell_name constant within cell_line_id', cells.groupby('cell_line_id', observed=True)['cell_name'].nunique().max() == 1)

# ---------------- pseudobulk consistency (also proves the pseudobulk holds DMSO wells only)
grp = pd.read_csv(DMSO + '/pseudobulk_groups.tsv', sep='\t')
check('pseudobulk_groups samples are DMSO_TF wells', grp['sample'].isin(set(smap)).all())
cg = cells.groupby(['sample', 'cell_line_id', 'stratum'], observed=True).agg(
    n_cells=('total_umi', 'size'), umi=('total_umi', 'sum')).reset_index()
m = cg.merge(grp, on=['sample', 'cell_line_id', 'stratum'], how='outer', suffixes=('_cells', '_grp'),
             indicator=True, validate='1:1')
check('pseudobulk_groups match cell groups 1:1', (m['_merge'] == 'both').all() and len(m) == 2910, len(m))
check('group n_cells equal', (m['n_cells_cells'] == m['n_cells_grp']).all())
check('group umi equal', (m['umi_cells'] == m['umi_grp']).all())

genes = pd.read_csv(DMSO + '/genes.tsv', sep='\t')
sym2tok = {s: t for s, t in zip(genes['gene_symbol'], genes['token_id']) if s in set(MARKERS)}
check('each marker has one token', len(sym2tok) == len(MARKERS) and
      all((genes['gene_symbol'] == s).sum() == 1 for s in MARKERS))
t1 = time.time()
pbt = pq.read_table(DMSO + '/pseudobulk_counts.parquet', columns=['sample', 'cell_line_id', 'stratum', 'token_id', 'count'])
pb = pbt.to_pandas()
del pbt
pb['count'] = pb['count'].astype(np.int64)
check('pseudobulk_counts samples are DMSO_TF wells', set(pb['sample'].astype(str).unique()) <= set(smap))
check('pseudobulk tokens are gene tokens (>= 3)', int(pb['token_id'].min()) >= 3, int(pb['token_id'].min()))
check('pseudobulk counts positive', int(pb['count'].min()) > 0, int(pb['count'].min()))
tot = pb.groupby(['sample', 'cell_line_id', 'stratum'], observed=True)['count'].sum().reset_index()
tot['sample'] = tot['sample'].astype(str)
tot['cell_line_id'] = tot['cell_line_id'].astype(str)
mm = cg.merge(tot, on=['sample', 'cell_line_id', 'stratum'], how='outer', indicator=True, validate='1:1')
check('pseudobulk count sum == cells total_umi per group', (mm['_merge'] == 'both').all() and (mm['umi'] == mm['count']).all())
tok2sym = {t: s for s, t in sym2tok.items()}
pm = pb[pb['token_id'].isin(set(tok2sym))].copy()
pm['gene'] = pm['token_id'].map(tok2sym)
pmw = pm.pivot_table(index=['sample', 'cell_line_id', 'stratum'], columns='gene', values='count',
                     aggfunc='sum', fill_value=0, observed=True).reset_index()
pmw['sample'] = pmw['sample'].astype(str)
pmw['cell_line_id'] = pmw['cell_line_id'].astype(str)
cm = cells.groupby(['sample', 'cell_line_id', 'stratum'], observed=True)[MARKERS].sum().reset_index()
j = cm.merge(pmw, on=['sample', 'cell_line_id', 'stratum'], how='left', suffixes=('', '_pb'), validate='1:1')
bad = [g for g in MARKERS if not (j[g] == j.get(g + '_pb', pd.Series(0, index=j.index)).fillna(0)).all()]
check('pseudobulk marker-token sums == cell marker sums (22 markers)', len(bad) == 0, bad)
del pb, pm, pmw
log(f'pseudobulk checks took {time.time()-t1:.1f} s', LOG)

# ---------------- annotation
lines = pd.read_csv(DMSO + '/lines.tsv', sep='\t')
check('lines.tsv has the 50 cell lines of the cells', set(lines['cell_line_id']) == set(cells['cell_line_id']))
C = pd.read_csv(INPUTS + '/C_cellline_sex_labels.tsv', sep='\t', usecols=['RRID', 'model_name', 'cellosaurus_sex'])
check('Task C RRID unique', C['RRID'].is_unique)
ann = lines[['cell_line_id', 'cell_name', 'organ']].merge(C, left_on='cell_line_id', right_on='RRID', how='left',
                                                           validate='1:1')
ann['annotation'] = ann['cellosaurus_sex']
ann['annotation_source'] = np.where(ann['RRID'].notna(), 'Task C cellosaurus_sex (join on RRID)', '')
miss = ann['annotation'].isna()
check('lines without a Task C row', miss.sum() == 1 and set(ann.loc[miss, 'cell_name']) == {'hTERT-HPNE'},
      list(ann.loc[miss, 'cell_name']))
txt = open(OUT + '/downloads/cellosaurus/CVCL_C466.txt').read().splitlines()
ac = [l for l in txt if l.startswith('AC   ')]
sx = [l[5:].strip() for l in txt if l.startswith('SX   ')]
check('Cellosaurus entry is CVCL_C466 with one SX line', ac and ac[0][5:].strip() == 'CVCL_C466' and len(sx) == 1, (ac, sx))
prov = pd.read_csv(OUT + '/provenance_downloads.tsv', sep='\t')
acc = prov.loc[prov['item'] == 'cellosaurus_CVCL_C466', 'accessed_utc'].iloc[-1]
ann.loc[miss, 'annotation'] = sx[0]
ann.loc[miss, 'annotation_source'] = f'Cellosaurus entry CVCL_C466, SX line (api.cellosaurus.org, accessed {acc})'
check('annotation values allowed', ann['annotation'].isin(['Male', 'Female', 'Unknown']).all(),
      ann['annotation'].value_counts().to_dict())
ann = ann[['cell_line_id', 'cell_name', 'organ', 'annotation', 'annotation_source']]

# ---------------- calling (primary: line x plate cut-offs)
agg = aggregate(cells)
check('aggregate keeps every cell', int(agg['n_cells'].sum()) == len(cells))
res = call(ann, 'lineplate', agg)
calls = res['calls']
n_before = cells.groupby('cell_line_id', observed=True).size()
calls['n_cells_before_filter'] = calls['cell_line_id'].map(n_before)
removed = calls['n_cells_before_filter'] - calls['n_cells']
cy = (cells['y8'] >= 3).groupby(cells['cell_line_id']).sum()
cx = (cells['XIST'] >= 3).groupby(cells['cell_line_id']).sum()
calls['n_removed_Y'] = [int(cy[l]) if yn else 0 for l, yn in zip(calls.cell_line_id, calls.y_negative_filtered)]
calls['n_removed_XIST'] = [int(cx[l]) if xn else 0 for l, xn in zip(calls.cell_line_id, calls.xist_negative_filtered)]
both = ((cells['y8'] >= 3) & (cells['XIST'] >= 3)).groupby(cells['cell_line_id']).sum()
calls['n_removed_both'] = [int(both[l]) if (yn and xn) else 0 for l, yn, xn in
                           zip(calls.cell_line_id, calls.y_negative_filtered, calls.xist_negative_filtered)]
check('removed cells add up', (removed == calls['n_removed_Y'] + calls['n_removed_XIST'] - calls['n_removed_both']).all())

table = pd.DataFrame({
    'cell_name': calls['cell_name'], 'cellosaurus_id': calls['cell_line_id'], 'organ': calls['organ'],
    'annotation': calls['annotation'], 'annotation_source': calls['annotation_source'],
    'n_vehicle_cells': calls['n_cells'].astype(int),
    'XIST': calls['XIST'].round(4), 'y_score': calls['y_score'].round(4),
    'yp_score': calls['yp_score'].round(4), 'yq_score': calls['yq_score'].round(4),
    'state': calls['state'], 'arm_state': calls['arm_state'], 'confidence': calls['confidence'],
    'plates_agree': calls['plates_agree'], 'grid_stable': calls['grid_stable'],
    'n_vehicle_cells_before_filter': calls['n_cells_before_filter'].astype(int),
    'n_removed_Y_ge3': calls['n_removed_Y'], 'n_removed_XIST_ge3': calls['n_removed_XIST'],
    'pseudobulk_umi': calls['umi'].astype(np.int64),
    'xist_high': calls['xist_high'], 'y_high': calls['y_high'], 'yp_high': calls['yp_high'], 'yq_high': calls['yq_high'],
    'n_plates': calls['n_plates'], 'n_plates_agree': calls['n_plates_agree'],
    'modal_plate_state': calls['modal_plate_state'], 'grid_states': calls['grid_states'],
    'arm_plates_agree': calls['arm_plates_agree'], 'plate_states': calls['plate_states'],
})
table = table.sort_values('cell_name', key=lambda s: s.str.lower()).reset_index(drop=True)
check('one row per line', len(table) == 50 and table['cellosaurus_id'].is_unique)
check('states allowed', (table['state'].isin(['XY', 'XX-like', 'Y-lost', 'XIST-lost', 'discordant', 'ambiguous'])
                         | table['state'].str.startswith('not called')).all(), table['state'].value_counts().to_dict())

# ---------------- write the state table and its hash FIRST
path = OUT + '/partA_local_states.tsv'
if os.path.exists(path):
    raise SystemExit('STOP: partA_local_states.tsv exists; it must not be overwritten (write a new version instead)')
table.to_csv(path, sep='\t', index=False)
written = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
h = hashlib.sha256(open(path, 'rb').read()).hexdigest()
with open(OUT + '/partA_local_states.sha256', 'w') as fh:
    fh.write(f'{h}  partA_local_states.tsv\n')
    fh.write(f'written_utc {written}\n')
    fh.write(f'hashed_utc {time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}\n')
os.chmod(path, 0o444)
log(f'WROTE {path} sha256 {h} at {written}', LOG)

# ---------------- auxiliary outputs (no comparison with any other table)
res['cuts'].to_csv(OUT + '/partA_local_cutoffs.tsv', sep='\t', index=False)
res['hist'].to_csv(OUT + '/partA_local_filter_iterations.tsv', sep='\t', index=False)
lpo = res['lineplate'].merge(ann[['cell_line_id', 'cell_name']], on='cell_line_id', validate='m:1')
lpo.to_csv(OUT + '/partA_local_lineplate_scores.tsv', sep='\t', index=False, float_format='%.5f')
pd.DataFrame(checks).to_csv(OUT + '/partA_local_input_checks.tsv', sep='\t', index=False)
ann.to_csv(OUT + '/partA_local_annotation.tsv', sep='\t', index=False)
consort = [
    dict(step='vehicle cells read (drug DMSO_TF, pass_filter full)', cells=len(cells), lines=n_lines),
    dict(step='foreign cells removed (y8>=3 in Y-negative lines; XIST>=3 in XIST-negative lines)',
         cells=int(removed.sum()), lines=int((removed > 0).sum())),
    dict(step='vehicle cells kept for calling', cells=int(calls['n_cells'].sum()), lines=n_lines),
    dict(step=f'lines with >= {MIN_CELLS_CALL} cells after the filter (called)', cells=int(calls.loc[calls.n_cells >= MIN_CELLS_CALL, 'n_cells'].sum()),
         lines=int((calls.n_cells >= MIN_CELLS_CALL).sum())),
    dict(step='lines with >= 1,000 cells after the filter (mosaic tests)', cells=int(calls.loc[calls.n_cells >= 1000, 'n_cells'].sum()),
         lines=int((calls.n_cells >= 1000).sum())),
]
pd.DataFrame(consort).to_csv(OUT + '/partA_local_consort.tsv', sep='\t', index=False)
log(f'done 02_call_states in {time.time()-T0:.1f} s', LOG)
