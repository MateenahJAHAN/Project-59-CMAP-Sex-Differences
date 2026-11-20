"""Step 4 (after the state table was hashed): agreement with Task C, Task A and the pilot table;
label audit under both cut-off levels; PIE context mapping and its "Sex:" values."""
import os, sys, json, hashlib, re, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from common import OUT, INPUTS, log, norm_name

T0 = time.time()
LOG = open(OUT + '/logs/04_agreement_audit.log', 'a')
log('start 04_agreement_audit', LOG)
# the state table must be the hashed one
h = hashlib.sha256(open(OUT + '/partA_local_states.tsv', 'rb').read()).hexdigest()
h0 = open(OUT + '/partA_local_states.sha256').read().split()[0]
assert h == h0, 'state table changed after hashing'
S = pd.read_csv(OUT + '/partA_local_states.tsv', sep='\t')
L = pd.read_csv(OUT + '/partA_local_states_linelevel.tsv', sep='\t')
assert len(S) == 50 and S.cellosaurus_id.is_unique and set(S.cellosaurus_id) == set(L.cellosaurus_id)

# expressed pattern, independent of state names
def pattern_tahoe(state, annotation):
    if state == 'XY':
        return 'Y high'
    if state == 'XX-like':
        return 'XIST high, Y low'
    if state in ('Y-lost', 'XIST-lost'):
        return 'both low'
    if state == 'discordant':
        return 'Y high' if annotation == 'Female' else 'XIST high, Y low'
    if state == 'ambiguous':
        return 'both low'  # ambiguous arises only for both low with no label (rule table)
    return None

def pattern_ref(state):
    return {'XY': 'Y high', 'XX-like': 'XIST high, Y low', 'Y-lost': 'both low',
            'XIST-lost': 'both low'}.get(state, None)

# ---------------- Task C
C = pd.read_csv(INPUTS + '/C_cellline_sex_labels.tsv', sep='\t')
assert C.RRID.is_unique
Cs = C[['RRID', 'model_name', 'cellosaurus_sex', 'expressed_state', 'expressed_state_source', 'cmp_state', 'ccle_state']]
ag = S[['cell_name', 'cellosaurus_id', 'annotation', 'state', 'confidence']].merge(
    L[['cellosaurus_id', 'state', 'confidence']].rename(columns={'state': 'state_linelevel', 'confidence': 'confidence_linelevel'}),
    on='cellosaurus_id', validate='1:1')
ag = ag.merge(Cs, left_on='cellosaurus_id', right_on='RRID', how='left', validate='1:1').drop(columns='RRID')
ag = ag.rename(columns={'expressed_state': 'C_state', 'expressed_state_source': 'C_source',
                        'cmp_state': 'C_cmp_state', 'ccle_state': 'C_ccle_state', 'model_name': 'C_model_name'})
# ---------------- Task A
A = pd.read_csv(INPUTS + '/A_cellline_sex_states.tsv', sep='\t')
assert A.cell_iname.is_unique
a_ac = A.dropna(subset=['cellosaurus_ac'])
# Task A lists NCI-H1975 (CVCL_1511) twice; it is not a Tahoe line. Uniqueness is asserted per Tahoe id.
amap, akey = {}, {}
for l, nm in zip(S.cellosaurus_id, S.cell_name):
    hit = a_ac[a_ac.cellosaurus_ac == l]
    assert len(hit) <= 1, (l, hit.cell_iname.tolist())
    if len(hit) == 1:
        amap[l] = hit.index[0]; akey[l] = 'cellosaurus_ac'
        continue
    n = norm_name(nm)
    hit = A[(A.cell_iname.map(norm_name) == n) | (A.cellosaurus_name.map(norm_name) == n)]
    assert len(hit) <= 1, (nm, hit.cell_iname.tolist())
    if len(hit) == 1:
        amap[l] = hit.index[0]; akey[l] = 'normalised name'
ag['A_cell_iname'] = ag.cellosaurus_id.map(lambda l: A.loc[amap[l], 'cell_iname'] if l in amap else np.nan)
ag['A_join_key'] = ag.cellosaurus_id.map(akey)
ag['A_state'] = ag.cellosaurus_id.map(lambda l: A.loc[amap[l], 'state'] if l in amap else np.nan)
ag['A_confidence'] = ag.cellosaurus_id.map(lambda l: A.loc[amap[l], 'confidence'] if l in amap else np.nan)
ag['A_cellosaurus_ac'] = ag.cellosaurus_id.map(lambda l: A.loc[amap[l], 'cellosaurus_ac'] if l in amap else np.nan)
bad_ac = ag[(ag.A_join_key == 'normalised name') & ag.A_cellosaurus_ac.notna() & (ag.A_cellosaurus_ac != ag.cellosaurus_id)]
assert bad_ac.empty, bad_ac
# ---------------- pilot
P = pd.read_csv(INPUTS + '/tahoe_pilot_states.tsv', sep='\t')
assert P.cell_line_id.is_unique and len(P) == 50
ag = ag.merge(P[['cell_line_id', 'n_cells', 'xist', 'y_score', 'tahoe_state', 'confidence', 'plates', 'plates_agree',
                 'grid_stable', 'pie_sex', 'ann']].rename(columns={
    'cell_line_id': 'cellosaurus_id', 'n_cells': 'pilot_n_cells', 'xist': 'pilot_xist', 'y_score': 'pilot_y_score',
    'tahoe_state': 'pilot_state', 'confidence': 'pilot_confidence', 'plates': 'pilot_plates',
    'plates_agree': 'pilot_plates_agree', 'grid_stable': 'pilot_grid_stable', 'pie_sex': 'pilot_pie_sex',
    'ann': 'pilot_ann'}), on='cellosaurus_id', how='left', validate='1:1')
assert ag.pilot_state.notna().all()

ag['pattern'] = [pattern_tahoe(s, a) for s, a in zip(ag.state, ag.annotation)]
ag['C_pattern'] = ag.C_state.map(pattern_ref)
ag['A_pattern'] = ag.A_state.map(pattern_ref)
ag['agree_C'] = np.where(ag.C_state.isin(['not_measured']) | ag.C_state.isna(), np.nan, (ag.state == ag.C_state).astype(float))
ag['agree_C_linelevel'] = np.where(ag.C_state.isin(['not_measured']) | ag.C_state.isna(), np.nan,
                                   (ag.state_linelevel == ag.C_state).astype(float))
ag['agree_A'] = np.where(ag.A_state.isin(['no_data']) | ag.A_state.isna(), np.nan, (ag.state == ag.A_state).astype(float))
ag['agree_pilot_state'] = (ag.state == ag.pilot_state).astype(float)
ag['agree_pilot_confidence'] = (ag.confidence == ag.pilot_confidence).astype(float)
ag.to_csv(OUT + '/partA_local_agreement.tsv', sep='\t', index=False)

summ = []
def add(comp, mask, agree, note=''):
    n = int(mask.sum()); k = int(agree[mask].sum())
    summ.append(dict(comparison=comp, n_lines_with_both=n, n_agree=k, share=round(k / n, 4) if n else np.nan, note=note))

C_has = ag.C_state.notna() & (ag.C_state != 'not_measured')
add('Task C, state (primary cut-offs), ambiguous counted as a state [A1]', C_has, ag.state == ag.C_state,
    'A1 passes at >= 0.90; below 0.80 states are culture-specific')
noamb = C_has & (ag.C_state != 'ambiguous') & (ag.state != 'ambiguous')
add('Task C, state, ambiguous removed on either side', noamb, ag.state == ag.C_state)
add('Task C, state with line-level cut-offs', C_has, ag.state_linelevel == ag.C_state)
pat = C_has & ag.C_pattern.notna() & ag.pattern.notna()
add('Task C, expressed pattern (Y high / XIST high, Y low / both low)', pat, ag.pattern == ag.C_pattern,
    'secondary; independent of state names (Task C has no discordant state)')
A_has = ag.A_state.notna() & (ag.A_state != 'no_data')
add('Task A, state', A_has, ag.state == ag.A_state)
add('Task A, state, ambiguous removed on either side', A_has & (ag.A_state != 'ambiguous') & (ag.state != 'ambiguous'),
    ag.state == ag.A_state)
patA = A_has & ag.A_pattern.notna() & ag.pattern.notna()
add('Task A, expressed pattern', patA, ag.pattern == ag.A_pattern)
allp = pd.Series(True, index=ag.index)
add('pilot table, state', allp, ag.state == ag.pilot_state)
add('pilot table, confidence', allp, ag.confidence == ag.pilot_confidence)
add('pilot table, state and confidence', allp, (ag.state == ag.pilot_state) & (ag.confidence == ag.pilot_confidence))
add('pilot table, state (line-level cut-offs here)', allp, ag.state_linelevel == ag.pilot_state)
sm = pd.DataFrame(summ)
sm.to_csv(OUT + '/partA_local_agreement_summary.tsv', sep='\t', index=False)
for r in summ:
    log(f"{r['comparison']}: {r['n_agree']}/{r['n_lines_with_both']} = {r['share']}", LOG)

# ---------------- PIE mapping
pie_path = OUT + '/downloads/pie/descriptions.json'
raw = open(pie_path, 'rb').read()
ph = hashlib.sha256(raw).hexdigest()
assert ph.startswith('2b8875977f559faf') and len(raw) == 62725, (ph, len(raw))
D = json.loads(raw)
assert len(D) == 88
recs = []
for k, t in D.items():
    f = {}
    for part in t.split(' ;\n'):
        if ': ' in part:
            a, b = part.split(': ', 1)
            f[a.strip()] = b.strip().rstrip(' ;')
    syn = [s.strip() for s in f.get('Synonyms', '').split(',') if s.strip()]
    ids = set(re.findall(r'CVCL_[A-Z0-9]{4}', k + ' ' + t))
    recs.append(dict(context=k, recommended=f.get('Recommended Name', ''), sex=f.get('Sex', ''),
                     names={norm_name(x) for x in [k, f.get('Context Name', ''), f.get('Recommended Name', '')] + syn if x},
                     ids=ids, key_is_id=bool(re.fullmatch(r'CVCL_[A-Z0-9]{4}', k))))
pm = []
for l, nm, an in zip(S.cellosaurus_id, S.cell_name, S.annotation):
    by_key = [r for r in recs if r['context'] == l]
    by_id = [r for r in recs if l in r['ids']]
    by_name = [r for r in recs if norm_name(nm) in r['names']]
    if by_key:
        hit, how = by_key, 'context key = Cellosaurus id'
    elif by_id:
        hit, how = by_id, 'Cellosaurus id in context text'
    elif by_name:
        hit, how = by_name, 'normalised name'
    else:
        hit, how = [], 'not mapped'
    assert len(hit) <= 1, (nm, [r['context'] for r in hit])
    r = hit[0] if hit else None
    pm.append(dict(cell_name=nm, cellosaurus_id=l, annotation=an, pie_context=r['context'] if r else '',
                   pie_recommended_name=r['recommended'] if r else '', map_key=how,
                   pie_sex=r['sex'] if r else '', pie_sex_equals_annotation=(r['sex'] == an) if r else np.nan))
pm = pd.DataFrame(pm)
pm.to_csv(OUT + '/partA_local_pie_mapping.tsv', sep='\t', index=False)
n_map = int((pm.map_key != 'not mapped').sum())
log(f"PIE sha256 {ph}; contexts {len(D)}; Tahoe lines mapped {n_map}/50; by key {(pm.map_key=='context key = Cellosaurus id').sum()}, "
    f"by id in text {(pm.map_key=='Cellosaurus id in context text').sum()}, by name {(pm.map_key=='normalised name').sum()}", LOG)

# ---------------- label audit, both cut-off levels
usual = {'Male': 'XY', 'Female': 'XX-like'}
la = S[['cell_name', 'cellosaurus_id', 'annotation', 'annotation_source', 'XIST', 'y_score', 'state', 'confidence', 'xist_high']].merge(
    L[['cellosaurus_id', 'state', 'confidence']].rename(columns={'state': 'state_linelevel', 'confidence': 'confidence_linelevel'}),
    on='cellosaurus_id', validate='1:1').merge(pm[['cellosaurus_id', 'pie_context', 'pie_sex', 'map_key']],
                                               on='cellosaurus_id', validate='1:1')
la['labelled'] = la.annotation.isin(['Male', 'Female'])
la['usual_pattern'] = la.annotation.map(usual)
la['mismatch_primary'] = la.labelled & (la.state != la.usual_pattern)
la['mismatch_linelevel'] = la.labelled & (la.state_linelevel != la.usual_pattern)
la['male_XY_with_XIST_high'] = (la.annotation == 'Male') & (la.state == 'XY') & la.xist_high
la['pie_labelled'] = la.pie_sex.isin(['Male', 'Female'])
la['pie_usual_pattern'] = la.pie_sex.map(usual)
la['mismatch_pie_label_primary'] = la.pie_labelled & (la.state != la.pie_usual_pattern)
la.to_csv(OUT + '/partA_local_label_audit.tsv', sep='\t', index=False)
au = []
for lev, col in [('line x plate cut-offs (primary)', 'mismatch_primary'), ('line-level cut-offs', 'mismatch_linelevel')]:
    stcol = 'state' if col == 'mismatch_primary' else 'state_linelevel'
    sub = la[la.labelled]
    au.append(dict(level=lev, n_labelled=int(len(sub)), n_without_usual_pattern=int(sub[col].sum()),
                   share=round(sub[col].mean(), 4),
                   by_state=';'.join(f'{k}:{v}' for k, v in sub.loc[sub[col], stcol].value_counts().sort_index().items()),
                   lines=';'.join(sub.loc[sub[col], 'cell_name'])))
sub = la[la.pie_labelled]
au.append(dict(level='PIE Sex: label, line x plate cut-offs', n_labelled=int(len(sub)),
               n_without_usual_pattern=int(sub.mismatch_pie_label_primary.sum()),
               share=round(sub.mismatch_pie_label_primary.mean(), 4) if len(sub) else np.nan,
               by_state=';'.join(f'{k}:{v}' for k, v in sub.loc[sub.mismatch_pie_label_primary, 'state'].value_counts().sort_index().items()),
               lines=';'.join(sub.loc[sub.mismatch_pie_label_primary, 'cell_name'])))
au = pd.DataFrame(au)
au.to_csv(OUT + '/partA_local_label_audit_summary.tsv', sep='\t', index=False)
for r in au.itertuples():
    log(f'label audit {r.level}: {r.n_without_usual_pattern}/{r.n_labelled} ({r.by_state})', LOG)
log(f"male XY lines with XIST high: {int(la.male_XY_with_XIST_high.sum())} {';'.join(la.loc[la.male_XY_with_XIST_high,'cell_name'])}", LOG)
log(f'done 04_agreement_audit in {time.time()-T0:.1f} s', LOG)
