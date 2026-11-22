"""Verify the public files, control calculations, design checks and archive comparisons."""
from pathlib import Path
import csv,hashlib,json,math
root=Path(__file__).resolve().parent
for line in (root/'SHA256SUMS').read_text().splitlines():
    checksum,name=line.split('  ',1)
    assert hashlib.sha256((root/name).read_bytes()).hexdigest()==checksum,name
with (root/'rebuild_checks.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
arrays=[r for r in rows if r['stage']=='outcome_tensors']
assert len(arrays)==45 and all(r['rebuild_sha256']==r['archived_sha256'] for r in arrays)
with (root/'model_design_comparison.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
assert len(rows)==30 and all(r['mask_lineage_equal']==r['matrix_equal']=='True' for r in rows)
checks=json.loads((root/'ploidy_control_comparison.json').read_text())
assert len(checks)==3 and all(r['median_split_equal'] and not r['differences'] for r in checks)
for r in checks:
    result=r['observation'];b=result['B'];p=result['panel_p']
    assert b==10000 and math.isclose(p,(1+round(p*(b+1)-1))/(b+1),rel_tol=1e-14)
with (root/'X2_result.tsv').open() as f: row=next(csv.DictReader(f,delimiter='\t'))
assert (row['n_Y1'],row['n_Y0'],row['testable'])==('2','18','False')
print('Verified all package hashes, 45 tensor objects, 30 design comparisons, 3 full ploidy controls and the high-dose eligibility count.')

primary = root/'direct_verification.json'
if primary.is_file():
    primary_results=json.loads(primary.read_text())
    assert {k:v['rows'] for k,v in primary_results.items()} == {'H1':2, 'H3':2, 'per_compound':232}
    assert all(not entry['diffs'] for entry in primary_results.values())
    with (root/'direct_H1.tsv').open() as f: h1=list(csv.DictReader(f, delimiter='\t'))
    assert {(r['term'], int(r['B']), int(r['n_lines']), int(r['n_tested'])) for r in h1} == {('factY',10000,41,114),('factX',10000,41,114)}
    print('Primary panel statistics, the pooled estimates and all 232 compound/control rows match.')

conditions = root/'idarubicin_conditions_check.json'
if conditions.is_file():
    values=json.loads(conditions.read_text())
    with (root/'X2all_idarubicin_conditions.tsv').open() as f:
        tested=list(csv.DictReader(f, delimiter='\t'))
    assert len(tested)==values['conditions']==3
    assert hashlib.sha256((root/'X2all_idarubicin_conditions.tsv').read_bytes()).hexdigest()==values['replayed_sha256']==values['archived_sha256']
    assert [row['testable'] for row in tested]==['True','True','False']
    print('All three exploratory idarubicin conditions match the archived table.')

distributions=root/'figure4_nulls_check.json'
if distributions.is_file():
    import gzip
    with gzip.open(root/'figure4_nulls.tsv.gz','rt') as fp:
        data=list(csv.DictReader(fp,delimiter='\t'))
    reference=json.loads(distributions.read_text())
    assert list(map(int,(row['draw'] for row in data)))==list(range(10000))
    for term,r in reference.items():
        assert r['draws']==10000 and r['unequal_draws']==0 and r['p_matches_archive'] and r['observed_matches_archive']
        observed = r['observed']
        # The frozen panel test counts numerical ties at the engine's relative tolerance.
        exceedances=sum(float(row[term])>=observed*(1-1e-5) for row in data)
        assert exceedances==r['at_or_above'] and r['panel_p']==(exceedances+1)/10001
    intervals=json.loads((root/'figure4_intervals_check.json').read_text())
    assert intervals['rows']==12 and not intervals['different_cells']
    print('Figure 4 reproduces 30,000 draws and 12 estimates with both confidence intervals.')
