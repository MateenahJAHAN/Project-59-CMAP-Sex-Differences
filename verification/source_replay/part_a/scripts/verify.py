"""Verify the entire aggregate package using the Python standard library."""
from pathlib import Path
import csv, hashlib, json
base=Path(__file__).resolve().parents[1]
def read(path):
    with (base/path).open(newline='') as file:return list(csv.DictReader(file,delimiter='\t'))
def receipt(name):return json.loads((base/'receipts'/name).read_text())
for row in read('SHA256SUMS.tsv'):
    path=base/row['file']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==row['sha256'],path
sources=receipt('source_verification.json')
assert len(sources['sciplex_GEO_files'])==16 and sources['emerald_exact_file_matches']==116
src=receipt('sciplex_source_coverage.json')
assert sum(x['matrix_rows'] for x in src.values())==1171291549
assert sum(x['n_vehicle'] for x in src.values())==23570
obs=receipt('tahoe_observations.json')
assert (obs['rows'],obs['full_vehicle_rows'],obs['missing_or_different'])==(100648790,2205786,0)
em=read('receipts/emerald_vehicle_coverage.tsv')
assert len(em)==116 and sum(int(r['vehicle_rows']) for r in em)==130498
lines=read('tables/F1_emerald_states.tsv')
assert len(lines)==52 and sum(int(r['n_vehicle_cells']) for r in lines)==100253
assert len(read('tables/F1_emerald_units.tsv'))==207
for fn, pos, tot in [('str_fresh_EmeraldBay.tsv',47,52),('str_fresh_sci-Plex.tsv',3,3)]:
    rr=read('tables/'+fn)
    assert sum(r['match']=='True' for r in rr)==pos and sum(r['comparable']=='True' for r in rr)==tot
assert len(read('receipts/cellosaurus_record_checks.tsv'))==1102
z=read('tables/sciplex_vs_fresh_bulk.tsv')
assert 'bulk_y' not in z[0] and 'bulk_xist' not in z[0]
assert len(z)==6 and [r['experiment']+'/'+r['line'] for r in z if r['agree']=='False']==['sciPlex3/A549']
assert all(v['new_sha256']==v['old_sha256'] for v in receipt('compare_amelogenin.json').values())
for name,expected in [('CCLE',(773,734,714)),('CMP',(1082,1031,1007))]:
    rr=[r for r in read('tables/str_fresh_bulk.tsv') if r['source']==name]
    assert (len(rr),sum(r['comparable']=='True' for r in rr),sum(r['agree']=='True' for r in rr))==expected
    assert all(r['call_changed']=='False' for r in rr)
primary=receipt('primary_caller_comparison.json')
assert len(primary)==33 and all(v['new_sha256']==v['old_sha256'] for v in primary.values())
agg=receipt('tahoe_aggregate_comparison.json')
assert agg['partial_tables']==7640 and not agg['different_partial_tables'] and not agg['per_shard_differences']
assert len(agg['aggregate_outputs'])==21 and all(z['byte_identical'] for z in agg['aggregate_outputs'].values())
raw=receipt('tahoe_aggregate_checks.json')['checks']
assert raw['cells']==2205786 and raw['line_plate_pairs_equal_obs']==raw['line_plate_pairs']==700
labels=read('tables/partA_local_label_audit_summary.tsv')
assert any(r['n_labelled']=='49' and r['n_without_usual_pattern']=='27' for r in labels)
bulk=read('tables/patterns_tahoe_bulk.tsv')
assert len(bulk)==49 and sum(r['agree']=='True' for r in bulk)==47
cross=read('tables/patterns_cross_atlas.tsv')
assert len(cross)==47 and sum(r['match']=='True' for r in cross)==42
assert sum(r['eligible']=='True' for r in cross)==41 and sum(r['eligible']==r['match']=='True' for r in cross)==40
d=receipt('str_figure_inputs.json')
assert sum(v['fields_compared'] for v in d['tables'].values())==14351
assert not any(v['changed_fields'] for v in d['tables'].values())
for r in read('tables/figure3_agreement_by_source.tsv'):
    dataset=r['source'];rows=read('tables/figure3_supp_table_amelogenin_by_line.tsv')
    counts={k:sum(x[dataset+'_vs_STR']==k for x in rows) for k in ['agree','type A','type B']}
    assert int(r['n_compared'])==sum(counts.values()) and int(r['n_agree'])==counts['agree']
print('All file checks and cohort, call, and STR counts passed.')
