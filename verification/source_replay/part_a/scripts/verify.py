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
print('All file checks and cohort, call, and STR counts passed.')
