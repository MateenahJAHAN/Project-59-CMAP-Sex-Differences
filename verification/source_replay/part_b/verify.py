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
