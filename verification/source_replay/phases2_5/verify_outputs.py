"""Audit compound coverage and cross-check the tables behind all three phase figures."""
import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent

def rows(name):
    with (HERE / name).open() as f:
        return list(csv.DictReader(f))

def output(phase, index):
    with gzip.open(HERE / f'outputs/phase{phase}/{index:02}.csv.gz','rt') as f:
        return list(csv.DictReader(f))

names = json.loads((HERE/'prepared/drugs.json').read_text())
assert len(names) == len(set(names)) == 50
r2 = rows('comparison/phase2_all_tissue_tests.csv')
r3 = {r['drug']:r for r in rows('comparison/phase3_compounds.csv')}
r4 = {(r['drug'],r['lineage']):r for r in rows('comparison/phase4_pair_counts.csv')}
r5 = json.loads((HERE/'comparison/phase5.json').read_text())
r4hits = rows('comparison/phase4_changed_threshold_calls.csv')
checks = Counter()
for i, drug in enumerate(names):
    for phase in range(2,6):
        data = output(phase,i)
        assert all(r['drug'] == drug for r in data), (phase,drug)
        idcols = ['tissue'] if phase==2 else ['lineage','tissue']
        if phase==4:
            idcols.append('sex_direction')
        keys = {tuple(r[k] for k in idcols) for r in data}
        assert len(keys)==len(data), (phase,drug,'repeat')
        checks[f'phase{phase}_rows'] += len(data)
        if phase == 2:
            assert len(data) == 44 and len([r for r in r2 if r['drug']==drug]) == 44
        elif phase == 3:
            count = Counter(r['is_home_tissue'] for r in data if float(r['FDR']) < .05)
            record = r3[drug]
            assert (count['True'],count['False'])==(int(record['replay_A']),int(record['replay_B']))
            assert len({r['lineage'] for r in data}) == int(record['tested_lineages'])
        elif phase == 4:
            by = defaultdict(Counter)
            for r in data:
                col = by[r['lineage']]
                col['sets'] += 1
                if float(r['FDR']) < .05: col[r['is_home_tissue']] += 1
            for lin, col in by.items():
                record = r4[drug, lin]
                assert col['sets']==88
                assert (col['True'],col['False'])==(int(float(record['replay_A'])),int(float(record['replay_B'])))
            checks['phase4_pairs'] += len(by)
        else:
            for r in data: checks['phase5_'+r['category']] += 1
for category, n in r5['replay'].items():
    assert checks['phase5_'+category] == n
phase4 = json.loads((HERE/'comparison/phase4.json').read_text())
assert len(r4hits) == phase4['archive_only']+phase4['replay_only']
assert sum(r['original_significant']=='True' for r in r4hits) == phase4['archive_only']
assert checks['phase4_pairs']==len(r4)==phase4['pairs']
assert checks['phase4_rows']==phase4['all_tested_sets']
assert checks['phase3_rows']==json.loads((HERE/'comparison/phase3.json').read_text())['tissue_results']
assert checks['phase2_rows']==len(r2)
for file in ['phase3_within_lineage_counts','phase4_tested_reference_sets','phase5_testable_tissues']:
    for ext in ['png','pdf']:
        assert (HERE/'figures'/f'{file}.{ext}').stat().st_size > 1000
print(json.dumps(dict(checks),indent=2))
