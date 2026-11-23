"""Compare every public Phases 2–5 screen row with the corresponding fresh run."""
import argparse
import csv
import gzip
import json
import math
from pathlib import Path

from run import outside_checkout

RELEASE = Path(__file__).resolve().parents[1] / 'verification/source_replay/phases2_5'
KEYS = {
    2: ('drug','tissue'),
    3: ('drug','lineage','tissue'),
    4: ('drug','lineage','tissue','sex_direction'),
    5: ('drug','lineage','tissue'),
}
FLOATS = {'NES','FDR','p_value','fold','q_all','q_original'}


def read(path, keys):
    with gzip.open(path,'rt') as file:
        rows = csv.DictReader(file)
        result = {}
        for row in rows:
            key = tuple(row[k] for k in keys)
            if key in result:
                raise ValueError(f'Duplicate comparison in {path}: {key}')
            result[key] = row
        return rows.fieldnames, result


def compare(run, phase, pairs, report):
    keys = KEYS[phase]
    output = dict(compounds=len(pairs), matched_rows=0, missing_rows=0, additional_rows=0,
                  changed_rows=0, greatest_numeric_difference=0.0)
    for name, new_i, old_i in pairs:
        path = f'outputs/phase{phase}'
        a, latest = read(run/path/f'{new_i:02}.csv.gz', keys)
        b, reference = read(RELEASE/path/f'{old_i:02}.csv.gz', keys)
        if a!=b:
            raise ValueError(f'Phase {phase}, {name}: column names differ from the release')
        matched = latest.keys() & reference.keys()
        output['matched_rows'] += len(matched)
        output['missing_rows'] += len(reference.keys()-latest.keys())
        output['additional_rows'] += len(latest.keys()-reference.keys())
        for key in sorted(reference.keys()-latest.keys()):
            report.writerow([phase, name, '/'.join(key[1:]), 'missing row', '', ''])
        for key in sorted(latest.keys()-reference.keys()):
            report.writerow([phase, name, '/'.join(key[1:]), 'additional row', '', ''])
        for key in sorted(matched):
            differences = []
            for field in a:
                x, y = latest[key][field], reference[key][field]
                if x==y:
                    continue
                if field in FLOATS and x and y:
                    delta = abs(float(x)-float(y))
                    output['greatest_numeric_difference'] = max(output['greatest_numeric_difference'],delta)
                    if math.isclose(float(x),float(y),abs_tol=1e-12,rel_tol=1e-12):
                        continue
                differences.append(field)
                report.writerow([phase,name,'/'.join(key[1:]),field,y,x])
            if differences:
                output['changed_rows'] += 1
    return output


def check_gene_tests(run, names):
    folder = run/'phase2_gene_tests'
    with (run/'comparison/phase2_compounds.csv').open() as f:
        expected = {r['drug']:int(r['n_strong_de']) for r in csv.DictReader(f)}
    total = 0
    max_error = 0.0
    for i, drug in enumerate(names):
        with gzip.open(folder/f'{i:02}.csv.gz','rt') as f:
            genes = list(csv.DictReader(f))
        if len(genes)!=12328 or len({r['gene_id'] for r in genes})!=len(genes) or {r['drug'] for r in genes}!={drug}:
            raise ValueError(f'Phase 2 gene table has a different cohort for {drug}')
        ranked = sorted(genes, key=lambda row:float(row['p_value']))
        q = 1.0
        for position in range(len(ranked),0,-1):
            row = ranked[position-1]
            p = float(row['p_value'])
            if not 0.0<=p<=1.0:
                raise ValueError(f'Invalid gene-test p value for {drug}')
            q = min(q,p*len(ranked)/position)
            max_error = max(max_error,abs(float(row['p_adj'])-q))
        selected = {r['ensembl_id'] for r in genes if r['ensembl_id'] and
                    float(r['p_adj'])<.05 and abs(float(r['difference']))>.5}
        if len(selected)!=expected[drug]:
            raise ValueError(f'Gene table and tissue tests disagree about the strong gene list for {drug}')
        total += len(genes)
    if max_error > 1e-12:
        raise ValueError(f'Adjusted gene-test p values differ by up to {max_error}')
    return {'gene_tests':total,'max_absolute_BH_difference':max_error}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, type=outside_checkout)
    parser.add_argument('--phases', nargs='+', type=int, choices=KEYS, default=sorted(KEYS))
    args=parser.parse_args()
    new=json.loads((args.out/'prepared/drugs.json').read_text())
    archived=json.loads((RELEASE/'prepared/drugs.json').read_text())
    if len(new)!=len(archived) or len(new)!=50 or set(new)!=set(archived) or len(set(new))!=len(new):
        raise ValueError('Expected the same 50 distinct compounds in both panels')
    pairs = [(drug,i,archived.index(drug)) for i,drug in enumerate(new)]
    comparison = args.out/'comparison'
    comparison.mkdir(exist_ok=True)
    results = {}
    for n in args.phases:
        with (comparison/f'phase{n}_release_differences.csv').open('w', newline='') as fp:
            report = csv.writer(fp)
            report.writerow(['phase','drug','comparison','field','released','rerun'])
            result = compare(args.out,n,pairs,report)
        if n==2 and (args.out/'phase2_gene_tests').exists():
            result.update(check_gene_tests(args.out,new))
        (comparison/f'phase{n}_release_counts.json').write_text(json.dumps(result,indent=2)+'\n')
        results[str(n)] = result
    print(json.dumps(results, indent=2))


if __name__=='__main__':
    main()
