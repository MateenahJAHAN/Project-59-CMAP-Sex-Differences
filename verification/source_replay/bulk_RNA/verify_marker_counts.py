"""Recount the GDSC2/PRISM marker comparisons from the classified bulk-RNA data."""
import argparse
import csv
from collections import Counter
from pathlib import Path

def rows(path):
    with path.open() as file:
        yield from csv.DictReader(file, delimiter='\t')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('task_c_results', type=Path, help='Directory with the bulk state table and expression cutoffs')
    args = parser.parse_args()
    source = args.task_c_results
    cuts = {(r['expr_source'],r['variable']):float(r['cutoff_log2']) for r in rows(source/'C_expression_cutoffs.tsv')}
    counts = {label:Counter() for label in ['Male','Female']}
    for line in rows(source/'C_cellline_sex_labels.tsv'):
        label=line['annotated_sex']
        if label not in counts or not line['lineage'] or line['sex_specific_lineage'] != 'False':
            continue
        c=counts[label]
        c['all_eligible']+=1
        assay=line['expressed_state_source']
        if assay not in {'CCLE','CMP'}:
            c['no_rna']+=1
            continue
        y=float(line[assay.lower()+'_yscore'])>=cuts[assay,'yscore']
        x=float(line[assay.lower()+'_xist'])>=cuts[assay,'xist']
        c['rna_measured']+=1
        c['y_low']+=int(not y)
        c['xist_low']+=int(not x)
        c[{(True,True):'both_high',(True,False):'y_only',(False,True):'xist_only',(False,False):'neither_high'}[y,x]]+=1
    for expected in rows(Path(__file__).with_name('marker_counts.tsv')):
        label=expected.pop('annotated_sex')
        for key,value in expected.items():
            if counts[label][key] != int(value):
                raise ValueError(f'{label} {key}: observed {counts[label][key]}, expected {value}')
        print(label,dict(counts[label]))
if __name__=='__main__':
    main()
