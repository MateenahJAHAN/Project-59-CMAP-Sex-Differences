"""Verify the saved evidence offline, including AUC and rank correlation from all scored lines."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as handle:
        return list(csv.DictReader(handle,delimiter='\t'))


def ranks(vals):
    ordered=sorted(range(len(vals)),key=vals.__getitem__)
    result=[0.0]*len(vals)
    first=0
    while first < len(vals):
        end=first+1
        while end<len(vals) and vals[ordered[end]]==vals[ordered[first]]:
            end+=1
        for i in ordered[first:end]:
            result[i]=(first+1+end)/2
        first=end
    return result


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'results')
    parser.add_argument('--reference',type=Path,default=HERE/'reference')
    parser.add_argument('--states',type=Path,default=HERE/'reference/line_state_inputs.tsv')
    args=parser.parse_args()
    result=json.loads((args.output/'receipt.json').read_text())
    for name,digest in result['verification']['output_sha256'].items():
        assert sha(args.output/name)==digest,name
    assert sha(args.states)==result['verification']['input_sha256'][args.states.name]
    for fname in ('A_l1000_samples_used.tsv.gz','A_l1000_line_medians.tsv'):
        assert sha(args.reference/fname)==result['verification']['input_sha256'][fname]
    cohort=read(args.output/'sampled_well_hashes.tsv.gz')
    original=read(args.reference/'A_l1000_samples_used.tsv.gz')
    assert [{k:v for k,v in row.items() if k!='gctx_row_sha256'} for row in cohort]==original
    assert len({r['sample_id'] for r in cohort})==len(cohort)==2800
    assert len({r['row_index'] for r in cohort})==len(cohort)
    assert all(len(r['gctx_row_sha256'])==64 and int(r['gctx_row_sha256'],16)>=0 for r in cohort)
    per_line=Counter(r['cell_iname'] for r in cohort)
    assert max(per_line.values())<=30 and len(per_line)==127
    inputs={r['cell_iname']:r for r in read(args.states)}
    scores=read(args.output/'inferred_scores.tsv')
    medians={r['cell_iname']:r for r in read(args.reference/'A_l1000_line_medians.tsv')}
    assert set(per_line)==set(medians)=={r['cell_iname'] for r in scores}
    assert all(abs(float(r['L1000_Y_score'])-float(medians[r['cell_iname']]['L1000_Y_score']))<.000051 for r in scores)
    stats=result['results']
    both=[r for r in scores if inputs[r['cell_iname']]['y_score_ccle']]
    x=[float(inputs[r['cell_iname']]['y_score_ccle']) for r in both]
    y=[float(r['L1000_Y_score']) for r in both]
    assert len(y)==stats['spearman_n']==124
    assert math.isclose(statistics.correlation(ranks(x),ranks(y)),stats['spearman_rho'])
    for field,positive,negative,metric in [
        ('state','XY','Y-lost','original_AUC_XY_vs_Y_lost'),
        ('Yhi','True','False','alternate_AUC_all_Y_high_vs_all_Y_low')]:
        p,n=[[float(r['L1000_Y_score']) for r in scores if inputs[r['cell_iname']][field]==value] for value in (positive,negative)]
        score=sum((x>y)+.5*(x==y) for x in p for y in n)/(len(p)*len(n))
        m=stats[metric]
        assert (len(p),len(n))==(m['positive_lines'],m['negative_lines'])
        assert math.isclose(score,m['auroc'])
    print('Verified all 2,800 sampled controls, 127 inferred line scores, both AUCs and Spearman correlation.')


if __name__=='__main__':
    main()
