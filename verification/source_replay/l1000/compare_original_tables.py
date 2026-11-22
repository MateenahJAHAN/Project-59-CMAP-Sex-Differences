"""Check the complete archived exploratory table against a replay of the original Task A pipeline."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path


def read(path):
    with (gzip.open(path,'rt') if path.suffix=='.gz' else path.open()) as f:
        return list(csv.DictReader(f,delimiter='\t'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive',type=Path,required=True)
    parser.add_argument('--replay',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    gene_file='A_l1000_all_gene_auc_exploratory.tsv.gz'
    a,b=[p/gene_file for p in (args.archive,args.replay)]
    u,v=read(a),read(b)
    x,y=[{r['gene_id']:r for r in rows} for rows in (u,v)]
    assert x.keys()==y.keys() and len(u)==len(x) and len(v)==len(y)
    changes={col:sum(x[i][col]!=y[i][col] for i in x) for col in u[0]}
    assert all(n==0 for n in changes.values())
    ranks_stable=all(i['rank_by_discrimination']==j['rank_by_discrimination'] and i['discrimination']==j['discrimination'] for i,j in zip(u,v))
    assert ranks_stable
    data=dict(source_table=gene_file,archive_sha256=sha(a),rerun_sha256=sha(b),unique_genes=len(x),
              changed_cells_by_field_after_join_on_gene_id=changes,
              row_positions_with_different_gene_ids=sum(i['gene_id']!=j['gene_id'] for i,j in zip(u,v)),
              all_reordered_positions_have_same_rank_and_discrimination=ranks_stable,
              distinct_discrimination_values=len({i['discrimination'] for i in u}))
    exploratory='A_l1000_exploratory.tsv'
    variants=[p/exploratory for p in (args.replay,args.archive)]
    seen=[next(r for r in read(p) if r['analysis']=='state_contrast_auc' and r['gene']=='L1000_Y_score'
               and r['contrast']=='XY (positive) vs Y-lost') for p in variants]
    assert seen[0]==seen[1]
    m=seen[0]
    data['original_XY_vs_Y_lost_bootstrap']=dict(contrast=m['contrast'],auc=float(m['auc']),
                ci_low=float(m['ci_low']),ci_high=float(m['ci_high']),positive=int(float(m['n_pos'])),
                negative=int(float(m['n_neg'])),source=exploratory,source_sha256=sha(variants[0]),
                archived_sha256=sha(variants[1]))
    args.output.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':
    main()
