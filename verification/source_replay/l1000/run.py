"""Rebuild the complete sampled LINCS Level 3 control-well Y-score check from the GCTX source."""
import argparse
import hashlib
import json
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import fsspec
import h5py
import numpy as np
import pandas as pd
import requests
from scipy.stats import rankdata, spearmanr
from fingerprints import cohort_digest, median_digest, digest_rows, MARKERS, Y_GENES

BASE = Path(__file__).resolve().parent
SEED = 20261114
URL = 'https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level3/level3_beta_ctl_n188708x12328.gctx'
YS = Y_GENES
RESOLVED = {'XY','XX-like','Y-lost','XIST-lost','discordant'}
LOCAL = threading.local()


def sha(path):
    with path.open('rb') as fh:
        return hashlib.file_digest(fh,'sha256').hexdigest()


def input_schema():
    response=requests.head(URL,timeout=90)
    response.raise_for_status()
    etag=response.headers['ETag']
    length=int(response.headers['Content-Length'])
    fs=fsspec.filesystem('https')
    with fs.open(URL,'rb',block_size=512*1024) as file, h5py.File(file,'r') as h:
        d=h['/0/DATA/0/matrix']
        assert d.chunks is None and d.compression is None
        samples=[x.decode() for x in h['/0/META/COL/id'][:]]
        genes=[x.decode() for x in h['/0/META/ROW/id'][:]]
        assert len(samples)==len(set(samples)) and len(genes)==len(set(genes))
        assert tuple(d.shape)==(len(samples),len(genes))
        schema=dict(url=URL,etag=etag,content_length=length,shape=list(d.shape),
                    data_offset=d.id.get_offset(),dtype=d.dtype.str,
                    gene_order_sha256=hashlib.sha256(('\n'.join(genes)+'\n').encode()).hexdigest(),
                    well_order_sha256=hashlib.sha256(('\n'.join(samples)+'\n').encode()).hexdigest())
        checks={i:hashlib.sha256(d[i,:].tobytes()).hexdigest() for i in [0,len(samples)-1]}
    return schema,samples,genes,checks


def choices(samples,states,data_root,reference):
    eligible=[]
    cols=['sample_id','pert_type','cmap_name','cell_iname','qc_pass']
    total=0
    for chunk in pd.read_csv(data_root/'instinfo_beta.txt',sep='\t',usecols=cols,chunksize=100000):
        total+=len(chunk)
        dmso=chunk[(chunk.pert_type=='ctl_vehicle')&(chunk.cmap_name=='DMSO')&(chunk.qc_pass==1)]
        eligible.append(dmso[['sample_id','cell_iname']])
    all_dmso=pd.concat(eligible)
    assert all_dmso.sample_id.is_unique
    row_ids={s:i for i,s in enumerate(samples)}
    wells=all_dmso[all_dmso.sample_id.isin(row_ids)]
    choices={name:g.sample_id.sort_values().to_numpy() for name,g in wells.groupby('cell_iname')}
    rng=np.random.default_rng(SEED)
    records=[]
    for r in states.itertuples():
        if r.state not in RESOLVED or r.cell_iname not in choices:
            continue
        choices_line=choices[r.cell_iname]
        selected=choices_line if len(choices_line)<=30 else np.sort(rng.choice(choices_line,30,replace=False))
        records.extend(dict(cell_iname=r.cell_iname,state=r.state,sample_id=w,
                            row_index=row_ids[w],n_available=len(choices_line)) for w in selected)
    cohort=pd.DataFrame(records).sort_values('row_index').reset_index(drop=True)
    original_match=cohort_digest(cohort)
    assert original_match==reference['cohort_sha256'], 'Selected sample cohort differs from the original'
    counts=dict(instance_metadata_rows=total,qc_pass_DMSO=len(all_dmso),qc_pass_DMSO_in_matrix=len(wells),
                cohort_rows=len(cohort),cohort_lines=cohort.cell_iname.nunique(),
                archived_cohort_matches_every_field=True,cohort_sha256=original_match)
    return cohort,counts


def read_known_rows(ids,schema,marker_idx):
    rowbytes=int(np.dtype(schema['dtype']).itemsize*schema['shape'][1])
    # Consecutive selected wells share one request without reading unused wells.
    groups=[]
    for index in sorted(set(ids)):
        if groups and index==groups[-1][-1]+1 and len(groups[-1])<25:
            groups[-1].append(index)
        else:
            groups.append([index])

    def one_group(group):
        lo,hi=group[0],group[-1]
        start=schema['data_offset']+lo*rowbytes
        stop=schema['data_offset']+(hi+1)*rowbytes-1
        expected_range=f'bytes {start}-{stop}/{schema["content_length"]}'
        if not hasattr(LOCAL,'http'):
            LOCAL.http=requests.Session()
        for attempt in range(6):
            try:
                response=LOCAL.http.get(URL,headers={'Range':f'bytes={start}-{stop}', 'If-Match':schema['etag']},
                                        timeout=(20,60))
                response.raise_for_status()
                assert response.status_code==206 and response.headers['Content-Range']==expected_range
                assert response.headers['ETag']==schema['etag'] and len(response.content)==stop-start+1
                arr=np.frombuffer(response.content,dtype=schema['dtype']).reshape(len(group),-1)
                return [(i,hashlib.sha256(row.tobytes()).hexdigest(),row[marker_idx].tolist())
                        for i,row in zip(group,arr)]
            except (requests.RequestException,AssertionError) as e:
                if attempt==5: raise
                time.sleep(1+2**attempt)
    rows={}
    with ThreadPoolExecutor(max_workers=6) as ex:
        for n,passed in enumerate(ex.map(one_group,groups),start=1):
            rows.update({idx:(digest,vector) for idx,digest,vector in passed})
            if n%200==0:
                print('Downloaded ranges',n,'of',len(groups),'well rows',len(rows),flush=True)
    return rows,rowbytes


def auc(scores,pos,neg):
    p=scores[pos].to_numpy(); n=scores[neg].to_numpy()
    ranks=rankdata(np.r_[p,n])
    estimate=(ranks[:len(p)].sum()-len(p)*(len(p)+1)/2)/(len(p)*len(n))
    direct=np.mean((p[:,None]>n).astype(float)+0.5*(p[:,None]==n))
    assert np.isclose(estimate,direct)
    return dict(positive_lines=len(p),negative_lines=len(n),auroc=float(estimate))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True,help='Directory with Broad geneinfo_beta.txt and instinfo_beta.txt')
    parser.add_argument('--states',type=Path,required=True,help='Task A line-state TSV, with full-precision CCLE Y scores')
    parser.add_argument('--reference',type=Path,default=BASE/'reference',help='Published Level 3 sampling and rounded marker medians')
    parser.add_argument('--output',type=Path,default=BASE/'results',help='Output directory')
    args=parser.parse_args()
    data_root,states_path,reference_root,out=args.data,args.states,args.reference,args.output
    out.mkdir(parents=True,exist_ok=True)
    ref_file=reference_root/'reference_fingerprints.json'
    reference=json.loads(ref_file.read_text())
    src,samples,genes,edge_checks=input_schema()
    print('Read GCTX metadata',src['shape'],src['etag'],flush=True)
    info=pd.read_csv(data_root/'geneinfo_beta.txt',sep='\t',dtype={'gene_id':str})
    assert set(info.gene_id)==set(genes)
    chosen=info[info.gene_symbol.isin(MARKERS)].set_index('gene_symbol')
    assert chosen.index.is_unique
    marker_columns=[genes.index(chosen.loc[s,'gene_id']) for s in MARKERS]
    assert len(set(marker_columns)) == len(MARKERS)
    marker_manifest={sym:dict(gene_id=genes[i],matrix_index=i,feature_space=chosen.loc[sym,'feature_space'])
                     for sym,i in zip(MARKERS,marker_columns)}
    states=pd.read_csv(states_path,sep='\t',dtype={'Yhi':str})
    cohort,counts=choices(samples,states,data_root,reference)
    print('Rebuilt selected cohort',counts,flush=True)
    fetched,rowbytes=read_known_rows(list(cohort.row_index)+list(edge_checks),src,marker_columns)
    assert all(digest==fetched[i][0] for i,digest in edge_checks.items())
    scored=cohort.copy()
    scored['gctx_row_sha256']=[fetched[i][0] for i in scored.row_index]
    for j,marker in enumerate(MARKERS):
        scored[marker]=[fetched[i][1][j] for i in scored.row_index]
    med=scored.groupby('cell_iname')[MARKERS].median()
    med['L1000_Y_score']=med[YS].mean(axis=1)
    line_digest=median_digest(med)
    assert line_digest==reference['rounded_marker_medians_sha256'], 'Archived four-decimal medians differ'
    joined=med.join(states.set_index('cell_iname')[['y_score_ccle','Yhi','state']],how='left',validate='1:1')
    both=joined.dropna(subset=['y_score_ccle','L1000_Y_score'])
    rho,p=spearmanr(both.y_score_ccle,both.L1000_Y_score)
    check_rho=np.corrcoef(rankdata(both.y_score_ccle),rankdata(both.L1000_Y_score))[0,1]
    assert np.isclose(rho,check_rho)
    scores=joined.L1000_Y_score
    metrics=dict(spearman_rho=float(rho),spearman_p=float(p),spearman_n=len(both),
                 original_AUC_XY_vs_Y_lost=auc(scores,joined.state=='XY',joined.state=='Y-lost'),
                 alternate_AUC_all_Y_high_vs_all_Y_low=auc(scores,joined.Yhi=='True',joined.Yhi=='False'))
    checks=dict(hdf5_first_last_checks={str(i):v for i,v in edge_checks.items()},
                selected_raw_bytes=len(cohort)*rowbytes,
                selected_well_rows_sha256=digest_rows(scored[['row_index','gctx_row_sha256']].itertuples(index=False,name=None)),
                rounded_marker_medians_sha256=line_digest,all_marker_medians_match_archive_to_four_decimals=True,
                filename_level='level3_beta_ctl_n188708x12328.gctx',sampling_seed=SEED,maximum_wells_per_line=30,
                aggregation='Per gene, compute the median across selected wells in each line; then average the eight Y-gene medians within line.',
                input_sha256={p.name:sha(p) for p in [states_path,data_root/'geneinfo_beta.txt',data_root/'instinfo_beta.txt',ref_file]})
    checks['metadata_sources']={name:'https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/'+name for name in ['geneinfo_beta.txt','instinfo_beta.txt']}
    receipt=dict(source=src,markers=marker_manifest,sampling=counts,verification=checks,results=metrics)
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(metrics,indent=2),flush=True)


if __name__=='__main__':
    main()
