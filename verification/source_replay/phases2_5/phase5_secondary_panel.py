"""Validate the separate panel in Phase 5 notebook cell 60 from its original selection criteria."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import os
import threading
import time

import fsspec
import h5py
import numpy as np
import pandas as pd
import requests
import replay as p

OUT=p.BASE/'phase5_secondary'
OUT.mkdir(exist_ok=True)
provider=json.loads((p.BASE/'input_verification/comparison.json').read_text())
URL=provider['source']

class ExtractTable(HTMLParser):
    def __init__(self):
        super().__init__()
        self.rows=[]; self.inside=False
    def handle_starttag(self,tag,attrs):
        if tag=='tr':
            self.rows.append([])
        if tag in ('td','th'):
            self.inside=True
            self.rows[-1].append('')
    def handle_endtag(self,tag):
        if tag in ('td','th'):
            self.inside=False
    def handle_data(self,data):
        if self.inside:
            self.rows[-1][-1]+=data

def get_panel():
    nb=json.loads(p.NB[5].read_text())
    h=next(o for o in nb['cells'][60]['outputs'] if 'text/html' in o.get('data',{}))
    rows=ExtractTable(); rows.feed(''.join(h['data']['text/html']))
    df=pd.DataFrame([r for r in rows.rows if r[0].isdecimal()],columns=['record','drug','tested_tissues','untested_tissues','category'])
    df.to_csv(OUT/'archived_notebook_panel.csv',index=False)
    return df

def prepare(archive):
    cell=pd.read_csv(p.META/'cellinfo_beta.txt',sep='\t')
    use=['sig_id','cmap_name','cell_iname','pert_type']
    # Keep just the fields and perturbation class needed for the selection.
    parts=[]
    for z in pd.read_csv(p.META/'siginfo_beta.txt',sep='\t',usecols=use,chunksize=220000):
        z=z[z.pert_type=='trt_cp'].merge(cell[['cell_iname','donor_sex','cell_lineage']],on='cell_iname',how='inner')
        parts.append(z[z.donor_sex.isin(['M','F']) & z.cell_lineage.ne('unknown')])
    sig=pd.concat(parts,ignore_index=True)
    ok=sig.groupby(['cmap_name','cell_lineage','donor_sex']).cell_iname.nunique().unstack(fill_value=0).reset_index()
    ok=ok[(ok.M>=5)&(ok.F>=5)].copy()
    ok['pair_total']=ok.M+ok.F
    ranked=ok.groupby('cmap_name').pair_total.sum().sort_values(ascending=False)
    top=ranked.head(50).index.tolist()
    stored=archive.drug.tolist()
    assert len(stored)==len(set(stored))==50
    assert all(ranked.loc[c]>=ranked.iloc[49] for c in stored)
    main=set(json.loads((p.BASE/'prepared/drugs.json').read_text()))
    extra=set(stored)-main
    sig=sig[sig.cmap_name.isin(extra)].copy()
    sig.to_parquet(OUT/'nine_compound_meta.parquet',index=False)
    need=sig.merge(ok[['cmap_name','cell_lineage']],on=['cmap_name','cell_lineage'])
    dest=OUT/'secondary_eligible_signature_ids.txt'
    dest.write_text(''.join(s+'\n' for s in sorted(set(need.sig_id))))
    detail={'eligible_compounds':len(ranked),'primary_overlap':len(set(stored)&main),'additional_compounds':sorted(extra),'eligible_signatures_requiring_source_lookup':need.sig_id.nunique(),'fresh_selection_same_set':set(top)==set(stored),'fresh_selection_same_order':top==stored,'archived_only':sorted(set(stored)-set(top)),'fresh_only':sorted(set(top)-set(stored)),'raw_score_at_rank_50':int(ranked.iloc[49])}
    (OUT/'selection.json').write_text(json.dumps(detail,indent=2)+'\n')
    print('Secondary panel',detail,flush=True)

def fetch():
    ids=set((OUT/'secondary_eligible_signature_ids.txt').read_text().splitlines())
    fs=fsspec.filesystem('https',block_size=4*1024**2)
    head=requests.head(URL,timeout=30)
    head.raise_for_status(); e=head.headers['ETag']
    assert e==provider['source_etag']
    with fs.open(URL,'rb') as f,h5py.File(f,'r') as h:
        rows=[v.decode() for v in h['0/META/ROW/id'][:]]
        allids=[v.decode() for v in h['0/META/COL/id'][:]]
        matrix=h['0/DATA/0/matrix'];off=matrix.id.get_offset()
        assert matrix.dtype==np.dtype('<f4') and matrix.chunks is None and list(matrix.shape)==provider['provider_shape']
    assert rows==(p.MATRIX/'row_ids.txt').read_text().splitlines()
    order=[(i,s) for i,s in enumerate(allids) if s in ids]
    assert len(order)==len(ids)
    names=[s for _,s in order]
    (OUT/'source_signature_order.txt').write_text(''.join(s+'\n' for s in names))
    n=len(order);stride=4*len(rows)
    # Memory-map a disposable raw cache to avoid retaining each HTTP body.
    cache_path=Path(os.environ.get('P59_ALTERNATE_CACHE',OUT/'eligible_source_signatures.npy'))
    data=np.lib.format.open_memmap(cache_path.with_suffix('.part'),mode='w+',dtype='<f4',shape=(n,len(rows)))
    groups=[]
    for i,(j,s) in enumerate(order):
        if not groups or j-groups[-1][0][0]>=128 or j-groups[-1][-1][0]>7:
            groups.append([])
        groups[-1].append((j,s,i))
    http=threading.local()
    def read(group):
        if not hasattr(http,'session'):
            http.session=requests.Session()
        lo=off+group[0][0]*stride;hi=off+(group[-1][0]+1)*stride-1
        for k in range(5):
            try:
                r=http.session.get(URL,headers={'Range':f'bytes={lo}-{hi}','If-Match':e},timeout=(20,60))
                assert r.status_code==206 and len(r.content)==hi-lo+1 and r.headers['Content-Range']==f"bytes {lo}-{hi}/{provider['source_file_bytes']}"
                break
            except (requests.RequestException,AssertionError):
                if k==4:raise
                time.sleep(2**k)
        raw=memoryview(r.content);rows_out=[]
        for j,s,i in group:
            offset=(j-group[0][0])*stride
            b=raw[offset:offset+stride];data[i]=np.frombuffer(b,dtype='<f4')
            rows_out.append([s,j,hashlib.sha256(b).hexdigest()])
        return rows_out
    with ThreadPoolExecutor(4) as workers:
        result=[r for block in workers.map(read,groups) for r in block]
    data.flush();del data; cache_path.with_suffix('.part').replace(cache_path)
    pd.DataFrame(result,columns=['sig_id','provider_position','sha256']).to_csv(OUT/'provider_row_checksums.csv',index=False)
    (OUT/'provider_receipt.json').write_text(json.dumps({'url':URL,'etag':e,'signatures':n,'gene_count':len(rows),'source_ranges':len(groups),'float32_values':n*len(rows)},indent=2)+'\n')
    print('Fetched and hashed',n,'signatures,',n*len(rows),'float32 entries',flush=True)

def analyze(archive):
    p.init()
    p.ids=(OUT/'source_signature_order.txt').read_text().splitlines();p.positions={s:i for i,s in enumerate(p.ids)}
    p.Z=np.load(os.environ.get('P59_ALTERNATE_CACHE',OUT/'eligible_source_signatures.npy'),mmap_mode='r')
    current=pd.read_parquet(OUT/'nine_compound_meta.parquet')
    library=set(p.union)
    main=json.loads((p.BASE/'prepared/drugs.json').read_text())
    (OUT/'extra_results').mkdir(exist_ok=True)
    rows=[]
    for old in archive.itertuples():
        drug=old.drug
        if drug in main:
            out=pd.read_csv(p.BASE/'sensitivity_gseapy_1.3.1/phase5'/f'{main.index(drug):02}.csv.gz')
        else:
            src=OUT/'extra_results'/(drug+'.csv.gz')
            if not src.exists():
                t=p.phase5(drug,current[current.cmap_name==drug])
                t.to_csv(src,index=False,compression='gzip')
            out=pd.read_csv(src)
        passed=set(out.loc[out.category!='D','tissue'])
        # Cell 60 includes each eligible lineage's home tissues regardless of GSEA set eligibility.
        testable=out.loc[(out.n_M>=5)&(out.n_F>=5)]
        for lin in testable.lineage.unique():
            passed.update(set(p.map5.get(lin,[]))&library)
        actual='A' if (out.category=='A').any() else 'B' if (out.category=='B').any() else 'C'
        orig_set=set(old.tested_tissues.split(', ')) if old.tested_tissues else set()
        orig_untested=set(old.untested_tissues.split(', ')) if old.untested_tissues else set()
        rows.append({'drug':drug,'archived_category':old.category,'raw_category':actual,'archived_n_tested':len(orig_set),'replay_n_tested':len(passed),'tested_set_agrees':orig_set==passed,'untested_set_agrees':orig_untested==library-passed})
    x=pd.DataFrame(rows);x.to_csv(OUT/'results_by_compound.csv',index=False)
    summary={'gseapy':p.gseapy.__version__,'original':x.archived_category.value_counts().to_dict(),'replay':x.raw_category.value_counts().to_dict(),'compounds_with_different_labels':x.loc[x.archived_category!=x.raw_category,'drug'].tolist(),'compounds_with_tested_tissue_differences':x.loc[~(x.tested_set_agrees&x.untested_set_agrees),'drug'].tolist()}
    (OUT/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n'); print(summary,flush=True)
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('steps',nargs='+',choices=['prepare','fetch','analyze']);args=ap.parse_args()
    a=get_panel()
    if 'prepare' in args.steps:prepare(a)
    if 'fetch' in args.steps:fetch()
    if 'analyze' in args.steps:analyze(a)
