"""Check every cached signature and gene used by Phases 2–5 against the provider's HDF5 data bytes."""
import csv, hashlib, json, os, threading, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import fsspec,h5py,numpy as np,requests
HERE=Path(__file__).resolve().parent
CACHE=Path(os.environ['P59_MATRIX'])
URL='https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/level5_beta_trt_cp_n720216x12328.gctx'
def digest(x):return hashlib.sha256(x).hexdigest()
def main():
    info=requests.head(URL,timeout=45);info.raise_for_status()
    etag=info.headers['ETag'];length=int(info.headers['Content-Length'])
    fs=fsspec.filesystem('https',block_size=4*1024**2)
    with fs.open(URL,'rb') as f,h5py.File(f,'r') as h:
        signatures={x.decode():i for i,x in enumerate(h['0/META/COL/id'][:])}
        features=[x.decode() for x in h['0/META/ROW/id'][:]]
        d=h['0/DATA/0/matrix'];offset=d.id.get_offset();shape=d.shape
        assert d.chunks is None and d.dtype==np.dtype('<f4') and shape==(720216,12328)
    gids=(CACHE/'row_ids.txt').read_text().splitlines()
    ordered=(CACHE/'sig50.sigs.txt').read_text().splitlines()
    ids=sorted((signatures[s],i) for i,s in enumerate(ordered))
    cache=np.load(CACHE/'sig50.npy',mmap_mode='r')
    assert gids==features and len(set(ordered))==len(ids)==36523
    assert cache.shape==(len(ordered),len(gids)) and cache.dtype==np.dtype('<f4')
    stride=len(gids)*4
    groups=[]
    for sid,pos in ids:
        if not groups or sid-groups[-1][0][0]>=128 or sid-groups[-1][-1][0]>7:groups.append([])
        groups[-1].append((sid,pos))
    print('Provider ETag',etag,'signature blocks',len(groups),'estimated transport GB',sum((g[-1][0]-g[0][0]+1)*stride for g in groups)/1e9,flush=True)
    local=threading.local()
    def block(group):
        if not hasattr(local,'http'):local.http=requests.Session()
        lo=offset+stride*group[0][0];hi=offset+stride*(group[-1][0]+1)-1
        for k in range(6):
            try:
                resp=local.http.get(URL,headers={'Range':f'bytes={lo}-{hi}','If-Match':etag},timeout=(20,75))
                if resp.status_code!=206 or resp.headers.get('Content-Range')!=f'bytes {lo}-{hi}/{length}' or len(resp.content)!=hi-lo+1:
                    raise IOError((resp.status_code,resp.headers.get('Content-Range'),lo,hi))
                buf=memoryview(resp.content); rows=[]
                for sig,pos in group:
                    start=(sig-group[0][0])*stride
                    piece=buf[start:start+stride]
                    orig=memoryview(cache[pos]).cast('B')
                    a,b=digest(piece),digest(orig)
                    rows.append((ordered[pos],sig,a,b,a==b))
                return rows,len(buf)
            except (requests.RequestException,IOError):
                if k==5:raise
                time.sleep(2**k)
    import gzip
    rows_done=mismatches=transport=0
    temp=HERE/'signature_sha256.tsv.gz.part'
    with gzip.open(temp,'wt',newline='') as f,ThreadPoolExecutor(16) as pool:
        writer=csv.writer(f,delimiter='\t',lineterminator='\n')
        writer.writerow(['sig_id','provider_row','provider_float32_sha256','cache_float32_sha256','identical'])
        for i,(rows,nbytes) in enumerate(pool.map(block,groups),1):
            writer.writerows(rows);rows_done+=len(rows);transport+=nbytes;mismatches+=sum(not r[-1] for r in rows)
            if i%1000==0:print(i,'of',len(groups),'blocks;',rows_done,'signatures compared',flush=True)
    assert rows_done==len(ordered)
    temp.replace(HERE/'signature_sha256.tsv.gz')
    report=dict(source=URL,source_etag=etag,source_file_bytes=length,provider_shape=shape,matrix_offset=offset,
                all_gene_ids_and_order_match=gids==features,signatures_compared=rows_done,genes_per_signature=len(gids),
                float32_values_compared=rows_done*len(gids),mismatched_signatures=mismatches,bytes_received=transport,
                http_ranges=len(groups))
    (HERE/'comparison.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report,flush=True)
if __name__=='__main__':main()
