"""Re-execute phases 2–5 using the original notebook functions and public Level 5 scores.
Set P59_META, P59_GTEX and P59_MATRIX to the metadata, GTEx and expression directories.
Outputs are separated by compound so completed jobs can be resumed.
"""
import argparse, ast, csv, hashlib, itertools, json, logging, os
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import scipy, gseapy
from scipy.stats import mannwhitneyu, fisher_exact
from statsmodels.stats.multitest import multipletests

BASE=Path(__file__).resolve().parent
REPO=Path(os.environ.get('P59_PUBLIC', BASE.parents[2]))
PREPARED=Path(os.environ.get('P59_PREPARED', BASE/'prepared'))
META=Path(os.environ['P59_META']); GTEX=Path(os.environ['P59_GTEX']); MATRIX=Path(os.environ['P59_MATRIX'])
NB={2:next(REPO.glob('Phase2*/*.ipynb')),3:next(REPO.glob('Phase3*/*.ipynb')),4:next(REPO.glob('Phase4*/*new.ipynb')),5:next(REPO.glob('Phase5*/*.ipynb'))}
warnings.filterwarnings('ignore',category=FutureWarning)
logging.disable(logging.WARNING)

def notebook_code(phase, cell):
    return ''.join(json.loads(NB[phase].read_text())['cells'][cell]['source'])

def prepare():
    cell=pd.read_csv(META/'cellinfo_beta.txt',sep='\t')
    fields=['sig_id','cell_iname','cmap_name','pert_type']
    sig=pd.read_csv(META/'siginfo_beta.txt',sep='\t',usecols=fields,low_memory=False)
    sig=sig.merge(cell[['cell_iname','donor_sex','cell_lineage']],on='cell_iname',how='left')
    sig=sig[(sig.pert_type=='trt_cp')&(sig.donor_sex.isin(['M','F']))]
    summary=sig.pivot_table(index='cmap_name',columns='donor_sex',values='cell_iname',aggfunc='nunique',fill_value=0)
    summary['total']=summary.M+summary.F
    names=summary[(summary.M>=5)&(summary.F>=5)].sort_values('total',ascending=False).head(50).index.tolist()
    hist=json.loads((REPO/'Phase5_PerSignature_Categories/reference/top50_names.json').read_text())
    assert set(names)==set(hist),(set(names)-set(hist),set(hist)-set(names))
    sig=sig[sig.cmap_name.isin(names)]
    assert len(sig)==len(set(sig.sig_id))==36523
    assert set(sig.sig_id)==set((MATRIX/'sig50.sigs.txt').read_text().splitlines())
    PREPARED.mkdir(exist_ok=True, parents=True)
    sig.to_parquet(PREPARED/'signatures.parquet',index=False)
    (PREPARED/'drugs.json').write_text(json.dumps(names,indent=2)+'\n')
    print('Selected the same 50 compounds and 36,523 signatures as the archive. Compound files follow prepared/drugs.json.', flush=True)

def init():
    global Z, ids, positions, gene_ids, known, names, gene_info, gene_map, union, split, simplemap, map5, universe
    Z=np.load(MATRIX/'sig50.npy',mmap_mode='r'); ids=(MATRIX/'sig50.sigs.txt').read_text().splitlines()
    positions={s:i for i,s in enumerate(ids)};gene_ids=(MATRIX/'row_ids.txt').read_text().splitlines()
    assert Z.shape==(len(ids),len(gene_ids))==(36523,12328) and len(positions)==len(ids)
    known=pd.read_parquet(PREPARED/'signatures.parquet')
    names=json.loads((PREPARED/'drugs.json').read_text())
    gene_info=pd.read_csv(META/'geneinfo_beta.txt',sep='\t')
    gene_map=dict(zip(gene_info.gene_id.astype(str),gene_info.ensembl_id));universe=set(gene_info.ensembl_id.dropna())
    sbg=pd.read_csv(GTEX/'signif.sbgenes.txt',sep='\t')
    sbg['gene_clean']=sbg.gene.str.split('.').str[0]
    union={t:set(x.gene_clean) for t,x in sbg.groupby('tissue')}
    split={f'{arm}_{t}':set(g.gene_clean[g.effsize.gt(0) if arm=='Female' else g.effsize.lt(0)]) for arm in ['Female','Male'] for t,g in sbg.groupby('tissue',sort=False)}
    assert len(union)==44 and len(split)==88
    ref=REPO/'Phase5_PerSignature_Categories/reference'
    assert union=={key:set(value) for key,value in json.loads((ref/'sb_library_sets.json').read_text()).items()}
    map5=json.loads((ref/'cmap_to_gtex.json').read_text())
    s=notebook_code(4,47)
    statement=s[s.index('cmap_to_gtex = {'):s.index('# Universe + translator')]
    env={'sbg':sbg};exec(statement,env);simplemap=env['cmap_to_gtex']

def expr(ids_):
    # cmapPy indexes the GCTX by source column position. sig50 is cached in that order.
    ix=sorted(positions[s] for s in ids_)
    # cmapPy 4.0.1 selects the full row list as a second advanced index and then transposes.
    # Its result has a contiguous signature axis, which affects float32 reduction order.
    subset=np.asarray(Z[ix])
    original_layout=subset[:,range(len(gene_ids))].T
    assert original_layout.flags.c_contiguous
    return pd.DataFrame(original_layout,index=pd.Index(gene_ids,name='rid'),columns=pd.Index([ids[i] for i in ix],name='cid'))

def mwu(expr_, sig):
    sex=sig.set_index('sig_id').donor_sex
    m=expr_.loc[:,sex.loc[expr_.columns].values=='M'].values
    f=expr_.loc[:,sex.loc[expr_.columns].values=='F'].values
    p=np.concatenate([mannwhitneyu(m[k:k+128],f[k:k+128],axis=1,alternative='two-sided').pvalue for k in range(0,len(m),128)])
    # The notebook computes a one-dimensional float32 mean separately for each gene.
    d=np.array([x.mean()-y.mean() for x,y in zip(m,f)])
    return pd.DataFrame({'gene_id':gene_ids,'difference':d,'p_value':p,'p_adj':multipletests(p,method='fdr_bh')[1]})

def gsea(rank,sets,perms):
    # GSEApy defaults to a multi-threaded Rust backend. Limit per-worker parallelism.
    gs=gseapy.prerank(rnk=rank.reset_index(),gene_sets=sets,outdir=None,permutation_num=perms,seed=42,threads=1,verbose=False,no_plot=True,min_size=15,max_size=10000)
    return gs.res2d

def minimum_sets(rank):
    available=set(rank.index)
    subsets={n:list(v & available) for n,v in union.items()}
    return {n:v for n,v in subsets.items() if len(v)>=15}

def phase2(drug,sub,gene_tests=None):
    de=mwu(expr(sub.sig_id),sub)
    de['ensembl_id']=de.gene_id.map(gene_map)
    if gene_tests is not None:
        gene_tests.parent.mkdir(parents=True, exist_ok=True)
        tmp = gene_tests.with_name(f'{gene_tests.stem}.{os.getpid()}.part')
        de.assign(drug=drug).to_csv(tmp,index=False,compression='gzip')
        tmp.replace(gene_tests)
    strong=set(de.ensembl_id[(de.p_adj<.05)&(de.difference.abs()>.5)].dropna())
    home=set(itertools.chain.from_iterable(simplemap.get(l,[]) for l in sub.cell_lineage.dropna().unique()))
    rows=[]
    for t,s in union.items():
        comparable=s&universe
        if len(comparable)<50:continue
        a=len(strong&comparable);b=len(strong)-a;c=len(comparable)-a
        pv=fisher_exact([[a,b],[c,len(universe)-a-b-c]],alternative='greater').pvalue
        rows.append(dict(drug=drug,tissue=t,overlap=a,p_value=pv,fold=a*len(universe)/(len(strong)*len(comparable)) if strong else 0,in_lineage=t in home,n_sigs=len(sub),n_strong_de=len(strong)))
    df=pd.DataFrame(rows);sel=df.overlap>0
    df['q_all']=multipletests(df.p_value,method='fdr_bh')[1]
    df['q_original']=np.nan
    if sel.any():df.loc[sel,'q_original']=multipletests(df.loc[sel,'p_value'],method='fdr_bh')[1]
    return df

def phase3(drug,sub):
    ans=[]
    for lin,g in sub.groupby('cell_lineage',sort=True):
        if lin=='unknown' or g.groupby('donor_sex').cell_iname.nunique().reindex(['F','M']).fillna(0).min()<5:continue
        if g.donor_sex.value_counts().reindex(['F','M']).fillna(0).min()<10:continue
        de=mwu(expr(g.sig_id),g)
        de['score']=np.sign(de.difference)*(-np.log10(de.p_value+1e-300))
        de['ensembl_id']=de.gene_id.map(gene_map)
        de=de.sort_values('score',ascending=False).dropna(subset=['ensembl_id','score']).drop_duplicates('ensembl_id')
        rank=de.set_index('ensembl_id').score.sort_values(ascending=False)
        result=gsea(rank,minimum_sets(rank),1000)
        for _,r in result.iterrows():
            ans.append(dict(drug=drug,lineage=lin,tissue=r.Term,NES=r.NES,FDR=r['FDR q-val'],is_home_tissue=r.Term in simplemap.get(lin,[])))
    return pd.DataFrame(ans,columns=['drug','lineage','tissue','NES','FDR','is_home_tissue'])

def phase4(drug,sub):
    ans=[]
    for lin,g in sub.groupby('cell_lineage',sort=True):
        if lin=='unknown' or len(g)<10: continue
        scores=expr(g.sig_id).mean(axis=1)
        ensg={gene_map[gid]:v for gid,v in scores.items() if pd.notna(gene_map.get(gid)) and gene_map.get(gid)}
        rank=pd.Series(ensg).sort_values(ascending=False)
        # The original counted |z| > 2 but ranked all mapped scores. The stored 9,260-row table came from this cell.
        result=gsea(rank,split,200)
        for _,r in result.iterrows():
            arm,t=r.Term.split('_',1)
            ans.append(dict(drug=drug,lineage=lin,tissue=t,sex_direction=arm,NES=r.NES,FDR=r['FDR q-val'],n_sigs=len(g),is_home_tissue=t in simplemap.get(lin,[])))
    return pd.DataFrame(ans)

def phase5(drug,sub):
    ans=[]
    for lin,g in sub.groupby('cell_lineage',sort=True):
        if lin=='unknown':continue
        ns=g.groupby('donor_sex').cell_iname.nunique().reindex(['M','F']).fillna(0)
        result={};reason='few_M_or_F'
        if ns.min()>=5 and len(g)>=10:
            scores=expr(g.sig_id).mean(axis=1);scores=scores[scores.abs()>2]
            de=pd.DataFrame({'gene_id':scores.index,'score':scores.values})
            de['ensembl_id']=de.gene_id.map(gene_map)
            de=de.sort_values('score',ascending=False).dropna(subset=['ensembl_id','score']).drop_duplicates('ensembl_id')
            rank=de.set_index('ensembl_id').score.sort_values(ascending=False)
            sets=minimum_sets(rank);reason='gsea_empty'
            if sets:
                result=gsea(rank,sets,1000).set_index('Term').to_dict('index');reason='gene_set_too_small'
        elif ns.min()>=5: reason='few_signatures'
        for t in sorted(union):
            h=t in map5.get(lin,[])
            row=dict(drug=drug,lineage=lin,tissue=t,n_M=int(ns.M),n_F=int(ns.F),is_home_tissue=h,category='D',skip_reason=reason,NES=np.nan,FDR=np.nan)
            if t in result:
                r=result[t];row.update(category=('A' if h else 'B') if r['FDR q-val']<.05 else 'C',skip_reason='',NES=r['NES'],FDR=r['FDR q-val'])
            ans.append(row)
    return pd.DataFrame(ans)

def work(job):
    phase,i=job;drug=names[i];d=Path(os.environ.get('P59_OUTPUT',BASE/'outputs'))/f'phase{phase}';d.mkdir(exist_ok=True,parents=True)
    path=d/f'{i:02}.csv.gz'
    gene_tests = None
    if phase == 2 and 'P59_GENE_TESTS' in os.environ:
        gene_tests = Path(os.environ['P59_GENE_TESTS']) / path.name
    if path.exists() and (gene_tests is None or gene_tests.exists()):
        return
    sub = known[known.cmap_name==drug]
    res = phase2(drug,sub,gene_tests) if phase == 2 else {3:phase3,4:phase4,5:phase5}[phase](drug,sub)
    tmp=d/f'{i:02}.{os.getpid()}.part';res.to_csv(tmp,index=False,compression='gzip');tmp.replace(path)
    print(f'Phase {phase} {i+1}/50 {drug} {len(res)} rows',flush=True)

if __name__=='__main__':
    import multiprocessing as mp
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prepare',action='store_true');p.add_argument('--phases',type=int,nargs='+');p.add_argument('--workers',type=int,default=1);p.add_argument('--first',type=int,default=1);p.add_argument('--last',type=int,default=50);a=p.parse_args()
    if a.prepare:prepare()
    if a.phases:
        print('versions',np.__version__,pd.__version__,scipy.__version__,gseapy.__version__,flush=True)
        assert 1<=a.first<=a.last<=50
        jobs=[(phase,i) for phase in a.phases for i in range(a.first-1,a.last)]
        with mp.get_context('spawn').Pool(a.workers,initializer=init) as pool:
            for _ in pool.imap_unordered(work,jobs):pass
