"""Replay the Phase 4 original function on the largest early discrepancy."""
import ast,json,tempfile
from pathlib import Path
import numpy as np,pandas as pd,h5py,gseapy
from cmapPy.pandasGEXpress.parse import parse
import replay as p

p.init()
out=p.BASE/'phase4_literal';out.mkdir(exist_ok=True)
drug='NVP-BEZ235';lin='placenta';term='Female_Minor_Salivary_Gland'
sub=p.known[p.known.cmap_name==drug];pair=sub[sub.cell_lineage==lin]
ix=sorted(p.positions[s] for s in sub.sig_id);options={}
class Capture:
    def prerank(self,**kw):
        options.update(kw)
        return gseapy.prerank(**kw)
with tempfile.TemporaryDirectory() as tmp:
    name=Path(tmp)/'drug.gctx'
    with h5py.File(name,'w') as f:
        f.attrs['version']='GCTX1.0'; f.create_dataset('0/DATA/0/matrix',data=p.Z[ix])
        f.create_dataset('0/META/ROW/id',data=np.array(p.gene_ids,dtype='S'))
        f.create_dataset('0/META/COL/id',data=np.array([p.ids[i] for i in ix],dtype='S'))
    env=dict(np=np,pd=pd,gseapy=Capture(),gctx_file=str(name),gctx_sigs=set(sub.sig_id),sig_known=p.known,parse=parse,gid_to_ensg=p.gene_map,sb_libs_split=p.split)
    src=p.notebook_code(4,49);defs=[ast.get_source_segment(src,n) for n in ast.parse(src).body if isinstance(n,ast.FunctionDef)]
    (out/'original_functions.py').write_text('\n\n'.join(defs)+'\n')
    for code in defs:exec(compile(code,'Phase4_cell49','exec'),env)
    real=parse(str(name),cid=list(pair.sig_id)).data_df;fast=p.expr(pair.sig_id)
    assert np.array_equal(real.values,fast.values) and real.values.strides==fast.values.strides
    res=env['analyze_drug_lineage_v4'](drug,lin)
    li=res['gsea']
    cols=['Term','ES','NES','NOM p-val','FDR q-val','FWER p-val']
    tests=[li[cols].assign(scheme='literal_seed42_200')]
    t=pd.read_csv(p.BASE/'outputs/phase4'/f'{p.names.index(drug):02}.csv.gz')
    t=t[t.lineage==lin].copy();t['Term']=t.sex_direction+'_'+t.tissue
    m=t.merge(li[['Term','NES','FDR q-val']],on='Term',suffixes=('_fast','_literal'))
    assert len(m)==88 and ((m.NES_fast-m.NES_literal).abs()<1e-12).all() and ((m.FDR-m['FDR q-val']).abs()<1e-12).all()
    rank=options['rnk'].copy();value=rank.columns[-1]
    keys=['drug','lineage','tissue','sex_direction']
    a,b=[pd.read_csv(x).sort_values(keys).reset_index(drop=True) for x in sorted(next(p.REPO.glob('Phase4_*')).glob('*.csv'))]
    assert a.equals(b)
    hist=a[(a.drug==drug)&(a.lineage==lin)&(a.category!='C')].copy();hist['Term']=hist.sex_direction+'_'+hist.tissue
    comp=t.merge(hist[['Term','NES','FDR']],on='Term',how='left',suffixes=('_replay','_archive'))
    comp.to_csv(out/'archived_and_replay.csv',index=False)
    experiments=[('seed'+str(s)+'_'+str(n),s,n,rank) for s,n in [(11,200),(29,200),(47,200),(42,1000),(42,2000)]]
    for s in [11,29]:
        re=rank.sample(frac=1,random_state=s).sort_values(value,ascending=False,kind='stable').reset_index(drop=True)
        assert np.array_equal(re[value].to_numpy(),rank[value].to_numpy())
        experiments.append(('tieorder'+str(s)+'_seed42_200',42,200,re))
    for label,seed,n,rnk in experiments:
        gs=gseapy.prerank(**{**options,'rnk':rnk,'seed':seed,'permutation_num':n,'threads':1,'outdir':None}).res2d
        tests.append(gs[cols].assign(scheme=label));print(label,gs[gs.Term==term][['NES','FDR q-val']].to_dict('records'),flush=True)
    x=pd.concat(tests,ignore_index=True);x.to_csv(out/'all_seed_tests.csv',index=False)
    x[x.Term==term].to_csv(out/'target_seed_tests.csv',index=False)
    old=hist[hist.Term==term].iloc[0]
    summary={'example_selection':'Largest archived-to-replay NES difference among the first 33 completed compounds (archived significant cells).','drug':drug,'lineage':lin,'n_sigs':len(pair),'n_male_lines':res['n_M'],'n_female_lines':res['n_F'],'target':term,'archived_NES':old.NES,'archived_FDR':old.FDR,'both_archived_exports_equal_after_sort':True,'archived_row_count':len(a),'archived_significant_sets_this_pair':len(hist),'replay_significant_sets_this_pair':int((t.FDR<.05).sum()),'ranked_genes':len(rank),'genes_in_exact_ties':int(rank[value].duplicated(False).sum()),'max_abs_literal_vs_portable_NES':float((m.NES_fast-m.NES_literal).abs().max()),'max_abs_literal_vs_portable_FDR':float((m.FDR-m['FDR q-val']).abs().max()),'original_reader_matrix_and_strides_verified':True}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(summary,flush=True)
