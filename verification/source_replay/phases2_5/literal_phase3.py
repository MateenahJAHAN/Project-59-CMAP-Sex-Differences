"""Re-execute the original three Phase 3 functions for the two compounds whose categories changed, including the original GCTX parser."""
from pathlib import Path
import argparse, ast, json, re, tempfile
import numpy as np
import pandas as pd
import scipy, gseapy, h5py
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from cmapPy.pandasGEXpress.parse import parse
import replay as p

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--drug', choices=['fostamatinib','DL-PDMP'], default='fostamatinib')
drug = parser.parse_args().drug
p.init()
out = p.BASE / ('phase3_literal' if drug == 'fostamatinib' else 'phase3_dl_pdmp')
out.mkdir(exist_ok=True)
sub=p.known[p.known.cmap_name==drug]
ix=sorted(p.positions[s] for s in sub.sig_id)
with tempfile.TemporaryDirectory() as tmp:
    gfile=Path(tmp)/'subset.gctx'
    with h5py.File(gfile,'w') as f:
        f.attrs['version']='GCTX1.0';f.create_dataset('0/DATA/0/matrix',data=p.Z[ix])
        f.create_dataset('0/META/ROW/id',data=np.array(p.gene_ids,dtype='S'))
        f.create_dataset('0/META/COL/id',data=np.array([p.ids[i] for i in ix],dtype='S'))
    env=dict(np=np,pd=pd,gseapy=gseapy,multipletests=multipletests,mannwhitneyu=mannwhitneyu,parse=parse,
             gctx_file=str(gfile),gctx_sigs=set(sub.sig_id),sig_known=p.known,
             genes_info=p.gene_info,gid_to_ensg=p.gene_map,sb_library_sets=p.union,cmap_to_gtex=p.simplemap)
    src={}
    for ci in [29,31,33]:
        orig=p.notebook_code(3,ci)
        func=next(v for v in ast.parse(orig).body if isinstance(v,ast.FunctionDef))
        code=ast.get_source_segment(orig,func);src[func.name]=code
        exec(compile(code,'phase3_notebook_cell_'+str(ci),'exec'),env)
    (out/'notebook_functions.py').write_text('\n\n'.join(src.values())+'\n')
    res=env['analyze_drug_full'](drug,verbose=False)
    counts=[];tables=[];ctx=[]
    for lin,g in sub.groupby('cell_lineage'):
        groups=g.groupby('donor_sex').agg(lines=('cell_iname','nunique'),signatures=('sig_id','size')).reindex(['F','M']).fillna(0)
        tested=groups.lines.min()>=5 and groups.signatures.min()>=10 and len(g)>=20 and lin!='unknown'
        counts.append(dict(lineage=lin,female_lines=int(groups.loc['F','lines']),male_lines=int(groups.loc['M','lines']),female_signatures=int(groups.loc['F','signatures']),male_signatures=int(groups.loc['M','signatures']),eligible=tested,home_tissues='|'.join(sorted(p.simplemap.get(lin,[])))))
        assert tested==(lin in res['lineage_results'])
        if tested:
            real=parse(str(gfile),cid=list(g.sig_id)).data_df
            surrogate=p.expr(g.sig_id)
            assert real.index.tolist()==surrogate.index.tolist() and real.columns.tolist()==surrogate.columns.tolist()
            assert real.values.flags.c_contiguous==surrogate.values.flags.c_contiguous
            assert np.array_equal(real.values,surrogate.values)
    pd.DataFrame(counts).to_csv(out/'lineage_eligibility.csv',index=False)
    print('literal',res['category'],res['cat_a_evidence'],res['cat_b_evidence'],flush=True)
    target_lines={h[0] for h in res['cat_a_evidence']}
    # If literal did not give the replayed home hit, still test its affected lineage.
    full=pd.read_csv(p.BASE/'outputs/phase3'/f'{p.names.index(drug):02}.csv.gz')
    target_lines.update(full.loc[full.is_home_tissue & (full.FDR<.05),'lineage'])
    if drug == 'DL-PDMP':
        target_lines.update(full.loc[full.FDR < .05, 'lineage'])
    for lin,r in res['lineage_results'].items():
        frame=r['gsea'][['Term','NES','FDR q-val']].copy()
        frame['lineage']=lin;frame['scheme']='original_notebook';tables.append(frame)
        if lin not in target_lines:continue
        # Capture exactly the ranked values and sets supplied by the literal notebook.
        vals={}
        class Tap:
            def prerank(self,**kwargs):
                vals.update(kwargs)
                return type('Result',(),{'res2d':r['gsea']})()
        env['gseapy']=Tap();env['run_gsea'](r['de'],p.union)
        rank=vals['rnk'].copy()
        ctx.append({'lineage':lin,'gene_count':len(rank),'count_in_tied_scores':int(rank.iloc[:,1].duplicated(False).sum()),'distinct_scores':rank.iloc[:,1].nunique(),'n_eligible_sets':len(vals['gene_sets'])})
        for seed in [11,29,47]:
            rr=rank.sample(frac=1,random_state=seed).sort_values(rank.columns[1],ascending=False,kind='stable').reset_index(drop=True)
            assert np.array_equal(rr.iloc[:,1].to_numpy(),rank.iloc[:,1].to_numpy())
            kw={**vals,'rnk':rr,'threads':1}
            gs=gseapy.prerank(**kw).res2d[['Term','NES','FDR q-val']].copy()
            gs['lineage']=lin;gs['scheme']='tied_order_'+str(seed);tables.append(gs)
            print('tie run',lin,seed,flush=True)
    allgs=pd.concat(tables,ignore_index=True)
    allgs.to_csv(out/'gsea_all_tissues.csv',index=False)
    original=allgs[allgs.scheme=='original_notebook'].rename(columns={'Term':'tissue','NES':'NES_literal','FDR q-val':'FDR_literal'})
    match=full.merge(original[['lineage','tissue','NES_literal','FDR_literal']],on=['lineage','tissue'],validate='one_to_one',how='outer',indicator=True)
    assert all(match._merge=='both')
    a=pd.read_csv(next(p.REPO.glob('Phase3*/phase3_classification_final.csv'))); row=a[a.drug==drug].iloc[0]
    old=[]
    for col in ['cat_a_top','cat_b_top']:
        for lin,t,nes,fdr in ast.literal_eval(re.sub(r'np.float64\(([^)]+)\)',r'\1',row[col])):
            old.append(dict(lineage=lin,tissue=t,archived_NES=nes,archived_FDR=fdr))
    match=match.merge(pd.DataFrame(old, columns=['lineage','tissue','archived_NES','archived_FDR']),on=['lineage','tissue'],how='left',validate='one_to_one')
    match.drop(columns='_merge').to_csv(out/'archived_and_literal_scores.csv',index=False)
    asc=[]
    for h in allgs.itertuples(index=False,name=None):
        t,nes,fdr,lin,scheme=h
        if t in p.simplemap.get(lin,[]):
            asc.append({'scheme':scheme,'lineage':lin,'home_tissue':t,'NES':nes,'FDR':fdr,'FDR_below_.05':bool(fdr<.05)})
    pd.DataFrame(asc).to_csv(out/'home_tissue_scores.csv',index=False)
    summary={'drug':drug,'raw_signatures':len(sub),'passing_lineages':res['n_lineages_tested'],'literal_category':res['category'],'archived_category':row.category,'archive_recorded_A_B_counts':[int(row.cat_a_evidence_count),int(row.cat_b_evidence_count)],'archived_individual_numeric_pairs':len(old),'n_literal_tissue_results':len(original),'max_abs_NES_literal_vs_fast':float((match.NES-match.NES_literal).abs().max()),'max_abs_FDR_literal_vs_fast':float((match.FDR-match.FDR_literal).abs().max()),'tie_diagnostics':ctx,'home_calls_by_scheme':pd.DataFrame(asc).groupby('scheme')['FDR_below_.05'].sum().astype(int).to_dict(),'permutation_num':1000,'gsea_seed':42}
    if drug == 'DL-PDMP':
        target = allgs[(allgs.lineage == 'large_intestine') & (allgs.Term == 'Cells_EBV.transformed_lymphocytes')]
        target.to_csv(out/'target_tissue_scores.csv', index=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary),flush=True)
