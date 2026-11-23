"""Compare keyed historical exports to each completed source-data replay."""
from pathlib import Path
import ast, json, os, re, sys
import numpy as np
import pandas as pd

BASE=Path(__file__).resolve().parent
REPO=Path(os.environ.get('P59_PUBLIC', BASE.parents[2]))
PREPARED=Path(os.environ.get('P59_PREPARED', BASE/'prepared'))
OUT=Path(os.environ.get('P59_COMPARISON',BASE/'comparison'))
OUT.mkdir(exist_ok=True,parents=True)
names=json.loads((PREPARED/'drugs.json').read_text())
def outputs(n):
    files=sorted((Path(os.environ.get('P59_OUTPUT',BASE/'outputs'))/f'phase{n}').glob('*.csv.gz'))
    if len(files)!=50: raise RuntimeError(f'Phase {n}: {len(files)} / 50 files')
    rows = [pd.read_csv(f) for f in files]
    return pd.concat([table for table in rows if not table.empty], ignore_index=True)
def save(x,name):
    x.to_csv(OUT/(name+'.csv'),index=False)
def counts(x):
    return {str(k):int(v) for k,v in x.value_counts().items()}
def dif(x,y):
    a=pd.to_numeric(x); b=pd.to_numeric(y)
    return {'compared':int((a.notna()&b.notna()).sum()),'missingness_changed':int((a.isna()!=b.isna()).sum()),'above_1e-8':int(((a-b).abs()>1e-8).sum()),'max_abs':float((a-b).abs().max()) if (a-b).notna().any() else None}
def compare2():
    t=outputs(2);arc=pd.read_csv(next(REPO.glob('Phase2*/top50_drug_classification.csv'))).set_index('drug')
    rows=[]
    for drug,g in t.groupby('drug'):
        r={'drug':drug,'n_strong_de':int(g.n_strong_de.iloc[0]),'n_sigs':int(g.n_sigs.iloc[0])}
        for q in ['q_original','q_all']:
            sig=g[g[q]<.05].sort_values('p_value');home=sig[sig.in_lineage]
            best=home.iloc[0] if len(home) else sig.iloc[0] if len(sig) else None
            r.update({q+'_call':'A' if len(home) else 'B' if len(sig) else 'C',q+'_best_tissue':best.tissue if best is not None else None,q+'_best_q':best[q] if best is not None else np.nan,q+'_fold':round(best.fold,2) if best is not None else np.nan})
        old=arc.loc[drug];r.update(archive_category=old.category,archive_n_strong=int(old.n_strong_de),archive_n_sigs=int(old.n_sigs))
        rows.append(r)
    ret=pd.DataFrame(rows);save(ret,'phase2_compounds')
    save(t,'phase2_all_tissue_tests')
    av=arc.reindex(ret.drug)
    assert set(ret.drug)==set(arc.index)
    return {'compounds':len(ret),'gene_tests':len(ret)*12328,'tissue_tests':len(t),'category_changes':int((ret.q_original_call!=ret.archive_category).sum()),'signature_count_changes':int((ret.n_sigs!=ret.archive_n_sigs).sum()),'DE_gene_count_changes':int((ret.n_strong_de!=ret.archive_n_strong).sum()),'original_calls':counts(ret.q_original_call),'all_44_calls':counts(ret.q_all_call),'best_q_comparison':dif(ret.q_original_best_q,av.best_p.reset_index(drop=True)),'best_fold_comparison':dif(ret.q_original_fold,av.best_fold.reset_index(drop=True)),'best_tissue_mismatches':int((ret.q_original_best_tissue.fillna('')!=av.best_tissue.reset_index(drop=True).fillna('')).sum())}
def compare3():
    t=outputs(3)
    arc=pd.read_csv(next(REPO.glob('Phase3*/phase3_classification_final.csv'))).set_index('drug')
    rows=[];hits=[]
    for drug in names:
        s=t[t.drug==drug];sig=s[s.FDR<.05];a=sig[sig.is_home_tissue];b=sig[~sig.is_home_tissue]
        old=arc.loc[drug]
        rows.append(dict(drug=drug,original=old.category,replay='A' if len(a) else 'B' if len(b) else 'C',original_lineages=old.n_lineages_tested,tested_lineages=s.lineage.nunique(),original_A=old.cat_a_evidence_count,replay_A=len(a),original_B=old.cat_b_evidence_count,replay_B=len(b)))
        for col in ['cat_a_top','cat_b_top']:
            for lin,ti,nes,fdr in ast.literal_eval(re.sub(r'np.float64\(([^)]+)\)',r'\1',old[col])):
                r=s[(s.lineage==lin)&(s.tissue==ti)].iloc[0]
                hits.append(dict(drug=drug,lineage=lin,tissue=ti,original_NES=nes,replay_NES=r.NES,original_FDR=fdr,replay_FDR=r.FDR))
    x=pd.DataFrame(rows);h=pd.DataFrame(hits)
    save(x,'phase3_compounds');save(h,'phase3_archived_leading_evidence')
    sig=t[t.FDR<.05].sort_values(['drug','lineage','FDR','tissue']).copy()
    sig['classification']=np.where(sig.is_home_tissue,'A','B')
    save(sig,'phase3_replayed_significant_tissues')
    return {'compounds':len(x),'pairs':int(x.tested_lineages.sum()),'compounds_with_no_eligible_lineage':int((x.tested_lineages==0).sum()),'tested_compounds_without_a_hit':int(((x.tested_lineages>0)&(x.replay=='C')).sum()),'tissue_results':len(t),'archived':counts(x.original),'replay':counts(x.replay),'original_total_evidence':int((x.original_A+x.original_B).sum()),'replay_total_evidence':int((x.replay_A+x.replay_B).sum()),'archived_individually_scored_tissue_results':len(h),'category_changes':int((x.original!=x.replay).sum()),'tested_lineage_changes':int((x.original_lineages!=x.tested_lineages).sum()),'compounds_with_different_evidence_counts':int(((x.original_A!=x.replay_A)|(x.original_B!=x.replay_B)).sum()),'archived_leading_NES':dif(h.original_NES,h.replay_NES),'archived_leading_FDR':dif(h.original_FDR,h.replay_FDR)}
def compare4():
    t=outputs(4);arc=pd.read_csv(next(REPO.glob('Phase4*/phase4_classifications_final.csv'))); a=arc[arc.category!='C'].copy()
    s=t[t.FDR<.05].copy(); s['category']=np.where(s.is_home_tissue,'A','B')
    key=['drug','lineage','tissue','sex_direction']
    orig=a.rename(columns={'category':'original_category','NES':'original_NES','FDR':'original_FDR'})
    every=t.merge(orig[key+['original_category','original_NES','original_FDR']],on=key,how='outer',validate='one_to_one',indicator=True)
    assert all(every._merge!='right_only')
    every['original_significant']=every.original_category.notna();every['replay_significant']=every.FDR<.05
    chg=every[every.original_significant!=every.replay_significant].copy();save(chg,'phase4_changed_threshold_calls')
    cur=s.groupby(['drug','lineage']).size();pai=t[['drug','lineage']].drop_duplicates()
    nc=len(pai)-len(cur)
    x=pai.merge(arc.groupby(['drug','lineage','category']).size().unstack(fill_value=0).reset_index(),on=['drug','lineage'],how='outer')
    rem=s.groupby(['drug','lineage','category']).size().unstack(fill_value=0).reset_index()
    for c in ['A','B']: rem=rem.rename(columns={c:'replay_'+c})
    x=x.merge(rem,how='left',on=['drug','lineage']).fillna(0);x['replay_C']=((x.replay_A+x.replay_B)==0).astype(int);save(x,'phase4_pair_counts')
    oldboth=every[every.original_significant]
    return {'pairs':len(pai),'all_tested_sets':len(t),'original_categories':counts(arc.category),'archived_pair_signature_count_changes':int((arc.groupby(['drug','lineage']).n_sigs.first().sort_index()!=t.groupby(['drug','lineage']).n_sigs.first().sort_index()).sum()),'replay_categories':{**counts(s.category),'C':nc},'archive_only':int((every.original_significant&~every.replay_significant).sum()),'replay_only':int((~every.original_significant&every.replay_significant).sum()),'no_hit_both':int(((x.C==1)&(x.replay_C==1)).sum()),'no_hit_archive_only':x.loc[(x.C==1)&(x.replay_C==0),['drug','lineage']].to_dict('records'),'no_hit_replay_only':x.loc[(x.C==0)&(x.replay_C==1),['drug','lineage']].to_dict('records'),'pairs_with_changed_counts':int(((x.A!=x.replay_A)|(x.B!=x.replay_B)|(x.C!=x.replay_C)).sum()),'original_signal_NES':dif(oldboth.original_NES,oldboth.NES),'original_signal_FDR':dif(oldboth.original_FDR,oldboth.FDR),'archive_lib_sex':counts(a.sex_direction),'replay_lib_sex':counts(s.sex_direction),'archived_lapatinib_lung':counts(a[(a.drug=='lapatinib')&(a.lineage=='lung')].sex_direction),'replayed_lapatinib_lung':counts(s[(s.drug=='lapatinib')&(s.lineage=='lung')].sex_direction)}
def compare5():
    t=outputs(5);arc=pd.read_csv(next(REPO.glob('Phase5*/results/per_signature_categories.csv')))
    ks=['drug','lineage','tissue'];v=t.merge(arc,on=ks,suffixes=('_replay','_old'),validate='one_to_one',how='outer',indicator=True)
    assert (v._merge=='both').all()
    diffs=v[v.category_old!=v.category_replay];save(diffs,'phase5_changed_calls')
    s=v[v.category_old!='D']; quant=s[ks+['FDR_old','FDR_replay','NES_old','NES_replay','category_old','category_replay']].copy();save(quant,'phase5_tested_comparison')
    return {'rows':len(v),'pairs':len(v[ks[:2]].drop_duplicates()),'archived':counts(v.category_old),'replay':counts(v.category_replay),'changed_categories':len(diffs),'changed_skip_reasons':int((v.skip_reason_old.fillna('')!=v.skip_reason_replay.fillna('')).sum()),'home_mapping_changes':int((v.is_home_tissue_replay!=v.is_home_tissue_old).sum()),'cell_line_count_changes':int(((v.n_M_replay!=v.n_M_old)|(v.n_F_replay!=v.n_F_old)).sum()),'NES':dif(v.NES_old,v.NES_replay),'FDR':dif(v.FDR_old,v.FDR_replay)}
if __name__=='__main__':
    for phase in (sys.argv[1:] or ['2','3','4','5']):
        result={2:compare2,3:compare3,4:compare4,5:compare5}[int(phase)]()
        (OUT/f'phase{phase}.json').write_text(json.dumps(result,indent=2)+'\n')
        print('Phase',phase,json.dumps(result))
