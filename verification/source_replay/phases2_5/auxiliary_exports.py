"""Check the additional CSV and workbook tables that describe the historical phase results."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import compare

ROOT=Path(__file__).resolve().parent
dest=ROOT/'comparison'
receipt=[]

def check(label, archive, replay, keys, cols):
    merged=archive.merge(replay, how='outer', on=keys, suffixes=('_old','_new'), validate='one_to_one', indicator=True)
    paired=merged[merged._merge=='both']
    col_diffs={}
    for col in cols:
        a,b=paired[col+'_old'],paired[col+'_new']
        if pd.api.types.is_numeric_dtype(a.dtype) and pd.api.types.is_numeric_dtype(b.dtype):
            unequal=~np.isclose(a.to_numpy(dtype=float),b.to_numpy(dtype=float),atol=1e-12,rtol=0,equal_nan=True)
        else:
            unequal=a.astype('string').fillna('<missing>')!=b.astype('string').fillna('<missing>')
        col_diffs[col]=int(unequal.sum())
    receipt.append({'file_or_sheet':label, 'archived_rows':len(archive),'replayed_rows':len(replay),'archive_only_keys':int((merged._merge=='left_only').sum()),'replay_only_keys':int((merged._merge=='right_only').sum()),'mismatches_by_field':col_diffs})

def csv(folder,name):
    return pd.read_csv(folder/(name+'.csv'))

def phase2():
    folder=next(compare.REPO.glob('Phase2_*'))
    p=compare.outputs(2)
    gg=pd.read_csv(Path(__import__('os').environ['P59_META'])/'geneinfo_beta.txt',sep='\t')
    gene_ids=set(gg.ensembl_id.dropna())
    gtex=pd.read_csv(Path(__import__('os').environ['P59_GTEX'])/'signif.sbgenes.txt',sep='\t')
    gtex.gene=gtex.gene.str.split('.').str[0]
    set_size=gtex.groupby('tissue').gene.agg(lambda x:len(set(x)&gene_ids))
    a=p.query("drug=='AS-605240'").copy()
    a['tissue_genes_in_universe']=a.tissue.map(set_size)
    a['universe_genes']=len(gene_ids)
    a=a.rename(columns={'overlap':'overlap_with_two_genes','p_value':'fisher_p','q_original':'bh_nonzero_20','q_all':'bh_all_44'})
    orig=csv(folder,'as605240_tissue_tests')
    check('Phase2/as605240_tissue_tests.csv',orig,a,['tissue'],orig.columns.difference(['tissue']))
    o=csv(folder,'vorinostat_tissue_enrichment')
    v=p.query("drug=='vorinostat'").copy();v['tissue_size']=v.tissue.map(set_size)
    v['expected']=(v.n_strong_de*v.tissue_size/len(gene_ids)).round(1)
    v['fold']=v.fold.round(2);v['p_adj']=v.q_original
    check('Phase2/vorinostat_tissue_enrichment.csv',o,v,['tissue'],o.columns.difference(['tissue']))
    # Compare the lineage membership column separately; its archived order is arbitrary.
    a=csv(folder,'top50_drug_classification')
    meta=pd.read_parquet(ROOT/'prepared/signatures.parquet')
    flags=meta.groupby('cmap_name').cell_lineage.agg(lambda x:set(x.dropna())).to_dict()
    changed=[r.drug for r in a.itertuples() if set(r.drug_lineages.split('|'))!=flags[r.drug]]
    receipt.append({'file_or_sheet':'Phase2/top50_drug_classification.csv drug_lineages','compounds':len(a),'lineage_membership_changes':changed})
    (dest/'auxiliary_exports.json').write_text(json.dumps(receipt,indent=2)+'\n')

def phase5():
    folder=next(compare.REPO.glob('Phase5_*'))/'results'
    # GSEApy 1.3.1 exactly reproduces all stored Phase 5 FDR and NES values.
    f=pd.concat([pd.read_csv(path) for path in sorted((ROOT/'sensitivity_gseapy_1.3.1/phase5').glob('*.csv.gz'))],ignore_index=True)
    k=['drug','lineage','tissue'];p=['drug','lineage']
    f['sex_lean']=pd.Series(np.where(f.NES>0,'female','male')).where(f.NES.notna())
    sig=pd.read_parquet(ROOT/'prepared/signatures.parquet').rename(columns={'cmap_name':'drug','cell_lineage':'lineage','cell_iname':'cell_line','donor_sex':'sex'})
    meta=sig.groupby(p).agg(n_cell_lines=('cell_line','nunique'),n_signatures=('sig_id','nunique')).reset_index()
    for sex in ['M','F']:
        grouped=sig[sig.sex==sex].groupby(p).cell_line.agg(lambda x:', '.join(sorted(x.unique()))).rename('cell_lines_'+sex)
        meta=meta.merge(grouped,on=p,how='left')
        meta['cell_lines_'+sex]=meta['cell_lines_'+sex].fillna('—')
    full=f.merge(meta,on=p,validate='many_to_one')
    hit=full[full.category.isin(['A','B'])]
    n=f.groupby(p+['category']).size().unstack(fill_value=0).add_prefix('n_').reset_index()
    n=n.merge(full.sort_values(['category','FDR',*k]).groupby(p).first().reset_index()[p+['n_M','n_F','skip_reason','n_cell_lines','cell_lines_M','cell_lines_F','n_signatures']],on=p)
    n['passes_5M5F']=(n.n_M>=5)&(n.n_F>=5)
    n['analysed']=n.skip_reason.isna()
    n['n_hits']=n.n_A+n.n_B
    cnt=n.groupby('drug').agg(Cat_A_home=('n_A','sum'),Cat_B_cross=('n_B','sum'),Cat_C_tested_no_hit=('n_C','sum'),Cat_D_not_tested=('n_D','sum'),n_lineages=('lineage','nunique'),n_hits=('n_hits','sum')).reset_index()
    cnt['total']=cnt.filter(like='Cat_').sum(axis=1)
    checklist=[
      ('per_drug_4category_counts.csv',cnt,['drug']),('per_drug_lineage.csv',n,p),('significant_hits.csv',hit,k),
      ('all50_drug_counts.xlsx',cnt,['drug']),('4_category_breakdown_FIXED.xlsx',cnt,['drug']),
      ('FULL_50drugs_table.xlsx:Full_table',full,k),('FULL_50drugs_table.xlsx:Significant_hits',hit,k),
      ('FULL_50drugs_table.xlsx:Lineage_summary',n,p),('FULL_50drugs_table.xlsx:Per_drug_counts',cnt,['drug']),
      ('FULL_50drugs_table.xlsx:Cell_line_detail',sig.groupby(p+['cell_line','sex']).sig_id.nunique().reset_index(name='n_signatures'),p+['cell_line','sex'])]
    for fname,replayed,keys in checklist:
        if fname.endswith('.csv'):
            archive=pd.read_csv(folder/fname)
        else:
            name,sep,sheet=fname.partition(':');archive=pd.read_excel(folder/name,sheet_name=sheet if sep else 0)
        cols=[c for c in archive if c in replayed and c not in keys]
        print('Checking',fname,flush=True)
        check('Phase5/'+fname,archive,replayed,keys,cols)
    a=pd.read_excel(folder/'cat_a_enrichment_scores.xlsx')
    check('Phase5/cat_a_enrichment_scores.xlsx shared numeric results',a,hit[hit.category=='A'],k,['NES','FDR'])
    raw=pd.read_csv(Path(__import__('os').environ['P59_GTEX'])/'signif.sbgenes.txt',sep='\t',nrows=5)
    gt=pd.read_excel(folder/'gtex_sbgenes_first5.xlsx')
    check('Phase5/gtex_sbgenes_first5.xlsx',gt,raw,['gene','tissue'],['effsize','effsize_se','lfsr'])
    (dest/'auxiliary_exports.json').write_text(json.dumps(receipt,indent=2)+'\n')
    for r in receipt:
        m=r.get('mismatches_by_field',{})
        print(r['file_or_sheet'],sum(m.values()),{k:v for k,v in m.items() if v},flush=True)

if __name__=='__main__':
    phase2();phase5()
