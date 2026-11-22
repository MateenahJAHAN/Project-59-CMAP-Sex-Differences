"""Compare the older Phase 5 workbooks at their own drug-level units."""
from collections import Counter
from pathlib import Path
import json
import numpy as np
import pandas as pd
import auxiliary_exports as audit
import compare

root=Path(__file__).resolve().parent
folder=next(compare.REPO.glob('Phase5_*'))/'results'
fp=root/'comparison/phase5_working_exports.json'
f=pd.concat([pd.read_csv(p) for p in sorted((root/'sensitivity_gseapy_1.3.1/phase5').glob('*.csv.gz'))],ignore_index=True)
home=json.loads((folder.parent/'reference/cmap_to_gtex.json').read_text())
k=['drug','tissue']; hit=f[f.category.isin(['A','B'])].copy()
drug_cat=f.groupby('drug').category.min().replace({'D':'C'})
f['compound_category']=f.drug.map(drug_cat)
c=f.copy(); c['NES']=c.NES.round(3);c['FDR']=c.FDR.round(4)
c['abs_score']=c.NES.abs()
c=c.sort_values('abs_score',ascending=False).drop_duplicates(k).drop(columns='category').rename(columns={'compound_category':'category'})
old=pd.read_excel(folder/'50_drugs_44_tissues_FULL.xlsx')
assert len(c)==len(old)==2200
audit.check('50_drugs_44_tissues_FULL.xlsx',old,c,k,['NES','FDR','category'])
n=f.groupby('drug').lineage.agg(lambda x:set(t for lin in x.unique() for t in home.get(lin,[]))).to_dict()
c['tissue_category']=[('D' if pd.isna(row.FDR) else 'C' if row.FDR>=.05 else 'A' if row.tissue in n[row.drug] else 'B') for row in c.itertuples()]
cnt=c.groupby(['drug','tissue_category']).size().unstack(fill_value=0).rename(columns={'A':'Cat_A_home','B':'Cat_B_cross','C':'Cat_C_tested_no_hit','D':'Cat_D_not_tested'})
cnt['total']=cnt.sum(axis=1);cnt['activity']=cnt.Cat_A_home+cnt.Cat_B_cross
old=pd.read_excel(folder/'4_category_breakdown.xlsx')
audit.check('4_category_breakdown.xlsx',old,cnt.reset_index(),['drug'],list(cnt.columns))
old=pd.read_excel(folder/'50_drugs_complete_table.xlsx')
cols=['drug','tissue','NES','FDR','category']
hit['NES']=hit.NES.round(3);hit['FDR']=hit.FDR.round(4);hit['category']=hit.drug.map(drug_cat)
missing=[{'drug':drug,'tissue':'No significant hit','NES':'—','FDR':'—','category':cat} for drug,cat in drug_cat.items() if drug not in set(hit.drug)]
complete=pd.concat([hit[cols],pd.DataFrame(missing)],ignore_index=True)
olds=Counter(old[cols].itertuples(index=False,name=None));news=Counter(complete[cols].itertuples(index=False,name=None))
observed=lambda h:[{'row':x,'count':cnt} for x,cnt in h.items()]
audit.receipt.append({'file_or_sheet':'50_drugs_complete_table.xlsx','archived_rows':len(old),'replayed_rows':len(complete),'archive_only':observed(olds-news),'replay_only':observed(news-olds)})
fp.write_text(json.dumps(audit.receipt,indent=2)+'\n')
for r in audit.receipt:
 print(r['file_or_sheet'],r.get('mismatches_by_field',{}),len(r.get('archive_only',[])),len(r.get('replay_only',[])),flush=True)
