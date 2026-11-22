"""Measure the change in ranks and calls for the first discrepant compound/lineage of the layout sensitivity pass."""
from pathlib import Path
import json
import replay as r
import pandas as pd
import numpy as np
r.init()
g=r.known[(r.known.cmap_name=='tozasertib')&(r.known.cell_lineage=='breast')]
ix=sorted(r.positions[s] for s in g.sig_id); mat=np.asarray(r.Z[ix])
ranks={};tables=[]
for name,layout in [('simple_transpose',mat.T),('original_cmapPy',mat[:,range(len(r.gene_ids))].T)]:
    e=pd.DataFrame(layout,index=r.gene_ids,columns=[r.ids[i] for i in ix])
    scores=e.mean(axis=1)
    ensg={r.gene_map[k]:v for k,v in scores.items() if pd.notna(r.gene_map.get(k)) and r.gene_map.get(k)}
    ranks[name]=pd.Series(ensg).sort_values(ascending=False)
    t=r.gsea(ranks[name],r.split,200);t['layout']=name;tables.append(t[['layout','Term','NES','FDR q-val']])
t=pd.concat(tables)
arc=pd.read_csv(next(r.REPO.glob('Phase4*/phase4_classifications_final.csv')))
arc=arc[(arc.drug=='tozasertib')&(arc.lineage=='breast')&(arc.category!='C')].copy()
arc['Term']=arc.sex_direction+'_'+arc.tissue
t=t.merge(arc[['Term','NES','FDR']],on='Term',how='left',suffixes=('_run','_archive'))
out=r.BASE/'sensitivity_array_layout'
t.to_csv(out/'tozasertib_breast_layout.csv',index=False)
x,y=ranks.values();z=x.index.intersection(y.index)
p={'reason_for_selection':'first completed compound from the first layout run, and its largest observed false-discovery-rate discrepancy','signatures':len(g),'unique_gene_ids':len(z),'genes_with_different_float32_means':int((x[z]!=y[z]).sum()),'max_mean_difference':float((x[z]-y[z]).abs().max()),'different_rank_positions':int((x.index!=y.index).sum()),'mean_z_tied_genes':int(y.duplicated(keep=False).sum())}
for name,table in t.groupby('layout'):
    a=set(arc.Term);b=set(table.loc[table['FDR q-val']<.05,'Term'])
    p[name]={'lost_vs_archive':sorted(a-b),'gained_vs_archive':sorted(b-a),'max_shared_NES_error':float((table.NES_run-table.NES_archive).abs().max())}
(out/'layout_check.json').write_text(json.dumps(p,indent=2)+'\n');print(json.dumps(p))
