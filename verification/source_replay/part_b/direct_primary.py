"""Primary Part B rerun explicitly reading the newly built labels, covariates, Level 5 and Level 3 arrays."""
import os,sys
for key in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","VECLIB_MAXIMUM_THREADS","NUMEXPR_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]='1'
from pathlib import Path
from types import SimpleNamespace
import json
import numpy as np
import pandas as pd

import argparse
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--engine', type=Path, required=True, help='Checkout of the original analysis engine with its versioned data inputs.')
parser.add_argument('--rebuild', type=Path, required=True, help='Rebuild directory containing run/D, run/E and run/D_labels.')
parser.add_argument('--output', type=Path, required=True, help='Directory for regenerated result tables.')
args=parser.parse_args()
root=args.rebuild.resolve()
archive=args.engine.resolve()
sys.path.insert(0,str(archive/'src/G'))
import gconfig,ginputs,gcore,gsets,ganalyses
out=args.output.resolve()
out.mkdir(parents=True,exist_ok=True)
cfg=gconfig.load(str(archive/'config/G_defaults.json'))
inputs=ginputs.Inputs(str(archive),cfg,{'tensors_D':str(root/'run/D/data/D/tensors'),'tensors_E':str(root/'run/E/data/E/tensors')},strict=True)
D=root/'run/D_labels/results/D'
cov=gcore.read_covariates(str(D/'D_line_covariates.tsv'))
lt=pd.read_csv(D/'D_line_table.tsv',sep='\t')
B=cfg.run('B_full')
activity=[]
for name in cfg.const('panels'):
    p=inputs.panel(name,cov,B)
    scores=[]
    for j,compound in enumerate(p.comps):
        Y=p.R[p.avail[:,j],j,:]
        scores.append((compound,float(np.corrcoef(Y)[np.triu_indices(len(Y),1)].mean()),len(Y)))
    tab=pd.DataFrame(scores,columns=['cmap_name','activity_mean_pairwise_r','n_lines'])
    tab['panel_median']=tab.activity_mean_pairwise_r.median()
    tab['above_panel_median']=tab.activity_mean_pairwise_r>tab.panel_median
    tab.insert(0,'panel',name)
    activity.append(tab)
    del p
act=pd.concat(activity,ignore_index=True)
act.to_csv(out/'D_activity_scores.tsv',sep='\t',index=False)
frozen_act=inputs.activity()
assert act[['panel','cmap_name','above_panel_median']].equals(frozen_act[['panel','cmap_name','above_panel_median']])
p=inputs.panel('P1',cov,B)
aligned=gsets.Aligned(p,lt)
act=act[act.panel=='P1']
ctx=SimpleNamespace(cfg=cfg,P={'P1':p},s_level5={'P1':gcore.response_size(p.R,p.avail)},activity={'P1':dict(zip(act.cmap_name,act.above_panel_median))})
ctx.auc={'P1':inputs.prism_auc(p,inputs.panel_definitions())[0]}
l3=inputs.level3('R3_P1_d111',p,order='D')[0]
pooled=np.nanmean(gcore.response_size(l3,~np.isnan(l3[:,:,0])),axis=1)
old=set(inputs.old_list().cmap_name)
h1=[];h3=[];cpds=[]
for term in ('factY','factX'):
    d=aligned.design(term,gsets.primary_settings(cfg))
    ans=ganalyses.run_main_term(ctx,'P1',term,d,B,keep_perm=True)
    r4=ganalyses.rule4(ctx,'P1',ans)
    compounds=ganalyses.apply_findings(ctx,'P1',ans['per_cpd'],r4,old)
    e3=ganalyses.e3_row(ctx,'P1',term,d,pooled,'level3','R3_P1_d111','primary')
    h1.append(ans['panel_row']);h3.append(e3);cpds.append(compounds)
    print(term,'panel',h1[-1]['panel_p'],'DMSO',h1[-1]['dmso_p'],'estimate',e3['estimate'],'90% interval',e3['lo90'],e3['hi90'],'q minimum',compounds['q'].min(),'findings',compounds.finding.sum(),flush=True)
checks={}
for name,a,template,ids in [('H1',pd.DataFrame(h1),'G_H1.tsv',['panel','term','analysis']),('H3',pd.DataFrame(h3),'G_H3.tsv',['panel','term','scale','role']),('per_compound',pd.concat(cpds,ignore_index=True),'G_H2_per_compound.tsv',['panel','term','analysis','compound'])]:
    a.to_csv(out/(name+'.tsv'),sep='\t',index=False)
    # Compare stored serialized values, matching the published method for each field.
    a=pd.read_csv(out/(name+'.tsv'),sep='\t',float_precision='round_trip').set_index(ids).sort_index()
    ref=pd.read_csv(archive/'results/G'/template,sep='\t',float_precision='round_trip').set_index(ids).loc[a.index,a.columns].sort_index()
    differences={}
    for k in a:
        eq=(a[k]==ref[k])|(a[k].isna()&ref[k].isna())
        if not eq.all():differences[k]=int((~eq).sum())
    checks[name]={'rows':len(a),'fields':len(a.columns),'diffs':differences}
(out/'verification.json').write_text(json.dumps(checks,indent=2)+'\n')
assert not any(r['diffs'] for r in checks.values()),checks
