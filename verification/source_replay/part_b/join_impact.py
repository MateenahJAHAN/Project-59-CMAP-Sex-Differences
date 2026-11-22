"""Compare every model design and rerun the three actual ploidy benchmark tests using rebuilt covariates."""
import json, os, sys
from pathlib import Path
from types import SimpleNamespace

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = '1'
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
sys.path.insert(0, str(archive/'src/G'))
import gcore, gconfig, ginputs, gsets, ganalyses, gsensitivity
from run_g import PanelMeta

cfg=gconfig.load(str(archive/'config/G_defaults.json'))
data=ginputs.Inputs(str(archive),cfg,{'tensors_D':str(root/'run/D/data/D/tensors')},strict=True)
refdir=archive/'data/G/inputs/taskD'
freshdir=root/'run/D_labels/results/D'
cov=[gcore.read_covariates(str(p/'D_line_covariates.tsv')) for p in (freshdir,refdir)]
labels=[pd.read_csv(p/'D_line_table.tsv',sep='\t') for p in (freshdir,refdir)]
designs=[]; controls=[]
prim=gsets.primary_settings(cfg); B=cfg.run('B_full')
saved=pd.read_csv(archive/'results/G/G_controls.tsv',sep='\t',float_precision='round_trip')
exposure='ploidy within-lineage median split'
for panel in cfg.const('panels'):
    base=[PanelMeta(panel,data,c) for c in cov]
    aligns=[gsets.Aligned(m,l) for m,l in zip(base,labels)]
    for term in ('factY','factX','annotated'):
        pairs=[al.design(term,prim) for al in aligns]
        variations=[('primary',prim)]
        if term in ('factY','factX'):
            variations += [(label+' '+desc,opts) for label,desc,opts in gsensitivity.variants(cfg,prim,pairs[0]['mask'])
                           if (panel=='P1' or (panel=='P2' and label=='S2')) and (term=='factY' or not opts.get('y_only'))]
        for name,opts in variations:
            a,b=[al.design(term,opts) for al in aligns]
            same=np.array_equal(a['mask'],b['mask']) and np.array_equal(a['lin'],b['lin'])
            # Off-set covariates excluded by the design mask cannot enter any fit.
            xx=a['mask']
            xmatch=np.array_equal(a['Xo'][xx],b['Xo'][xx],equal_nan=True)
            designs.append(dict(panel=panel,term=term,variant=name,n_lines=int(xx.sum()),mask_lineage_equal=same,matrix_equal=xmatch))
            assert same and xmatch,(panel,term,name)
    a,b=[al.design('factY',prim) for al in aligns]
    median=[ganalyses.ploidy_split(p,d['mask'],d['lin']) for p,d in zip(base,(a,b))]
    assert np.array_equal(median[0],median[1],equal_nan=True)
    if panel=='P1':
        for term in ('factY','factX'):
            mask=aligns[0].design(term,prim)['mask']
            assert np.array_equal(base[0].cov.ploidy[mask].to_numpy(),base[1].cov.ploidy[mask].to_numpy(),equal_nan=True)
    p=data.panel(panel,cov[0],B)
    tp,_=gsets.tp53_columns(p.cov,gsets.tp53_code(cfg))
    ctx=SimpleNamespace(P={panel:p},cfg=cfg)
    row=ganalyses.exposure_run(ctx,panel,exposure,np.c_[tp,median[0]],tp.shape[1],[tp.shape[1]],a['mask'],a['lin'],prim['rule'],B)
    frozen=saved[(saved.panel==panel)&(saved.exposure==exposure)].iloc[0]
    differences={}
    for field,value in row.items():
        if field not in frozen.index or not isinstance(value,(float,int,np.number)): continue
        other=frozen[field]
        if not (value==other or (pd.isna(value) and pd.isna(other))):differences[field]=[value,other]
    controls.append(dict(panel=panel,median_split_equal=True,observation=row,differences=differences))
    print(panel,'ploidy benchmark',row['panel_p'],differences,flush=True)
    del p,ctx
out=args.output.resolve()
out.mkdir(parents=True,exist_ok=True)
pd.DataFrame(designs).to_csv(out/'model_design_comparison.tsv',sep='\t',index=False)
(out/'ploidy_control_comparison.json').write_text(json.dumps(controls,indent=2,default=lambda x:x.item())+'\n')
assert not any(r['differences'] for r in controls), controls
