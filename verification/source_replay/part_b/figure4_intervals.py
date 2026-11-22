"""Recompute Figure 4D's 12 estimates and intervals from rebuilt expression and covariates."""
import argparse
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--engine',required=True,type=Path)
parser.add_argument('--rebuild',required=True,type=Path)
parser.add_argument('--output',required=True,type=Path)
args=parser.parse_args()
engine=args.engine.resolve()
rebuild=args.rebuild.resolve()
args.output.mkdir(exist_ok=True,parents=True)
sys.path.insert(0,str(engine/'src/G'))
import gconfig
import ginputs
import gcore
import gsets
import ganalyses

cfg=gconfig.load(str(engine/'config/G_defaults.json'))
data=ginputs.Inputs(str(engine),cfg,{'tensors_D':str(rebuild/'run/D/data/D/tensors'),
                                    'tensors_E':str(rebuild/'run/E/data/E/tensors')},strict=True)
join=rebuild/'run/D_labels/results/D'
labels=pd.read_csv(join/'D_line_table.tsv',sep='\t')
cov=gcore.read_covariates(str(join/'D_line_covariates.tsv'))
rows=[]
for panel in cfg.const('panels'):
    p=data.panel(panel,cov,cfg.run('B_full'))
    aligned=gsets.Aligned(p,labels)
    tensor=cfg.const('level3_tensors')[panel]
    if panel=='P1':
        tensor=tensor['7dA']
    level3=data.level3(tensor,p,order='D')[0]
    sizes={'level3':(np.nanmean(gcore.response_size(level3,~np.isnan(level3[:,:,0])),axis=1),tensor,'primary'),
           'level5':(np.nanmean(gcore.response_size(p.R,p.avail),axis=1),f'{panel}_resp','S12')}
    for marker in ('factY','factX'):
        d=aligned.design(marker,gsets.primary_settings(cfg))
        for scale,(pooled,source,role) in sizes.items():
            rows.append(ganalyses.e3_row(SimpleNamespace(cfg=cfg),panel,marker,d,pooled,scale,source,role))
fields=['panel','term','scale','tensor','role','unit','interval','n_lines','df','estimate','lo90','hi90','lo95','hi95']
key=['panel','term','scale','role']
out=args.output/'figure4_intervals.tsv'
pd.DataFrame(rows)[fields].to_csv(out,index=False,sep='\t')
new=pd.read_csv(out,sep='\t',float_precision='round_trip').set_index(key)
old=pd.read_csv(engine/'results/G/G_H3.tsv',sep='\t',float_precision='round_trip').set_index(key)
assert len(new)==12 and new.index.is_unique and old.index.is_unique
old=old.loc[new.index,new.columns]
different=(new.ne(old) & ~(old.isna() & new.isna())).sum()
result={'rows':len(new),'compared_nonkey_fields':len(new.columns),'different_cells':{k:int(v) for k,v in different.items() if v}}
(args.output/'figure4_intervals_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(result)
