"""Recreate each 10,000-draw distribution in Figure 4 from the rebuilt P1 panel and covariates."""
import argparse
import gzip
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

parser=argparse.ArgumentParser(description=__doc__)
for name in ('engine','rebuild','archive-capture','output'):
    parser.add_argument('--'+name,required=True,type=Path)
args=parser.parse_args()
engine=args.engine.resolve()
rebuild=args.rebuild.resolve()
args.output.mkdir(exist_ok=True,parents=True)
sys.path.insert(0,str(engine/'src/G'))
import gcore
import gconfig
import ginputs
import gsets
import ganalyses

cfg=gconfig.load(str(engine/'config/G_defaults.json'))
inputs=ginputs.Inputs(str(engine),cfg,{'tensors_D':str(rebuild/'run/D/data/D/tensors')},strict=True)
join=rebuild/'run/D_labels/results/D'
p=inputs.panel('P1',gcore.read_covariates(str(join/'D_line_covariates.tsv')),cfg.run('B_full'))
a=gsets.Aligned(p,pd.read_csv(join/'D_line_table.tsv',sep='\t'))
primary=gsets.primary_settings(cfg)
manifest=pd.read_csv(args.archive_capture/'capture_manifest.tsv',sep='\t')
B=cfg.run('B_full')
assert B==10000
outcomes={}
draws={}
for term in ('factY','factX','TP53'):
    if term=='TP53':
        d=a.design('factY',primary)
        X=gsets.tp53_columns(p.cov,'7gB' if gsets.tp53_code(cfg)=='7gB' else 'cmp')[0]
        k=0
        mask=gcore.analysis_set(d['lin'],X,[k],d['mask'])
        extra=None
        context='exposure|TP53 status (comparison exposure)|P1'
    else:
        d=a.design(term,primary)
        X,k,mask=d['Xo'],d['k'],d['mask']
        extra={ganalyses.DMSO:p.dmso}
        context='main|primary|P1|'+term
    df,perms,_=gcore.run_term(p,X,k,mask,B,extra=extra,lin=d['lin'],rule=d['rule'])
    df,observed,pvalue,null=gcore.finish_table(df,perms,compounds_only=p.comps)
    archived=manifest[(manifest.context==context)&(manifest.kind=='finish_table')].iloc[0]
    with np.load(args.archive_capture/archived.file) as frozen:
        changes=int(np.count_nonzero(null!=frozen['null']))
        observed_equal=bool(frozen['panel_stat']==observed)
        p_equal=bool(frozen['panel_p']==pvalue)
    hits=int(np.count_nonzero(null>=observed*(1-gcore.TIE_TOL)))
    direct=(1+hits)/(B+1)
    assert direct==pvalue
    outcomes[term]=dict(lines=int(mask.sum()),tested=int(df.in_ranking.sum()),draws=B,
                        at_or_above=hits,observed=float(observed),panel_p=float(pvalue),
                        unequal_draws=changes,observed_matches_archive=observed_equal,p_matches_archive=p_equal)
    draws[term]=null
    print(term,outcomes[term],flush=True)
with (args.output/'figure4_nulls.tsv.gz').open('wb') as raw:
    with gzip.GzipFile(fileobj=raw,filename='',mtime=0,mode='wb') as compressed:
        frame=pd.DataFrame(draws)
        frame.index.name='draw'
        compressed.write(frame.to_csv(sep='\t',index=True).encode())
(args.output/'figure4_nulls_check.json').write_text(json.dumps(outcomes,indent=2)+'\n')
