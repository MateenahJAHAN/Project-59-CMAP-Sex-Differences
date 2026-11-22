"""Verify whole-line bootstrap and same-plate controls using independently validated well totals."""
from pathlib import Path
import json

import numpy as np
import pandas as pd

base=Path(__file__).resolve().parent
raw=pd.read_csv(base/"well_counts.tsv.gz",sep="\t")
rep=json.loads((base/"results.json").read_text())
keys=["plate","cell_line_id"]
y=["RPS4Y1","DDX3Y","KDM5D","UTY","EIF1AY","ZFY","USP9Y","NLGN4Y"]
genes=["XIST",*y]
thr=pd.read_csv(base/"reference_cutoffs.tsv",sep="\t").set_index("score").cutoff
v=raw[raw.drug.eq("DMSO_TF")]
t=raw[~raw.drug.eq("DMSO_TF") & raw.cells.ge(200)].copy()
s=t.merge(v.groupby(keys)[genes+["total_umi","cells"]].sum(),on=keys,suffixes=("","_c"),validate="m:1")
s=s[s.cells_c.ge(200)]
def matrix(suffix):
    m=s[[g+suffix for g in genes]].to_numpy(dtype=np.int64)
    return np.column_stack([m,s["total_umi"+suffix].to_numpy()-m.sum(axis=1)])
def label(a):
    obs=np.log2(1+1e6*a[:,:9]/a.sum(axis=1)[:,None])
    return (obs[:,0]>=thr.XIST)*2+(obs[:,1:].mean(axis=1)>=thr.y_score)
def equal_line_rate(a,b):
    return float(pd.Series(label(a)==label(b),index=s.cell_line_id).groupby(level=0).mean().mean())
m1,m2=matrix(""),matrix("_c")
unmatched=equal_line_rate(m1,m2)
assert len(s)==rep["same_plate"]["line_wells"] and np.isclose(unmatched,rep["same_plate"]["unmatched"])
thins=[]
for seed in (20261121,591122,591123):
    rng=np.random.default_rng(seed)
    a,b=m1.copy(),m2.copy()
    for i in range(len(s)):
        n=a[i].sum();m=b[i].sum()
        if n>m: a[i]=rng.multivariate_hypergeometric(a[i],m,method="marginals")
        if n<m: b[i]=rng.multivariate_hypergeometric(b[i],n,method="marginals")
    assert np.array_equal(a.sum(axis=1),b.sum(axis=1))
    thins.append(equal_line_rate(a,b))
assert np.isclose(thins[0],rep["same_plate"]["matched"])
rawcalls=pd.read_csv(base/"line_well_calls.tsv.gz",sep="\t")
q=rawcalls[rawcalls.cells.ge(200)].copy()
q["treated"]=q.drug.ne("DMSO_TF")
g=q.groupby(["cell_line_id","treated"]).joint_agreement.mean().unstack().dropna()
diff=(g[True]-g[False]).to_numpy()
rng=np.random.default_rng(20261121)
boot=diff[rng.integers(len(diff),size=(10000,len(diff)))].mean(axis=1)
limits=np.quantile(boot,[.025,.975]).tolist()
assert np.allclose(limits,rep["equal_line_comparison"]["interval"])
out=dict(common_lines=len(g),heldout_equal_line_delta=float(diff.mean()),bootstrap_limits=limits,
         same_plate_wells=len(s),same_plate_lines=int(s.cell_line_id.nunique()),raw_same_plate_rate=unmatched,
         thinned_same_plate_rates=dict(zip([20261121,591122,591123],thins)))
(base/"independent_check_controls.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out,indent=2))
