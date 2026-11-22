"""Apply two archived sets of cutoffs to identical uncensored well totals and per-plate references."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

ROOT=Path.home()/"p59_data/tahoe/partA_local"
HERE=Path(__file__).resolve().parent
df=pd.read_csv(HERE/"well_counts.tsv.gz",sep="\t").merge(pd.read_csv(HERE/"sample_key.tsv",sep="\t")[["sample","dose_um"]],on="sample",validate="m:1")
df["drug"]=df.drug.str.strip()
genes=["XIST","RPS4Y1","DDX3Y","KDM5D","UTY","EIF1AY","ZFY","USP9Y","NLGN4Y"]
cols=["total_umi",*genes]
veh=df.drug.eq("DMSO_TF")
all_veg=df.loc[veh].groupby("cell_line_id")[cols].sum()
refs=pd.concat({p:df.loc[veh&df.plate.ne(p)].groupby("cell_line_id")[cols].sum() for p in sorted(df.plate.unique())},names=["plate","cell_line_id"])
s=df[df.cells.ge(200)].copy()
out={}
for title,cutfile,statefile in [
 ("lineplate","partA_local_cutoffs.tsv","partA_local_states.tsv"),
 ("line","partA_local_cutoffs_linelevel.tsv","partA_local_states_linelevel.tsv")
]:
    cut=pd.read_csv(ROOT/cutfile,sep="\t").set_index("score").cutoff
    state=pd.read_csv(ROOT/statefile,sep="\t").set_index("cellosaurus_id")
    x_frozen=state.XIST>=cut.XIST
    y_frozen=state.y_score>=cut.y_score
    # Archived printed scores of every line are farther than rounding precision from both cutoffs.
    assert abs(state.XIST-cut.XIST).min()>0.00005 and abs(state.y_score-cut.y_score).min()>0.00005
    frozen=2*x_frozen.astype(int)+y_frozen.astype(int)
    def call(m):
        v=np.log2(m[genes].div(m.total_umi,axis=0)*1e6+1)
        return (2*(v.XIST>=cut.XIST).astype(int)+(v[genes[1:]].mean(axis=1)>=cut.y_score).astype(int))
    a=s.copy()
    a["test"]=call(a)
    a=a.join(call(refs).rename("heldout"),on=["plate","cell_line_id"]).join(frozen.rename("frozen"),on="cell_line_id")
    a["agree"]=a.test==a.heldout
    stats={}
    for name,g in a.groupby(np.where(a.drug=="DMSO_TF","vehicle","treated")):
        c=g[g.cell_line_id=="CVCL_0218"]
        stats[name]=dict(n=len(g),lines=g.cell_line_id.nunique(),equal_line_agreement=float(g.groupby("cell_line_id").agree.mean().mean()),
                         pooled_agreement=float(g.agree.mean()),discordant=int((~g.agree).sum()),
                         COLO_205=dict(n=len(c),matching=int(c.agree.sum()),discordant=int((~c.agree).sum()),agreement=float(c.agree.mean()),
                                       reference_codes=sorted(c.heldout.unique().tolist())))
    gg=a[a.drug!="DMSO_TF"]
    rr=gg[gg.plate.isin(["plate6","plate14"])].pivot(index=["drug","dose_um","cell_line_id"],columns="plate",values="test").dropna()
    rr=rr.join(frozen.rename("frozen"),on="cell_line_id")
    rr["same"]=rr.plate6==rr.plate14
    rr["new"]=rr.same & rr.plate6.ne(rr.frozen)
    c=rr.xs("CVCL_0218",level="cell_line_id")
    common=a.groupby(["cell_line_id",a.drug.ne("DMSO_TF")]).agree.mean().unstack().dropna()
    check=call(all_veg).reindex(frozen.index)==frozen
    mismatch=a[a.frozen!=a.heldout]
    stats["reference"]=dict(agreement_all_vehicle_with_frozen=int(check.sum()),n_lines=len(check),mismatch_all_vehicle=state.loc[~check,"cell_name"].tolist(),
                            eligible_well_line_refs_differ_from_frozen=len(mismatch),COLO_archived_vehicle=int(frozen["CVCL_0218"]))
    stats["repeat"]=dict(pairs=len(rr),same=int(rr.same.sum()),same_new=int(rr.new.sum()),
                         same_new_by_line=rr[rr.new].reset_index().cell_line_id.map(state.cell_name).value_counts().to_dict(),
                         COLO_pairs=len(c),COLO_same=int(c.same.sum()),COLO_same_new=int(c.new.sum()),
                         COLO_both_same_other=dict(c[c.new].plate6.astype(int).value_counts().to_dict()))
    stats["common_lines_treated_minus_vehicle"]=float((common[True]-common[False]).mean())
    stats["cutoffs"]=dict(XIST=float(cut.XIST),Y=float(cut.y_score))
    stats["input_hashes"]={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in [cutfile,statefile]}
    out[title]=stats
(HERE/"independent_check_thresholds.json").write_text(json.dumps(out,indent=2)+"\n")
print(json.dumps(out["line"],indent=2))
