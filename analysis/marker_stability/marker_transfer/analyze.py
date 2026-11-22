"""Estimate agreement of two RNA marker calls across all eligible Tahoe line-wells."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from aggregate import MARKERS

OUT = Path(__file__).resolve().parent
KEY = "cell_line_id"
measure = ["cells","total_umi"] + MARKERS
thresholds = pd.read_csv(OUT/"reference_cutoffs.tsv",sep="\t").set_index("score").cutoff
ref = pd.read_csv(OUT/"reference_lines.tsv",sep="\t").rename(columns={"cellosaurus_id":KEY}).set_index(KEY)
ref["frozen"] = ref.xist_high.astype(int)*2 + ref.y_high.astype(int)
meta = pd.read_csv(OUT/"sample_key.tsv",sep="\t")
raw = pd.read_csv(OUT/"well_counts.tsv.gz",sep="\t")
# Seven atlas metadata names contain a trailing space absent from the cell table.
meta["drug"] = meta.drug.str.strip()
raw["drug"] = raw.drug.str.strip()
raw = raw.merge(meta,on=["sample","plate","drug"],validate="many_to_one")
assert raw["sample"].nunique() == len(meta) == 1344 and raw[KEY].nunique() == 50
assert raw.cells.sum() == 95624334 and (raw[MARKERS].sum(axis=1) <= raw.total_umi).all()

def call(frame):
    expr = np.log2(1 + frame[MARKERS].to_numpy(dtype=float) * (1.e6/frame["total_umi"].to_numpy())[:,None])
    return 2*(expr[:,0]>=thresholds.XIST).astype(int) + (expr[:,1:].mean(axis=1)>=thresholds.y_score)

vehicle = raw[raw.drug.eq("DMSO_TF")]
all_veh = vehicle.groupby(KEY)[measure].sum()
same_veh = vehicle.groupby(["plate",KEY])[measure].sum()
assert len(same_veh)==700
outside = all_veh.reindex(same_veh.index.get_level_values(KEY)).set_axis(same_veh.index) - same_veh
assert (outside.cells >= 200).all()
outside["reference_call"] = call(outside)
same_veh["same_plate_call"] = call(same_veh)
raw["call"] = call(raw)
raw["treated"] = raw.drug.ne("DMSO_TF")
raw = raw.merge(outside[["reference_call"]],on=["plate",KEY],validate="many_to_one")
raw["joint_agreement"] = raw.call.eq(raw.reference_call)
raw["xist_agreement"] = (raw.call//2).eq(raw.reference_call//2)
raw["y_agreement"] = (raw.call%2).eq(raw.reference_call%2)
rates = ["joint_agreement","xist_agreement","y_agreement"]
primary = raw[raw.cells.ge(200)].copy()
def eligible_summary(cutoff):
    kept = raw[raw.cells.ge(cutoff)]
    rows = []
    for treated,group in kept.groupby("treated"):
        by_line = group.groupby(KEY)[rates].mean()
        rows.append(dict(minimum_cells=cutoff,treated=bool(treated),lines=len(by_line),line_wells=len(group),
                         possible_line_wells=(1316 if treated else 28)*50,pooled=float(group.joint_agreement.mean()),
                         **{col:float(by_line[col].mean()) for col in rates}))
    return rows
detail = [row for n in [50,200,500] for row in eligible_summary(n)]
pd.DataFrame(detail).to_csv(OUT/"coverage_and_agreement.tsv",sep="\t",index=False)
def table(keys):
    return primary.groupby(keys,as_index=False).agg(eligible_line_wells=("sample","size"), **{k:(k,"mean") for k in rates})
line_table = table(["treated",KEY]).merge(ref[["cell_name","organ","frozen"]],on=KEY)
line_table.to_csv(OUT/"by_line.tsv",sep="\t",index=False)
table(["treated","plate"]).to_csv(OUT/"by_plate.tsv",sep="\t",index=False)
pairs = line_table.pivot(index=KEY,columns="treated",values="joint_agreement").dropna()
delta = (pairs[True]-pairs[False]).to_numpy()
rng = np.random.default_rng(20261121)
boot = rng.choice(delta,(10000,len(delta)),replace=True).mean(axis=1)
agree = (pd.Series(call(all_veh), index=all_veh.index)==ref.loc[all_veh.index,"frozen"])
summary = dict(conditions=detail, archived_vehicle_pattern_agree=int(agree.sum()),total_lines=len(agree),
               equal_line_comparison=dict(lines=len(delta),vehicle=float(pairs[False].mean()),treated=float(pairs[True].mean()),
                   difference=float(delta.mean()),interval=np.quantile(boot,[0.025,0.975]).tolist()))
# Use the exact same line-wells for both same-plate summaries and the thinning.
tx = primary[primary.treated].copy()
same = tx.merge(same_veh.add_suffix("_v"),on=["plate",KEY],validate="many_to_one")
same = same[same.cells_v.ge(200)].copy()
def molecule_counts(frame,ending):
    specific = frame[[g+ending for g in MARKERS]].to_numpy(dtype=np.int64)
    total = frame["total_umi"+ending].to_numpy(dtype=np.int64)
    return np.column_stack([specific, total-specific.sum(axis=1)])
left,right = molecule_counts(same,""),molecule_counts(same,"_v")
rng = np.random.default_rng(20261121)
def thin(a,b):
    for i,(x,y) in enumerate(zip(a,b)):
        n,m = int(x.sum()),int(y.sum())
        if n>m: a[i] = rng.multivariate_hypergeometric(x,m,method="marginals")
        if m>n: b[i] = rng.multivariate_hypergeometric(y,n,method="marginals")
thin(left,right)
def call_array(a):
    counts = pd.DataFrame(a[:,:-1],columns=MARKERS)
    counts["total_umi"] = a.sum(axis=1)
    return call(counts)
same["matched"] = call_array(left)==call_array(right)
same["unmatched"] = same.call.eq(same.same_plate_call_v)
summary["same_plate"] = {"line_wells":len(same),**same.groupby(KEY)[["joint_agreement","unmatched","matched"]].mean().mean().to_dict()}
# Across biological replicate plates, match exact drug and concentration within each line.
repeat = tx[tx.plate.isin(["plate6","plate14"])].copy()
idx=[KEY,"drug","dose_um"]
assert not repeat.duplicated(idx+["plate"]).any()
repeat = repeat.pivot(index=idx,columns="plate",values="call").dropna().astype(int).reset_index()
repeat = repeat.join(ref[["frozen","cell_name"]],on=KEY)
repeat["agree"] = repeat.plate6.eq(repeat.plate14)
repeat["both_same_new"] = repeat.agree & repeat.plate6.ne(repeat.frozen)
repeat.to_csv(OUT/"replicate_pairs.tsv.gz",sep="\t",index=False,compression={"method":"gzip","mtime":0})
byline = repeat.groupby([KEY,"cell_name"]).agg(pairs=("agree","size"),agree=("agree","sum"),both_same_new=("both_same_new","sum"))
byline.to_csv(OUT/"replicates_by_line.tsv",sep="\t")
summary["replicate_plates"] = dict(pairs=len(repeat),drugdoses=len(repeat[["drug","dose_um"]].drop_duplicates()),lines=len(byline),
    matching=int(repeat.agree.sum()),same_new=int(repeat.both_same_new.sum()),lines_same_new=int(byline.both_same_new.gt(0).sum()))
raw[["plate","sample",KEY,"drug","dose_um","cells","total_umi","call","reference_call","joint_agreement"]].to_csv(OUT/"line_well_calls.tsv.gz",
    sep="\t",index=False,compression={"method":"gzip","mtime":0})
(OUT/"results.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
