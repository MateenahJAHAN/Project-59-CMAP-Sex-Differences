"""Audit the marker transfer results using gene-token pseudobulks, separately from the per-cell marker tables."""
import argparse
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
Y_GENES = ["RPS4Y1","DDX3Y","KDM5D","UTY","EIF1AY","ZFY","USP9Y","NLGN4Y"]
GENES = ["XIST",*Y_GENES]
KEY = ["plate","sample","cell_line_id"]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream,"sha256").hexdigest()


def rebuild(root, baselines):
    genes = pd.read_parquet(root/"index/meta/gene_metadata.parquet")
    selected = genes[genes.gene_symbol.isin(GENES)].set_index("token_id").gene_symbol
    assert len(selected) == len(GENES) and selected.is_unique
    hashes = dict(row.split(maxsplit=1)[::-1] for row in (root/"stream_SHA256SUMS.txt").read_text().splitlines())
    plates=[]
    for plate in range(1,15):
        base=f"plates/plate{plate}/"
        path=root/"stream"/base/"pseudobulk_counts.parquet"
        assert digest(path)==hashes[base+"pseudobulk_counts.parquet"]
        grp_path=path.with_name("pseudobulk_groups.tsv")
        assert digest(grp_path)==hashes[base+"pseudobulk_groups.tsv"]
        # The sparse archive stores all gene tokens by marker-defined stratum.
        ft=pq.read_table(path,columns=["sample","cell_line_id","stratum","token_id","count"],filters=[("token_id","in",selected.index.tolist())])
        sp=ft.to_pandas()
        sp["gene"]=sp.token_id.map(selected)
        keys=["sample","cell_line_id","stratum"]
        assert not sp.duplicated(keys+["gene"]).any()
        totals = pd.read_csv(grp_path,sep="\t").set_index(keys)[["n_cells","umi"]]
        marks=sp.pivot_table(index=keys,columns="gene",values="count",aggfunc="sum",fill_value=0,observed=True)
        frame=totals.join(marks).fillna(0).reset_index()
        frame["plate"]=f"plate{plate}"
        frame=frame.merge(baselines[["xist_high","y_high"]],left_on="cell_line_id",right_index=True,validate="m:1")
        # The original state-dependent filter removed cells with >=3 discordant marker counts.
        frame["dependent_drop"] = ((~frame.y_high)&frame.stratum.isin([2,3]))|((~frame.xist_high)&frame.stratum.isin([1,3]))
        plates.append(frame)
        print(f"Recounted original token pseudobulks for plate {plate}",flush=True)
    return pd.concat(plates,ignore_index=True)


def evaluate(frame,meta,cut):
    cols=["n_cells","umi",*GENES]
    totals=frame.groupby(KEY,as_index=False)[cols].sum().merge(meta,left_on="sample",right_index=True,validate="m:1")
    assert (totals.plate==totals.plate_meta).all()
    totals["drug"]=totals.drug.str.strip()
    veh=totals.drug.eq("DMSO_TF")
    # Reference by explicit aggregation on the thirteen remaining plates.
    refs=[]
    for p in sorted(meta.plate_meta.unique()):
        ref=totals[veh & totals.plate.ne(p)].groupby("cell_line_id",as_index=False)[cols].sum()
        ref["plate"]=p
        refs.append(ref)
    def score(f):
        x=np.log2(1e6*f.XIST/f.umi + 1)
        y=np.array([np.log2(1+1e6*f[g]/f.umi) for g in Y_GENES]).mean(axis=0)
        return x,y,(np.asarray(x>=cut["XIST"],dtype=int)*2 + np.asarray(y>=cut["y_score"],dtype=int))
    refs=pd.concat(refs,ignore_index=True)
    refs["x_ref"],refs["y_ref"],refs["ref_call"]=score(refs)
    totals["x_score"],totals["y_score"],totals["joint"]=score(totals)
    out=totals.merge(refs[["plate","cell_line_id","x_ref","y_ref","ref_call"]],on=["plate","cell_line_id"],validate="m:1")
    out["match"]=out.joint.eq(out.ref_call)
    return out


def summarize(df):
    rate=[]
    for minimum in [50,200,500]:
        for label,part in df[df.n_cells>=minimum].groupby(df.drug.ne("DMSO_TF")):
            per=part.groupby("cell_line_id").match.mean()
            rate.append(dict(min_cells=minimum,treated=bool(label),n_lines=len(per),n_line_wells=len(part),
                             discordant=int((~part.match).sum()),equal_line=float(per.mean())))
    return rate


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=Path.home()/"p59_data/tahoe")
    root=p.parse_args().root
    state=pd.read_csv(root/"partA_local/partA_local_states.tsv",sep="\t").set_index("cellosaurus_id")
    cut=pd.read_csv(root/"partA_local/partA_local_cutoffs.tsv",sep="\t").set_index("score").cutoff
    orig_cut=pd.read_csv(HERE/"reference_cutoffs.tsv",sep="\t").set_index("score").cutoff
    assert cut.equals(orig_cut)
    frozen=(2*state.xist_high.astype(int)+state.y_high.astype(int)).rename("frozen")
    meta=pd.read_parquet(root/"index/meta/sample_metadata.parquet").set_index("sample").rename(columns={"plate":"plate_meta"})
    meta["dose"]=meta.drugname_drugconc.map(lambda v:ast.literal_eval(v)[0][1])
    meta=meta[["drug","dose","plate_meta"]]
    counts=rebuild(root,state)
    df=evaluate(counts,meta,cut)
    published=pd.read_csv(HERE/"well_counts.tsv.gz",sep="\t").set_index(KEY).sort_index()
    equivalent=df.set_index(KEY).sort_index().rename(columns={"n_cells":"cells","umi":"total_umi"})
    assert published.index.equals(equivalent.index)
    mismatches={}
    for c in ["cells","total_umi",*GENES]:
        n=int((equivalent[c]!=published[c]).sum())
        assert n==0
        mismatches[c]=n
    previously=pd.read_csv(HERE/"line_well_calls.tsv.gz",sep="\t").set_index(KEY).sort_index()
    assert equivalent.index.equals(previously.index)
    assert np.array_equal(equivalent.joint,previously.call) and np.array_equal(equivalent.ref_call,previously.reference_call)
    rates=summarize(df)
    seen=pd.read_csv(HERE/"coverage_and_agreement.tsv",sep="\t")
    for row in rates:
        orig=seen[(seen.minimum_cells==row["min_cells"])&(seen.treated==row["treated"])].iloc[0]
        assert (row["n_lines"],row["n_line_wells"])==(int(orig.lines),int(orig.line_wells))
        assert abs(row["equal_line"]-orig.joint_agreement)<1e-12
    df["admitted"]=df.n_cells.ge(200)
    tx=df[df.admitted & df.drug.ne("DMSO_TF")].copy()
    ve=df[df.admitted & df.drug.eq("DMSO_TF")].copy()
    vt=set(ve.cell_line_id); tt=set(tx.cell_line_id)
    outside=state.index.difference(tt)
    tx["y_flip"]=tx.joint.mod(2)!=tx.ref_call.mod(2)
    tx["x_flip"]=tx.joint.floordiv(2)!=tx.ref_call.floordiv(2)
    flips=tx[~tx.match].groupby("cell_line_id").size().sort_values(ascending=False)
    paired=tx[tx.plate.isin(["plate6","plate14"])].copy()
    pairs=paired.merge(paired,on=["drug","dose","cell_line_id"],suffixes=("_6","_14"))
    pairs=pairs[(pairs.plate_6=="plate6")&(pairs.plate_14=="plate14")].copy()
    assert not pairs.duplicated(["drug","dose","cell_line_id"]).any()
    pairs["frozen"]=pairs.cell_line_id.map(frozen)
    pairs["same_new"]=pairs.joint_6.eq(pairs.joint_14)&pairs.joint_6.ne(pairs.frozen)
    expected=pd.read_csv(HERE/"replicate_pairs.tsv.gz",sep="\t")
    j=pairs.merge(expected,left_on=["drug","dose","cell_line_id"],right_on=["drug","dose_um","cell_line_id"],validate="1:1")
    assert len(pairs)==len(j)==len(expected) and (j.same_new==j.both_same_new).all()
    both=pairs[pairs.same_new]
    # This sensitivity applies the earlier state-dependent exclusion to both test and reference, retrospectively.
    old=evaluate(counts[~counts.dependent_drop],meta,cut)
    oldrate=summarize(old)
    common=tx.merge(old[old.n_cells.ge(200)],on=KEY,suffixes=("_raw","_filtered"),validate="1:1")
    exact_holdout_change=tx[tx.ref_call!=tx.cell_line_id.map(frozen)].groupby("cell_line_id").size().to_dict()
    paired_rates=tx.groupby("cell_line_id").match.mean().rename("tx").to_frame().join(ve.groupby("cell_line_id").match.mean().rename("vehicle"),how="inner")
    colo=state.index[state.cell_name=="COLO 205"][0]
    tc=tx[tx.cell_line_id==colo]
    other=paired_rates.drop(index=colo)
    repeat_meta=meta[meta.plate_meta.isin(["plate6","plate14"]) & meta.drug.ne("DMSO_TF")].copy()
    repeat_meta["drug"]=repeat_meta.drug.str.strip()
    samplepairs=repeat_meta.groupby(["drug","dose"]).plate_meta.nunique()
    missing=set(samplepairs[samplepairs==2].index)-set(pairs[["drug","dose"]].drop_duplicates().itertuples(index=False,name=None))
    result=dict(compared_present_line_wells=len(equivalent),all_possible_line_wells=1344*50,
                verified_cells=int(equivalent.cells.sum()),mismatches=mismatches,
                checks_identical="per-well calls, leave-plate-out references, coverage and repeat pair keys",
                full_filter_uncensored=rates,thresholds={k:float(cut[k]) for k in ["XIST","y_score"]},
                minimum_200_line_groups=dict(eligible_both=len(vt&tt),treated_only=state.loc[list(tt-vt),"cell_name"].tolist(),vehicle_only=state.loc[list(vt-tt),"cell_name"].tolist(),no_eligible_treated=state.loc[outside,"cell_name"].tolist()),
                primary_disc_by_line={state.loc[i,"cell_name"]:int(n) for i,n in flips.items()},
                marker_failures=dict(y=int(tx.y_flip.sum()),xist=int(tx.x_flip.sum()),both=int((tx.y_flip&tx.x_flip).sum())),
                heldout_reference_differs_from_frozen_by_line={state.loc[i,"cell_name"]:int(n) for i,n in exact_holdout_change.items()},
                colo=dict(eligible_treated=len(tc),eligible_vehicle=len(ve[ve.cell_line_id==colo]),
                          raw_score_median=float(tc.y_score.median()),raw_score_quartiles=tc.y_score.quantile([.25,.75]).tolist(),
                          reference_y_across_plates=tc.groupby("plate").y_ref.first().to_dict(),
                          group_minus_COLO_205_equal_line_treated_minus_vehicle=float((other.tx-other.vehicle).mean())),
                repeat=dict(source_pairs=len(samplepairs[samplepairs==2]),missing_from_eligible=sorted(missing),eligible=len(pairs),
                            dual_matching=int((pairs.joint_6==pairs.joint_14).sum()),both_same_new=len(both),
                            names=both.cell_line_id.map(state.cell_name).value_counts().to_dict(),
                            ref_to_new=[{"from":int(a),"to":int(b),"n":int(k)} for (a,b),k in both.groupby(["frozen","joint_6"]).size().items()],
                            new_y_score_p6_range=[float(both.y_score_6.min()),float(both.y_score_6.max())],
                            new_y_score_p14_range=[float(both.y_score_14.min()),float(both.y_score_14.max())]),
                retrospective_conditional_filter=oldrate,
                identical_wells_filter_sensitivity=dict(eligible_treated_under_both=len(common),
                    uncensored_misses=int((~common.match_raw).sum()),conditional_filtered_misses=int((~common.match_filtered).sum())),
                cells_removed_by_retrospective_filter=int(counts.loc[counts.dependent_drop,"n_cells"].sum()))
    (HERE/"independent_check.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
