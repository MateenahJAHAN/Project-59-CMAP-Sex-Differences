"""Reconstruct the entire analysis matrix from individual cell annotations; no imports from analyse.py."""
from pathlib import Path
import argparse
import ast
import itertools
import json
import math

import numpy as np
import pandas as pd
import pyarrow.parquet as pq


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.home()/"p59_data/tahoe")
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent/"results")
    args = parser.parse_args()
    root, out = args.root, args.out
    meta = pd.read_parquet(root/"index/meta/sample_metadata.parquet").set_index("sample")
    states = pd.read_csv(root/"partA_local/partA_local_states.tsv",sep="\t").set_index("cellosaurus_id")
    ids = list(sorted(states.index))
    samples = list(sorted(meta.index))
    slot = np.full(max(map(lambda s:int(s[4:]),samples))+1,-1)
    for i,s in enumerate(samples):
        slot[int(s[4:])] = i
    mat = np.zeros((len(samples),len(ids)), dtype=np.int64)
    f = pq.ParquetFile(root/"index/obs_compact.parquet")
    for i in range(f.num_row_groups):
        rows = f.read_row_group(i, columns=["pass_full","cell_line","sample"])
        use = rows["pass_full"].to_numpy().astype(bool)
        idx = rows["cell_line"].combine_chunks()
        to_axis = np.array([ids.index(k) for k in idx.dictionary.to_pylist()])
        line = to_axis[idx.indices.to_numpy()[use]]
        sample = slot[rows["sample"].to_numpy()[use]]
        assert min(sample)>=0
        hist = np.bincount(sample*len(ids)+line,minlength=mat.size)
        mat += hist.reshape(mat.shape)
    df = pd.DataFrame(mat,index=samples,columns=ids)
    n_bad = 0
    for i in range(1,15):
        original = pd.read_csv(root/f"stream/plates/plate{i}/pseudobulk_groups.tsv",sep="\t")
        from_stream=original.groupby(["sample","cell_line_id"]).n_cells.sum().unstack(fill_value=0)
        n_bad+=int(np.count_nonzero(from_stream.to_numpy()!=df.loc[from_stream.index,from_stream.columns].to_numpy()))
    assert n_bad==0
    meta=meta.loc[samples]
    meta["dose"]=meta.drugname_drugconc.map(lambda x:ast.literal_eval(x)[0][1])
    tests = meta.drug.ne("DMSO_TF")
    passed = np.min(mat[~tests],axis=0)>=20
    sg = (states.y_high | states.xist_high).loc[ids].to_numpy()
    organ=states.organ.loc[ids].to_numpy()
    # The scalar well denominator cancels between marker groups. Omit it deliberately.
    logged=np.log2(mat+.5)
    for plate,group in meta.groupby("plate"):
        inds=meta.index.get_indexer(group.index)
        controls=meta.index.get_indexer(group.index[group.drug=="DMSO_TF"])
        logged[inds]-=logged[controls].mean(axis=0)
    drug_delta=[]
    for drug,m in meta[tests].groupby("drug"):
        hi=meta.index.get_indexer(m.index[m.dose==5.])
        lo=meta.index.get_indexer(m.index[m.dose==.05])
        drug_delta.append(logged[hi].mean(axis=0)-logged[lo].mean(axis=0))
    y=np.mean(drug_delta,axis=0)
    diff=[]
    distributions=[]
    for o in sorted(set(organ)):
        use=(organ==o)&passed
        n_pos=int((use&sg).sum()); n_neg=int((use&~sg).sum())
        if min(n_pos,n_neg)<2:
            continue
        yp=y[use&sg]; yn=y[use&~sg]
        diff.append(yp.mean()-yn.mean())
        values=np.r_[yp,yn]
        sums=np.array([values[list(idx)].sum() for idx in itertools.combinations(range(len(values)),n_pos)])
        distributions.append(sums/n_pos-(values.sum()-sums)/n_neg)
    null=np.zeros(1)
    for v in distributions:
        null=(null[:,None]+v[None,:]).reshape(-1)
    null/=len(distributions)
    estimate=float(np.mean(diff))
    main_summary=json.loads((out/"summary.json").read_text())
    assert np.isclose(estimate,main_summary["primary"]["estimate"],atol=1e-11)
    assert int(mat.sum())==main_summary["counts"]["cells"]
    exact=float(np.mean(np.abs(null)+1e-12>=abs(estimate)))
    assert abs(exact-main_summary["primary"]["permutation_p"])<.01
    info=dict(raw_annotation_rows=int(f.metadata.num_rows),captured_rows=int(mat.sum()),
              well_line_entries=int(mat.size),mismatched_count_entries=n_bad,estimate=estimate,
              exact_two_sided_p=exact,unique_labelings=len(null),main_monte_carlo_p=main_summary["primary"]["permutation_p"])
    (out/"independent_verification.json").write_text(json.dumps(info,indent=2)+"\n")
    print(json.dumps(info,indent=2))


if __name__=="__main__":
    main()
