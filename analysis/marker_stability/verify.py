"""Reconstruct the entire analysis matrix from individual cell annotations; no imports from analyse.py."""
from pathlib import Path
import argparse
import ast
import hashlib
import itertools
import json
import math
import sys

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]


def public_inputs(work):
    sys.path.insert(0, str(REPO))
    from analysis.part_a.run import acquire, REPLAY
    receipt = json.loads((REPLAY / "part_a/receipts/tahoe_observations.json").read_text())
    path = work / "tahoe/index/meta/obs_metadata.parquet"
    acquire(dict(url=receipt["url"], bytes=receipt["file_bytes"], sha256=receipt["sha256"]), path)
    return path, receipt


def counts_by_well(observations, ids, samples, raw):
    matrix = np.zeros((len(samples), len(ids)), dtype=np.int64)
    def sample_num(s):
        return int(s.removeprefix("smp_"))
    slot = np.full(max(map(sample_num, samples)) + 1, -1)
    for i, name in enumerate(samples):
        slot[sample_num(name)] = i
    source = pq.ParquetFile(observations)
    scanned = 0
    filtered = "pass_filter" if raw else "pass_full"
    for batch in source.iter_batches(batch_size=262144, columns=[filtered, "cell_line", "sample"]):
        scanned += batch.num_rows
        flag, cell, sample = batch.columns
        if any(column.null_count for column in batch.columns):
            raise ValueError("Missing filter, cell line or well in the Tahoe annotations")
        if raw:
            flag = pc.equal(flag, "full")
            sample = pc.cast(pc.utf8_slice_codeunits(sample, 4), pa.int32())
        use = flag.to_numpy(zero_copy_only=False)
        to_axis = np.array([ids.index(k) for k in cell.dictionary.to_pylist()])
        line = to_axis[cell.indices.to_numpy()[use]]
        codes = sample.to_numpy()[use]
        if np.any(codes < 0) or np.any(codes >= len(slot)) or np.any(slot[codes] < 0):
            raise ValueError("Full-filter observation with unknown sample ID")
        matrix += np.bincount(slot[codes] * len(ids) + line, minlength=matrix.size).reshape(matrix.shape)
    if scanned != source.metadata.num_rows:
        raise ValueError(f"Read {scanned} observation rows, expected {source.metadata.num_rows}")
    return matrix, scanned


def compare_plates(counts, meta, root):
    with (root / "stream_SHA256SUMS.txt").open() as f:
        hashes = dict(row.rstrip().split(maxsplit=1)[::-1] for row in f)
    mismatches = 0
    for plate in range(1, 15):
        rel = f"plates/plate{plate}/pseudobulk_groups.tsv"
        path = root / "stream" / rel
        with path.open("rb") as f:
            if hashlib.file_digest(f, "sha256").hexdigest() != hashes[rel]:
                raise ValueError(f"Unexpected digest for {rel}")
        expected = pd.read_csv(path, sep="\t")
        plate_id = f"plate{plate}"
        wells = meta.index[meta.plate == plate_id]
        if (len(wells) != 96 or set(expected.plate) != {plate_id} or set(expected["sample"]) != set(wells)
                or not set(expected.cell_line_id).issubset(counts.columns)):
            raise ValueError(f"Sample or line mismatch in {rel}")
        released = expected.groupby(["sample", "cell_line_id"]).n_cells.sum().unstack(fill_value=0)
        released = released.reindex(index=wells, columns=counts.columns, fill_value=0)
        mismatches += int(np.count_nonzero(released.to_numpy() != counts.loc[wells].to_numpy()))
    if mismatches:
        raise ValueError(f"{mismatches} reconstructed well–line counts differ from the fourteen plate tables")
    return mismatches


def main():
    parser = argparse.ArgumentParser()
    sources = parser.add_mutually_exclusive_group()
    sources.add_argument("--work", type=Path, help="Public work directory for the pinned Tahoe provider file")
    sources.add_argument("--root", type=Path, help="Historical stream directory with an existing compact observation table")
    parser.add_argument("--out", type=Path, help="Results directory, also read for the primary effect summary if available")
    args = parser.parse_args()
    receipt = None
    if args.work:
        work = args.work.expanduser().resolve()
        if work.is_relative_to(REPO):
            parser.error("Choose a work directory outside the source checkout")
        observations, receipt = public_inputs(work)
        root = HERE / "inputs"
        out = args.out or work / "marker_stability/results"
    else:
        root = (args.root or Path.home() / "p59_data/tahoe").expanduser()
        observations = root / "index/obs_compact.parquet"
        out = args.out or HERE / "results"
    out.mkdir(parents=True, exist_ok=True)
    meta = pd.read_parquet(root/"index/meta/sample_metadata.parquet").set_index("sample")
    states = pd.read_csv(root/"partA_local/partA_local_states.tsv",sep="\t").set_index("cellosaurus_id")
    if args.work:
        current_states = work / "results/tahoe/latest/partA_local_states.tsv"
        if current_states.is_file():
            current = pd.read_csv(current_states, sep="\t").set_index("cellosaurus_id")
            keys = ["organ", "y_high", "xist_high"]
            pd.testing.assert_frame_equal(states[keys].sort_index(), current[keys].sort_index())
            states = current
    if not states.index.is_unique or len(states) != 50 or not meta.index.is_unique or len(meta) != 1344:
        raise ValueError("Expected 50 unique cell lines and 1,344 unique wells")
    ids = list(sorted(states.index))
    samples = list(sorted(meta.index))
    mat, observed = counts_by_well(observations, ids, samples, raw=receipt is not None)
    if receipt and (observed != receipt["rows"] or int(mat.sum()) != receipt["full_filter_rows"]):
        raise ValueError("Observed row counts differ from the provider receipt")
    df = pd.DataFrame(mat,index=samples,columns=ids)
    n_bad = compare_plates(df, meta, root)
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
    summary = out / "summary.json"
    main_summary=json.loads((summary if summary.is_file() else HERE/"results/summary.json").read_text())
    assert np.isclose(estimate,main_summary["primary"]["estimate"],atol=1e-11)
    assert int(mat.sum())==main_summary["counts"]["cells"]
    exact=float(np.mean(np.abs(null)+1e-12>=abs(estimate)))
    assert abs(exact-main_summary["primary"]["permutation_p"])<.01
    info=dict(raw_annotation_rows=observed,captured_rows=int(mat.sum()),
              well_line_entries=int(mat.size),mismatched_count_entries=n_bad,estimate=estimate,
              exact_two_sided_p=exact,unique_labelings=len(null),main_monte_carlo_p=main_summary["primary"]["permutation_p"])
    if receipt:
        df.to_csv(out/"independent_well_line_counts.tsv.gz", sep="\t", index_label="sample",
                  compression={"method":"gzip","mtime":0})
    (out/"independent_verification.json").write_text(json.dumps(info,indent=2)+"\n")
    print(json.dumps(info,indent=2))


if __name__=="__main__":
    main()
