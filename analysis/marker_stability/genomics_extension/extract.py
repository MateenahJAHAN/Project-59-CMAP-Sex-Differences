"""Extract preregistered feature counts from all pseudobulk strata of both replicate plates."""

import argparse
import ast
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
CONTROL = ["FKBP5", "DDIT4", "SGK1", "PER1"]


def sha(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--tahoe", type=Path, default=Path.home() / "p59_data/tahoe")
    a.add_argument("--geneinfo", type=Path, required=True, help="Broad LINCS 2020 geneinfo_beta.txt")
    a.add_argument("--chrom", type=Path, required=True, help="Tahoe token-to-chromosome table from the fellowship caller")
    args = a.parse_args()
    src = args.tahoe
    chrom = args.chrom
    names = pd.read_csv(chrom, sep="\t", dtype={"chrom": str})
    info = pd.read_csv(args.geneinfo, sep="\t")
    lm = info.loc[info.feature_space.eq("landmark"), ["ensembl_id", "gene_symbol"]]
    assert len(lm) == 978
    auto = names[names.chrom.isin(map(str, range(1, 23)))].copy()
    assert not auto.ensembl_id.duplicated().any()
    lm = lm.merge(
        names[["ensembl_id", "token_id", "chrom", "gene_symbol"]],
        how="left",
        on="ensembl_id",
        suffixes=("_broad", "_tahoe"),
    )
    lm["retained"] = lm.chrom.isin(map(str, range(1, 23))) & ~lm.ensembl_id.duplicated(
        keep=False
    )
    lm.to_csv(HERE / "broad_landmark_mapping.tsv", sep="\t", index=False)
    selected = set(lm[lm.retained].ensembl_id)
    genes = (
        auto[auto.ensembl_id.isin(selected) | auto.gene_symbol.isin(CONTROL)]
        .copy()
        .sort_values("token_id")
        .reset_index(drop=True)
    )
    genes["landmark"] = genes.ensembl_id.isin(selected)
    assert genes.token_id.is_unique and all(
        genes.gene_symbol.eq(g).sum() == 1 for g in CONTROL
    )
    genes[["gene_symbol", "ensembl_id", "token_id", "chrom", "landmark"]].to_csv(
        HERE / "features.tsv", sep="\t", index=False
    )
    samples = src / "index/meta/sample_metadata.parquet"
    s = pd.read_parquet(samples)
    s = (
        s[s.plate.isin(["plate6", "plate14"])]
        .sort_values(["plate", "sample"])
        .reset_index(drop=True)
    )
    s["drug"] = s.drug.str.strip()
    s["micromolar"] = s.drugname_drugconc.map(lambda v: ast.literal_eval(v)[0][1])
    s["condition"] = s.drug + " | " + s.micromolar.astype(str)
    s.to_csv(HERE / "wells.tsv", sep="\t", index=False)
    sid = {v: i for i, v in enumerate(s["sample"])}
    states_path = src / "partA_local/partA_local_states.tsv"
    states = (
        pd.read_csv(states_path, sep="\t")
        .sort_values("cellosaurus_id")
        .reset_index(drop=True)
    )
    alts_path = src / "partA_local/partA_local_states_linelevel.tsv"
    th_path = src / "partA_local/partA_local_cutoffs_linelevel.tsv"
    alt = pd.read_csv(alts_path, sep="\t").set_index("cellosaurus_id")
    th = pd.read_csv(th_path, sep="\t").set_index("score").cutoff
    states["marker"] = states.y_high | states.xist_high
    states["alternate"] = states.cellosaurus_id.map(
        (alt.y_score >= th.y_score) | (alt.XIST >= th.XIST)
    )
    prior = HERE.parent / "results/lines.tsv"
    primary = pd.read_csv(prior, sep="\t").set_index("cellosaurus_id")
    states["primary"] = states.cellosaurus_id.map(primary.in_primary)
    states[
        [
            "cellosaurus_id",
            "cell_name",
            "organ",
            "annotation",
            "marker",
            "alternate",
            "primary",
        ]
    ].to_csv(HERE / "lines.tsv", sep="\t", index=False)
    lineids = states.cellosaurus_id.tolist()
    out = np.zeros((len(s), len(states), len(genes)), np.int32)
    nc = np.zeros(out.shape[:2], np.int64)
    umi = nc.copy()
    prov = {
        p.name: sha(p)
        for p in [args.geneinfo, chrom, samples, states_path, alts_path, th_path, prior]
    }
    manifest = {
        line.split()[1]: line.split()[0]
        for line in (src / "stream_SHA256SUMS.txt").read_text().splitlines()
    }
    tok = np.full(int(names.token_id.max()) + 1, -1, dtype=np.int32)
    tok[genes.token_id] = np.arange(len(genes))
    for plate in ("plate6", "plate14"):
        directory = src / "stream/plates" / plate
        for fn in ["pseudobulk_counts.parquet", "pseudobulk_groups.tsv"]:
            path = directory / fn
            h = sha(path)
            assert h == manifest[str(path.relative_to(src / "stream"))]
            prov[str(path.relative_to(src))] = h
        n = (
            pd.read_csv(directory / "pseudobulk_groups.tsv", sep="\t")
            .groupby(["sample", "cell_line_id"], as_index=False)[["n_cells", "umi"]]
            .sum()
        )
        ii = n["sample"].map(sid)
        jj = n.cell_line_id.map({v: i for i, v in enumerate(lineids)})
        nc[ii, jj] = n.n_cells
        umi[ii, jj] = n.umi
        f = pq.ParquetFile(directory / "pseudobulk_counts.parquet")
        for b in f.iter_batches(
            batch_size=800000, columns=["sample", "cell_line_id", "token_id", "count"]
        ):
            idx = tok[b["token_id"].to_numpy()]
            ix = np.flatnonzero(idx >= 0)
            sel = b.take(pa.array(ix))
            sw = pc.index_in(sel["sample"], pa.array(s["sample"].tolist())).to_numpy()
            ln = pc.index_in(sel["cell_line_id"], pa.array(lineids)).to_numpy()
            np.add.at(out, (sw, ln, idx[ix]), sel["count"].to_numpy())
        assert (
            nc[s.index[s.plate == plate]].sum()
            == json.loads((directory / "MERGE_OK.json").read_text())["cells"]
        )
        print(
            plate,
            "cells",
            nc[s.index[s.plate == plate]].sum(),
            "selected UMIs",
            out[s.index[s.plate == plate]].sum(),
            flush=True,
        )
    assert np.all(out.sum(axis=2) <= umi)
    np.savez_compressed(HERE / "counts.npz", gene_counts=out, cells=nc, umi=umi)
    (HERE / "sources.json").write_text(
        json.dumps(
            {
                "sha256": prov,
                "n_broad": 978,
                "n_landmarks": int(genes.landmark.sum()),
                "features_sha256": sha(HERE / "features.tsv"),
                "counts_sha256": sha(HERE / "counts.npz"),
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
