#!/usr/bin/env python3
"""Stream the two RNA-seq files and keep only XIST and the eight Y genes.

Local outputs:
  results/C/C_expr_sanger_cmp_long.tsv.gz  (Cell Model Passports rnaseq_latest, rsem_tpm)
  results/C/C_expr_ccle_legacy_long.tsv.gz (CCLE legacy RPKM)
"""
import gzip
import os
import sys
import time

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import Step  # noqa: E402

ROOT = os.environ["P59_WORK"]
DATA = os.path.join(ROOT, "sources")
OUT = os.path.join(ROOT, "results", "C")
GENES = ["XIST", "RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y", "NLGN4Y"]


def sanger():
    t0 = time.time()
    keep = []
    n_rows = 0
    cols = ["dataset_id", "model_id", "gene_symbol", "ensembl_gene_id", "model_name",
            "rsem_tpm", "rsem_fpkm", "data_source", "duplicate"]
    for chunk in pd.read_csv(os.path.join(DATA, "rnaseq_latest.csv.gz"), usecols=cols,
                             chunksize=2_000_000, low_memory=False):
        n_rows += len(chunk)
        keep.append(chunk[chunk.gene_symbol.isin(GENES)])
    df = pd.concat(keep, ignore_index=True)
    assert set(df.gene_symbol) <= set(GENES)
    df.to_csv(os.path.join(OUT, "C_expr_sanger_cmp_long.tsv.gz"), sep="\t", index=False)
    print(f"sanger: scanned {n_rows} rows, kept {len(df)} rows, genes={sorted(df.gene_symbol.unique())}, "
          f"models={df.model_id.nunique()}, {time.time()-t0:.0f}s", flush=True)


def ccle():
    t0 = time.time()
    path = os.path.join(DATA, "CCLE_RNAseq_genes_rpkm_20180929.gct.gz")
    with gzip.open(path, "rt") as fh:
        assert fh.readline().startswith("#1.2")
        dims = fh.readline().strip().split("\t")
        header = fh.readline().rstrip("\n").split("\t")
        rows = []
        n = 0
        for line in fh:
            n += 1
            parts = line.rstrip("\n").split("\t")
            if parts[1] in GENES:
                rows.append(parts)
    assert n == int(dims[0]), (n, dims)
    assert len(header) == int(dims[1]) + 2
    wide = pd.DataFrame(rows, columns=header)
    long = wide.melt(id_vars=["Name", "Description"], var_name="ccle_name", value_name="rpkm")
    long["rpkm"] = long.rpkm.astype(float)
    long = long.rename(columns={"Name": "ensembl_gene_id", "Description": "gene_symbol"})
    long.to_csv(os.path.join(OUT, "C_expr_ccle_legacy_long.tsv.gz"), sep="\t", index=False)
    print(f"ccle: {n} genes scanned, kept {len(wide)} gene rows {sorted(wide.Description)}, "
          f"{long.ccle_name.nunique()} lines, {time.time()-t0:.0f}s", flush=True)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    with Step("extract expression: CCLE legacy"):
        ccle()
    with Step("extract expression: CMP rnaseq_latest"):
        sanger()
