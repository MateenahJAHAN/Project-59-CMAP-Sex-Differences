#!/usr/bin/env python3
"""Step 1: build compact derived tables from the raw downloads.

Outputs (data/A/derived unless noted):
  results/A/A_marker_feature_space.tsv
  ccle_markers.tsv        one row per CCLE column, log2(RPKM+1) per marker
  ccle_columns.tsv        all CCLE column names
  sanger_markers.tsv      one row per Sanger model (data_source == Sanger), log2(TPM+1)
  sanger_source_counts.tsv
  cellosaurus_compact.tsv AC, ID, synonyms, sex, human flag
  trt_cp_counts.tsv       trt_cp signatures per cell_iname
  p59_signatures.tsv      Project-59 signature set (trt_cp, M/F, top-50 drugs)
  dmso_samples.tsv        DMSO vehicle samples with qc_pass == 1
"""
import gzip
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
from common import (DERIVED, MARKERS, P59, RAW, RES, XIST, Y_GENES, log, timed)  # noqa: E402


def marker_feature_space():
    g = pd.read_csv(os.path.join(RAW, "geneinfo_beta.txt"), sep="\t", dtype={"gene_id": str})
    assert g["gene_id"].is_unique and len(g) == 12328, len(g)
    assert set(g["feature_space"]) == {"landmark", "best inferred", "inferred"}
    m = g[g["gene_symbol"].isin(MARKERS)][["gene_symbol", "gene_id", "ensembl_id", "feature_space"]]
    assert len(m) == len(MARKERS) and m["gene_symbol"].is_unique
    m = m.set_index("gene_symbol").loc[MARKERS].reset_index()
    m.to_csv(os.path.join(RES, "A_marker_feature_space.tsv"), sep="\t", index=False)
    counts = g["feature_space"].value_counts().rename_axis("feature_space").reset_index(name="n_genes")
    counts.to_csv(os.path.join(RES, "A_geneinfo_feature_space_counts.tsv"), sep="\t", index=False)
    log("extract", f"markers feature_space: {dict(zip(m.gene_symbol, m.feature_space))}")
    return m


def ccle(markers):
    ens = dict(zip(markers["ensembl_id"], markers["gene_symbol"]))
    path = os.path.join(RAW, "CCLE_RNAseq_genes_rpkm_20180929.gct.gz")
    rows = {}
    with gzip.open(path, "rt") as fh:
        version = fh.readline().strip()
        dims = fh.readline().strip().split("\t")
        header = fh.readline().rstrip("\n").split("\t")
        n = 0
        for line in fh:
            n += 1
            name, rest = line.split("\t", 1)
            e = name.split(".")[0]
            if e in ens:
                parts = line.rstrip("\n").split("\t")
                assert ens[e] not in rows, f"duplicate CCLE row for {e}"
                rows[ens[e]] = (parts[1], np.array(parts[2:], dtype=float))
    assert version == "#1.2", version
    assert int(dims[0]) == n and int(dims[1]) == len(header) - 2, (dims, n, len(header))
    assert header[:2] == ["Name", "Description"]
    cols = header[2:]
    assert len(set(cols)) == len(cols)
    missing = [s for s in MARKERS if s not in rows]
    log("extract", f"CCLE: {n} genes x {len(cols)} lines; markers found {sorted(rows)}; missing {missing}")
    df = pd.DataFrame({"ccle_name": cols})
    for s in MARKERS:
        if s in rows:
            desc, v = rows[s]
            assert desc == s, (desc, s)
            df[f"{s}_rpkm"] = v
            df[s] = np.log2(v + 1.0)
    df.to_csv(os.path.join(DERIVED, "ccle_markers.tsv"), sep="\t", index=False)
    pd.DataFrame({"ccle_name": cols}).to_csv(os.path.join(DERIVED, "ccle_columns.tsv"), sep="\t", index=False)
    with open(os.path.join(DERIVED, "ccle_meta.json"), "w") as fh:
        json.dump({"gct_version": version, "n_genes": n, "n_lines": len(cols), "markers_missing": missing}, fh)


def sanger(markers):
    ens = dict(zip(markers["ensembl_id"], markers["gene_symbol"]))
    path = os.path.join(RAW, "rnaseq_latest.csv.gz")
    use = ["dataset_id", "model_id", "ensembl_gene_id", "gene_symbol", "model_name", "rsem_tpm",
           "data_source", "duplicate"]
    keep = []
    total = 0
    src_counts = {}
    model_src = set()
    for ch in pd.read_csv(path, usecols=use, chunksize=2_000_000, dtype={"dataset_id": str}):
        total += len(ch)
        for k, v in ch["data_source"].value_counts(dropna=False).items():
            src_counts[k] = src_counts.get(k, 0) + int(v)
        model_src.update(zip(ch["model_id"], ch["data_source"]))
        keep.append(ch[ch["ensembl_gene_id"].isin(ens)])
    k = pd.concat(keep, ignore_index=True)
    k["gene"] = k["ensembl_gene_id"].map(ens)
    log("extract", f"Sanger rows read {total}; marker rows {len(k)}; data_source counts {src_counts}")
    ms = pd.DataFrame(sorted(model_src), columns=["model_id", "data_source"])
    ms.to_csv(os.path.join(DERIVED, "sanger_model_sources.tsv"), sep="\t", index=False)
    pd.DataFrame([{"data_source": a, "rows": b} for a, b in src_counts.items()]).to_csv(
        os.path.join(DERIVED, "sanger_source_counts.tsv"), sep="\t", index=False)
    k.to_csv(os.path.join(DERIVED, "sanger_marker_rows_all_sources.tsv"), sep="\t", index=False)
    pivot_source(k, "Sanger", "sanger_markers.tsv")
    # Deviation 1 (see report): Broad-sourced rows kept as a fallback source only
    pivot_source(k, "Broad", "cmp_broad_markers.tsv")


def pivot_source(k, source, outname):
    s = k[k["data_source"] == source].copy()
    s["log2tpm"] = np.log2(s["rsem_tpm"].astype(float) + 1.0)
    # average across datasets per model and gene
    n_ds = s.groupby("model_id")["dataset_id"].nunique().rename("n_sanger_datasets")
    w = s.pivot_table(index="model_id", columns="gene", values="log2tpm", aggfunc="mean")
    w = w.reindex(columns=[g for g in MARKERS if g in w.columns])
    w = w.join(n_ds)
    names = s.groupby("model_id")["model_name"].first()
    w.insert(0, "sanger_model_name", names)
    w = w.reset_index()
    assert w["model_id"].is_unique
    w.to_csv(os.path.join(DERIVED, outname), sep="\t", index=False)
    log("extract", f"Sanger file, data_source=={source}: models with marker data: {len(w)}; columns {list(w.columns)}")


def cellosaurus():
    path = os.path.join(RAW, "cellosaurus.txt")
    recs = []
    cur = {}
    version = None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if version is None and line.startswith(" Version:"):
                version = line.strip()
            tag = line[:5]
            val = line[5:].rstrip("\n")
            if tag == "ID   ":
                cur = {"ID": val, "AC": None, "SY": "", "SX": "", "human": False}
            elif tag == "AC   ":
                cur["AC"] = val
            elif tag == "SY   ":
                cur["SY"] = val
            elif tag == "SX   ":
                cur["SX"] = val
            elif tag == "OX   ":
                if "NCBI_TaxID=9606;" in val:
                    cur["human"] = True
            elif line.startswith("//"):
                if cur.get("AC"):
                    recs.append(cur)
                cur = {}
    df = pd.DataFrame(recs)
    assert df["AC"].is_unique
    df.to_csv(os.path.join(DERIVED, "cellosaurus_compact.tsv"), sep="\t", index=False)
    with open(os.path.join(DERIVED, "cellosaurus_meta.json"), "w") as fh:
        json.dump({"version_line": version, "n_entries": len(df)}, fh)
    log("extract", f"Cellosaurus entries {len(df)}; version line: {version}")


def siginfo():
    cell = pd.read_csv(os.path.join(RAW, "cellinfo_beta.txt"), sep="\t")
    assert cell["cell_iname"].is_unique and len(cell) == 240
    sig = pd.read_csv(os.path.join(RAW, "siginfo_beta.txt"), sep="\t",
                      usecols=["sig_id", "pert_type", "cell_iname", "cmap_name"], low_memory=False)
    assert sig["sig_id"].is_unique
    assert "trt_cp" in set(sig["pert_type"])
    cp = sig[sig["pert_type"] == "trt_cp"]
    cnt = cp.groupby("cell_iname").size().rename("n_trt_cp_signatures").reset_index()
    not_in_cellinfo = sorted(set(cnt["cell_iname"]) - set(cell["cell_iname"]))
    log("extract", f"siginfo rows {len(sig)}, trt_cp {len(cp)}, lines with trt_cp {len(cnt)}; "
                   f"trt_cp cell_iname not in cellinfo: {not_in_cellinfo}")
    cnt.to_csv(os.path.join(DERIVED, "trt_cp_counts.tsv"), sep="\t", index=False)
    # Project-59 signature set (Phase 5 notebook cell 6 + top50 drugs)
    top50 = json.load(open(os.path.join(P59, "Phase5_PerSignature_Categories", "reference", "top50_names.json")))
    assert len(top50) == 50 and len(set(top50)) == 50
    n0 = len(sig)
    m = sig.merge(cell[["cell_iname", "donor_sex", "cell_type", "cell_lineage"]], on="cell_iname", how="left")
    assert len(m) == n0
    known = m[(m["pert_type"] == "trt_cp") & (m["donor_sex"].isin(["M", "F"]))]
    log("extract", f"sig_known rows {len(known)} (Project-59 dataset_facts.csv says 596,117)")
    s50 = known[known["cmap_name"].isin(top50)].copy()
    log("extract", f"Project-59 50-drug signatures: {len(s50)}; by sex {s50['donor_sex'].value_counts().to_dict()}")
    s50.to_csv(os.path.join(DERIVED, "p59_signatures.tsv"), sep="\t", index=False)
    pd.DataFrame([{"item": "sig_known_rows", "value": len(known)},
                  {"item": "p59_50drug_signatures", "value": len(s50)},
                  {"item": "p59_50drug_F", "value": int((s50.donor_sex == "F").sum())},
                  {"item": "p59_50drug_M", "value": int((s50.donor_sex == "M").sum())},
                  {"item": "siginfo_rows", "value": len(sig)},
                  {"item": "siginfo_trt_cp_rows", "value": len(cp)}]).to_csv(
        os.path.join(RES, "A_signature_set_counts.tsv"), sep="\t", index=False)


def instinfo():
    use = ["sample_id", "pert_type", "cmap_name", "cell_iname", "qc_pass", "pert_id"]
    ins = pd.read_csv(os.path.join(RAW, "instinfo_beta.txt"), sep="\t", usecols=use, low_memory=False)
    assert ins["sample_id"].is_unique
    v = ins[ins["pert_type"] == "ctl_vehicle"]
    assert "DMSO" in set(v["cmap_name"])
    assert set(v["qc_pass"].dropna().unique()) <= {0.0, 1.0}
    d = v[(v["cmap_name"] == "DMSO") & (v["qc_pass"] == 1)]
    log("extract", f"instinfo rows {len(ins)}; ctl_vehicle {len(v)}; DMSO {int((v.cmap_name == 'DMSO').sum())}; "
                   f"DMSO qc_pass==1 {len(d)}")
    d.to_csv(os.path.join(DERIVED, "dmso_samples.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    if sys.argv[1:] == ["pivot_only"]:
        k = pd.read_csv(os.path.join(DERIVED, "sanger_marker_rows_all_sources.tsv"), sep="\t", dtype={"dataset_id": str})
        pivot_source(k, "Sanger", "sanger_markers.tsv")
        pivot_source(k, "Broad", "cmp_broad_markers.tsv")
        sys.exit(0)
    what = sys.argv[1:] or ["genes", "ccle", "cellosaurus", "siginfo", "instinfo", "sanger"]
    with timed("extract:geneinfo"):
        m = marker_feature_space()
    if "ccle" in what:
        with timed("extract:ccle"):
            ccle(m)
    if "cellosaurus" in what:
        with timed("extract:cellosaurus"):
            cellosaurus()
    if "siginfo" in what:
        with timed("extract:siginfo"):
            siginfo()
    if "instinfo" in what:
        with timed("extract:instinfo"):
            instinfo()
    if "sanger" in what:
        with timed("extract:sanger"):
            sanger(m)
