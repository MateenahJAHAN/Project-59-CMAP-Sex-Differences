"""Reproduce the RNA state analyses in a separate working directory."""
import argparse
import csv
import datetime as dt
import json
import os
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
REPLAY = REPO / "verification/source_replay"
sys.path.insert(0, str(REPLAY / "part_a/scripts"))
from fetch import download, sha256_file


def acquire(row, target):
    target = Path(target)
    if target.exists():
        if target.stat().st_size != int(row["bytes"]) or sha256_file(target) != row["sha256"]:
            raise ValueError(f"Source snapshot mismatch at {target}. See {row['url']}.")
        return
    print(f"Downloading {target.name}", flush=True)
    url = row["url"]
    if "figshare.com/" in url:
        import time
        url += ("&" if "?" in url else "?") + f"download=1&replay={int(time.time())}"
    download(url, str(target), row["sha256"], int(row["bytes"]))


def execute(script, work, *args, **env):
    print(f"Running {Path(script).name}", flush=True)
    threads = {name: "1" for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                                     "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")}
    subprocess.run([sys.executable, str(script), *map(str, args)], cwd=work,
                   env={**os.environ, **threads, "P59_WORK": str(work), "P59_DATA": str(work),
                        "P59_ANALYSIS": str(work), **env}, check=True)


def labels_c(work):
    for row in json.loads((REPLAY / "bulk_RNA/sources.json").read_text()):
        acquire(row, work / "sources" / row["file"])
    for name in ("extract_expression.py", "build_labels.py"):
        execute(HERE / "bulk" / name, work)


SOURCES = json.loads((HERE / "references/sources.json").read_text())
SCRIPTS = REPLAY / "part_a/scripts"
TAHOE = HERE.parent / "tahoe"


def link_input(source, target):
    source = source.resolve(strict=True)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.samefile(source):
        return
    if target.exists() or target.is_symlink():
        raise ValueError(f"Conflicting workflow input {target}")
    try:
        target.symlink_to(source)
    except FileExistsError:
        if not target.samefile(source):
            raise


def prepare_labels(work):
    for label, name in (("A", "A_cellline_sex_states.tsv"), ("C", "C_cellline_sex_labels.tsv")):
        target = work / "results" / label / name
        if not target.is_file():
            if label == "C":
                labels_c(work)
            else:
                execute(HERE.parent / "lincs/reproduce.py", work, "--work", work, "--stage", "labels-a")
        link_input(target, work / "inputs" / name)
    link_input(HERE / "references/tahoe_pilot_states.tsv", work / "inputs/tahoe_pilot_states.tsv")


def tahoe_metadata(work):
    meta = work / "tahoe/index/meta"
    for row in SOURCES["tahoe_metadata"]:
        acquire(row, meta / row["file"])


def reduce_tahoe(work):
    tahoe_metadata(work)
    obs = json.loads((REPLAY / "part_a/receipts/tahoe_observations.json").read_text())
    acquire(dict(url=obs["url"], bytes=obs["file_bytes"], sha256=obs["sha256"]),
            work / "tahoe/index/meta/obs_metadata.parquet")
    for name in ("t0_obs_compact", "t1_footer_fetch", "t1_build_index", "t2_dmso_select", "t4_dmso_pass"):
        execute(HERE / "reduce" / (name + ".py"), work)


def annotations(work, out):
    out.mkdir(parents=True, exist_ok=True)
    timestamp = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    fields = ["item", "url", "local_path", "bytes", "sha256", "accessed_utc", "notes"]
    with (out / "provenance_downloads.tsv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for item in SOURCES["annotations"]:
            dest = out / "downloads" / item["file"]
            if dest.exists():
                acquire(item, dest)
            else:
                source = work / "sources" / Path(item["file"]).name
                acquire(item, source)
                link_input(source, dest)
            writer.writerow({k: item[k] for k in ("item", "url", "bytes", "sha256")} |
                            dict(local_path=str(dest), accessed_utc=timestamp, notes="Verified source snapshot"))


def tahoe_call(work):
    prepare_labels(work)
    tahoe_metadata(work)
    out = work / "results/tahoe"
    out.mkdir(parents=True, exist_ok=True)
    # An interrupted caller can restart in a new directory without changing a completed result.
    import tempfile
    fresh = Path(tempfile.mkdtemp(dir=out, prefix="run-"))
    for folder in ("logs", "figures"):
        (fresh / folder).mkdir()
    annotations(work, fresh)
    for name in ("01_gene_chrom", "02_call_states", "03_sensitivity", "04_agreement_audit",
                 "05_ambient_foreign", "06_mosaic", "07_identity", "08_mosaic_rule", "09_figures", "10_gene_level"):
        execute(TAHOE / (name + ".py"), work, P59_PARTA_OUT=str(fresh))
    compare_states(fresh / "partA_local_states.tsv", "cellosaurus_id")
    pointer = out / "latest"
    temp = out / (fresh.name + ".link")
    temp.symlink_to(fresh.name, target_is_directory=True)
    temp.replace(pointer)
    print(f"Tahoe caller output: {pointer}", flush=True)


def compare_states(actual, *keys):
    import pandas as pd
    reference = REPLAY / "part_a/tables" / actual.name
    calculated = pd.read_csv(actual, sep="\t").set_index(list(keys)).sort_index()
    published = pd.read_csv(reference, sep="\t").set_index(list(keys)).sort_index()
    # The single Cellosaurus accession timestamp belongs to this run.
    if actual.name == "partA_local_states.tsv":
        for table in (calculated, published):
            table["annotation_source"] = table["annotation_source"].str.replace(
                r"(Cellosaurus entry CVCL_C466, SX line \(api.cellosaurus.org, accessed )[^)]+", r"\1<date>", regex=True)
    extra = calculated.columns.difference(published.columns).tolist()
    if extra:
        print(f"Additional columns from full pipeline in {actual.name}: {extra}")
    pd.testing.assert_frame_equal(calculated[published.columns], published, check_exact=False, rtol=1e-10, atol=1e-10)
    print(f"Matched {len(published)} rows and all {len(published.columns)} measured/annotation columns in {actual.name}.")


def external_reduce(work):
    execute(SCRIPTS / "prepare_emerald.py", work, work)
    execute(SCRIPTS / "process_emerald.py", work, work, "2")
    manifest = json.loads((REPLAY / "part_a/receipts/source_verification.json").read_text())
    for row in manifest["sciplex_GEO_files"]:
        source = row["source"]
        acquire(dict(url=source, bytes=row["size"], sha256=row["sha256"]),
                work / "data/F1/sciplex" / source.rsplit("/", 1)[-1])
    for tag in ("sciPlex2", "sciPlex3", "sciPlex4"):
        execute(SCRIPTS / "process_sciplex.py", work, work, tag)


def external_call(work):
    import pandas as pd
    prepare_labels(work)
    root = work / "data/F1"
    shards = sorted((root / "emerald_proc").glob("meta_*.parquet"))
    if len(shards) != 116:
        raise ValueError(f"Expected metadata from 116 EmeraldBay shards; found {len(shards)}")
    frame = pd.concat([pd.read_parquet(f) for f in shards], ignore_index=True)
    if len(frame) != 1831648 or not frame.BARCODE_SUB_LIB_ID.is_unique:
        raise ValueError(f"Unexpected EmeraldBay barcode inventory: {len(frame)}")
    frame.to_parquet(root / "emerald_cell_meta.parquet", index=False)
    del frame
    for script in (HERE / "external/analyze.py", SCRIPTS / "explore_similarity.py",
                   HERE / "external/figures.py", HERE / "external/summarize.py"):
        execute(script, work, work)
    for file, key in (("F1_emerald_states.tsv", ("cellosaurus_id",)),
                      ("F1_sciplex_states.tsv", ("sample", "cell_name"))):
        compare_states(work / "results/F1" / file, *key)


def second_reduce(work):
    tahoe_metadata(work)
    dest = work / "tahoe/index/meta/obs_metadata.parquet"
    obs = json.loads((REPLAY / "part_a/receipts/tahoe_observations.json").read_text())
    acquire(dict(url=obs["url"], bytes=obs["file_bytes"], sha256=obs["sha256"]), dest)
    for folder in ("data/F2", "results/F2", "logs", "provenance"):
        (work / folder).mkdir(parents=True, exist_ok=True)
    link_input(dest.parent, work / "data/F2/meta")
    for script in ("s00_filelist", "s01_obs", "s02_index", "s02b_index_check", "s03_reduce"):
        execute(HERE / "second" / (script + ".py"), work)


def second_call(work):
    prepare_labels(work)
    pie = SOURCES["annotations"][0]
    source = work / "sources" / Path(pie["file"]).name
    acquire(pie, source)
    link_input(source, work / "data/F2/pie/descriptions.json")
    for folder in ("figures/F2", "results/F2", "logs"):
        (work / folder).mkdir(parents=True, exist_ok=True)
    for script in ("s04_states", "s05_agreement", "s06_ambient", "s07_mosaic", "s08_baseline",
                   "s09_figures", "s10_consort_hashes"):
        execute(HERE / "second" / (script + ".py"), work)
    compare_states(work / "results/F2/F2_states.tsv", "cell_line_id")


def second(work):
    second_reduce(work)
    second_call(work)


def depth_matching(work):
    import pandas as pd
    from importlib.util import module_from_spec, spec_from_file_location
    spec = spec_from_file_location("external_constants", SCRIPTS / "common.py")
    constants = module_from_spec(spec)
    spec.loader.exec_module(constants)
    frames = []
    for tag in ("sciPlex2", "sciPlex3", "sciPlex4"):
        v = pd.read_parquet(work / f"data/F1/sciplex_proc/{tag}_vehicle_cells.parquet")
        unit = v.well if tag == "sciPlex2" else v.batch
        frame = pd.DataFrame(dict(dataset=tag, cell_barcode=v.cell, line=v.line, unit=unit.astype(str),
                                  bg_batch=v.bg_batch, depth=v.depth.astype("int64")))
        for marker in constants.MARKERS:
            frame[marker] = v[marker].astype("int64")
        frames.append(frame)
    path = work / "results/F1/sciplex_vehicle_cells.tsv.gz"
    pd.concat(frames, ignore_index=True).to_csv(path, sep="\t", index=False, compression={"method":"gzip","mtime":0})
    execute(HERE / "secondary/depth_matching.py", work, work)
    out = work / "results/depth_matching"
    compare_states(out / "A3_pairs.tsv", "pair", "estimator")
    compare_states(out / "A3_replicates.tsv", "pair", "replicate")
    compare_states(out / "A3_unmatched.tsv", "dataset", "line", "marker")


def str_comparison(work):
    out = work / "results/str"
    out.mkdir(parents=True, exist_ok=True)
    execute(HERE / "secondary/parse_str.py", work, out, work / "sources/cellosaurus.txt",
            HERE / "references/accessions.tsv")
    execute(HERE / "secondary/str_comparison.py", work, work)
    compare_states(out / "figure3_supp_table_amelogenin_by_line.tsv", "accession")
    compare_states(out / "figure3_agreement_by_source.tsv", "source")


def tahoe(work):
    reduce_tahoe(work)
    tahoe_call(work)


def external(work):
    external_reduce(work)
    external_call(work)


def consensus(work):
    for name in ("consensus", "mosaic_resolution", "cross_atlas"):
        execute(HERE / "secondary" / (name + ".py"), work, work)
    out = work / "results/consensus"
    compare_states(out / "tahoe_partA_states_consensus.tsv", "cellosaurus_id")
    compare_states(out / "A2_emerald_vs_full_agreement.tsv", "comparison", "subset")
    compare_states(out / "A2_emerald_vs_full_lines.tsv", "cellosaurus_id")


def transfer_s16(work):
    for name in ("transfer_calls", "transfer"):
        execute(HERE / "secondary" / (name + ".py"), work, work)


def table_s1(work):
    execute(HERE / "secondary/table_s1.py", work, work)


def all_stages(work):
    prepare_labels(work)
    tahoe(work)
    second(work)
    external(work)
    consensus(work)
    depth_matching(work)
    str_comparison(work)
    table_s1(work)


STAGES = {"all": all_stages, "transfer-s16": transfer_s16, "table-s1": table_s1, "consensus": consensus, "labels-c": labels_c, "tahoe-reduce": reduce_tahoe, "tahoe-call": tahoe_call, "tahoe": tahoe,
          "second": second, "second-reduce": second_reduce, "second-call": second_call,
          "depth-matching": depth_matching, "str": str_comparison,
          "external-reduce": external_reduce, "external-call": external_call, "external": external}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--stage", choices=STAGES, required=True)
    args = p.parse_args()
    work = args.work.expanduser().resolve()
    if work == REPO or REPO in work.parents:
        p.error("Choose a work directory outside the Git checkout for sources, results and checkpoints.")
    work.mkdir(parents=True, exist_ok=True)
    STAGES[args.stage](work)


if __name__ == "__main__":
    main()
