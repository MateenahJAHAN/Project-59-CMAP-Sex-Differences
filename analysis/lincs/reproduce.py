"""Rebuild the LINCS manuscript analysis from provider files and the released study design."""
import argparse
import csv
import hashlib
import fcntl
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
from urllib.request import Request, urlopen
from urllib.error import URLError
from http.client import HTTPException

HERE = Path(__file__).resolve().parent
REPOSITORY = HERE.parents[1]
PREPARE = HERE / "prepare"
CELL = "cellinfo_beta.txt"
SIG = "siginfo_beta.txt"
GENE = "geneinfo_beta.txt"
INST = "instinfo_beta.txt"
GTF = "Homo_sapiens.GRCh38.110.gtf.gz"
CCLE = "CCLE_RNAseq_genes_rpkm_20180929.gct.gz"
MAF = "CCLE_DepMap_18q3_maf_20180718.txt"
CMP = "model_list_20240110.csv"
MUT = "mutations_summary_latest.csv.gz"
CELLO = "cellosaurus.txt"
RNA = "rnaseq_latest.csv.gz"
PRISM = "prism_secondary_dose_response_curve_parameters.csv"
COLLECTRI = "collectri_omnipath.tsv"
WORKBOOK = "mmc5.xlsx"
ENV = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS")


def file_hash(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def download(work, names):
    with (HERE / "config/source_files.tsv").open() as f:
        sources = {row["filename"]: row for row in csv.DictReader(f, delimiter="\t")}
    for name in sorted(set(names)):
        row = sources[name]
        dest = work / "sources" / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if dest.is_file():
            if file_hash(dest) != row["sha256"]:
                raise ValueError(f"{dest} differs from the study SHA-256 {row['sha256']}. Supply the matching source revision.")
            continue
        partial = dest.with_name(dest.name + ".partial")
        lock = dest.with_name(dest.name + ".lock")
        with lock.open("a") as guard:
            fcntl.flock(guard, fcntl.LOCK_EX)
            if dest.is_file() and file_hash(dest) == row["sha256"]:
                continue
            size = int(row["bytes"])
            print("Downloading", name, flush=True)
            for attempt in range(5):
                offset = partial.stat().st_size if partial.exists() else 0
                if offset == size:
                    break
                headers = {"User-Agent": "Project-59-source-replay/1.0", "Accept-Encoding": "identity"}
                if offset:
                    headers["Range"] = f"bytes={offset}-"
                url = row["url"]
                if name == PRISM:
                    # Avoid intermediaries returning an expired signed Figshare redirect.
                    url += "?download=1&replay=" + str(int(time.time()))
                try:
                    with urlopen(Request(url, headers=headers), timeout=180) as src:
                        partial_response = src.status == 206
                        if partial_response and not src.headers.get("Content-Range", "").startswith(f"bytes {offset}-"):
                            raise ValueError(f"Provider returned a different range for {name}")
                        with partial.open("ab" if offset and partial_response else "wb") as out:
                            shutil.copyfileobj(src, out, length=1024 * 1024)
                    if partial.stat().st_size == size:
                        break
                except (URLError, HTTPException, TimeoutError, OSError) as exc:
                    if attempt == 4:
                        raise
                    print(f"{name}: transfer interrupted, retry {attempt + 1}: {exc}", flush=True)
                    time.sleep(2 ** attempt)
            if partial.stat().st_size != size or file_hash(partial) != row["sha256"]:
                partial.unlink()
                raise ValueError(f"Source size or SHA-256 differs for {name}. Obtain the revision with SHA-256 {row['sha256']} and place it at {dest}.")
            partial.replace(dest)


def link(src, dst):
    if not src.exists():
        raise FileNotFoundError(src)
    if dst.exists() and dst.samefile(src):
        return
    if dst.is_symlink() and (not dst.exists() or
            (src.is_file() and dst.is_file() and file_hash(src) == file_hash(dst))):
        # A work area can move between checkouts of the same released inputs.
        with tempfile.TemporaryDirectory(dir=dst.parent) as tmp:
            replacement = Path(tmp) / "link"
            replacement.symlink_to(src, target_is_directory=src.is_dir())
            replacement.replace(dst)
        return
    if dst.exists():
        if src.is_file() and dst.is_file() and file_hash(src) == file_hash(dst):
            return
        raise FileExistsError(f"Existing path points elsewhere: {dst}")
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.symlink_to(src, target_is_directory=src.is_dir())


def inputs(work, relative, names):
    download(work, names)
    for name in names:
        link(work / 'sources' / name, work / relative / name)


def run(work, script, *args, workers=None):
    print(f"Running {script.name}", flush=True)
    env = {**os.environ, **{k: "1" for k in ENV}, "PYTHONHASHSEED": "0",
           "P59_LINCS_WORK": str(work), "P59_PUBLIC_REPOSITORY": str(REPOSITORY)}
    if workers:
        env["P59_CAL_WORKERS"] = str(workers)
        env["P59_CALIBRATION_CACHE"] = str(work / "calibration_checkpoints")
    subprocess.run([sys.executable, str(script), *map(str, args)], cwd=work, env=env, check=True)


def labels_a(work):
    inputs(work, "data/A/raw", [CELL, GENE, SIG, CCLE, CMP, CELLO, RNA])
    run(work, PREPARE / "A/01_extract.py", "ccle", "cellosaurus", "siginfo", "sanger")
    run(work, PREPARE / "A/02_states.py")


def level5(work):
    inputs(work, "data/D/cmap", [GENE, SIG])
    run(work, PREPARE / "D/01_inspect_gctx.py")
    run(work, PREPARE / "D/02_tensors.py")


def wells(work):
    inputs(work, "data/E/meta", [CELL, GENE, SIG, INST])
    inputs(work, "data/E/cmp", [CMP, MUT])
    data = work / "data/E"
    run(work, PREPARE / "E/inspect_gctx.py", data / "gctx_layout.json")
    run(work, PREPARE / "E/fetch_ids.py", data / "ids")
    for name in ["01_wells", "02_lines", "03_read", "04_responses"]:
        run(work, PREPARE / "E" / (name + ".py"), work)


def design(work):
    inputs(work, "data/D/cmap", [CELL, SIG])
    inputs(work, "data/D/cmp", [CMP, MUT])
    inputs(work, "data/D/cello", [CELLO])
    inputs(work, "data/D/ccle", [CCLE, MAF])
    inputs(work, "data/D/reg", [COLLECTRI, WORKBOOK])
    run(work, PREPARE / "D/01_line_table.py")
    run(work, PREPARE / "D/03_regulons.py")
    download(work, [GTF])
    run(work, PREPARE / "D/annotations.py", work)
    link(work / "results/D", work / "data/D2/inputs/taskD")
    run(work, PREPARE / "D2/02_sets.py", work)
    for group in ("D", "D2"):
        run(work, PREPARE / group / "04_pseudolabels.py")
    run(work, HERE / "prepare_model.py", work)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work", type=Path, required=True, help="Data and output directory outside the repository.")
    stages = ["download", "labels-a", "labels-c", "level5", "wells", "design", "model", "exploratory", "figure4", "check", "all"]
    p.add_argument("--stage", choices=stages, default="all")
    p.add_argument("--workers", type=int, default=4)
    a = p.parse_args()
    if a.workers < 1:
        p.error("workers must be at least one")
    work = a.work.resolve()
    if work.is_relative_to(REPOSITORY):
        p.error("choose a work directory outside the checkout for the provider data")
    for d in ["sources", "logs", "data/E", "results/A", "results/D", "results/D2", "results/E"]:
        (work / d).mkdir(parents=True, exist_ok=True)
    steps = stages[:-1] if a.stage == "all" else [a.stage]
    for step in steps:
        if step == "download":
            download(work, [CELL, SIG, GENE, INST, CCLE, MAF, CMP, MUT, CELLO, RNA, PRISM, COLLECTRI, WORKBOOK, GTF])
        elif step == "labels-a":
            labels_a(work)
        elif step == "labels-c":
            run(work, HERE.parent / "part_a/run.py", "--work", work, "--stage", "labels-c")
        elif step == "level5":
            level5(work)
        elif step == "wells":
            wells(work)
        elif step == "design":
            design(work)
        elif step == "model":
            download(work, [PRISM, SIG])
            run(work, HERE / "model.py", "--work", work, "--workers", a.workers)
        elif step == "exploratory":
            for name in ("exploratory_X1_vehicle", "exploratory_X2_idarubicin", "exploratory_all", "diag_X1b_common_vector"):
                run(work, HERE / "exploratory" / (name + ".py"), workers=a.workers)
        elif step == "figure4":
            scripts = REPOSITORY / "verification/source_replay/part_b"
            for name in ("figure4_nulls", "figure4_intervals"):
                run(work, scripts / (name + ".py"), "--work", work, "--output", work / "results/figure4")
            run(work, HERE / "render_figure4.py", work)
        else:
            run(work, HERE / "check_model.py", "--work", work)


if __name__ == "__main__":
    main()
