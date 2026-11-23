#!/usr/bin/env python3
"""Run the Project 59 analyses and build the fellowship paper."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
REPLAY = ROOT / "verification/source_replay"
THREADS = ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "OMP_NUM_THREADS")
ENV = {**os.environ, **{v: "1" for v in THREADS}, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONHASHSEED": "0", "MPLBACKEND": "Agg"}


def command(*args, cwd=ROOT):
    subprocess.run(list(map(str, args)), cwd=cwd, env=ENV, check=True)


def verify():
    for script in ("verification/check_links.py", "verification/verify_results.py", "verification/source_replay/part_a/scripts/verify.py",
                   "verification/source_replay/part_b/verify.py", "verification/source_replay/l1000/verify.py",
                   "verification/source_replay/power/replay_data/replay_s3.py"):
        command(sys.executable, "-B", ROOT / script)
    from collections import Counter
    fields = {"source", "source_file", "row_source", "col_source", "source_line", "source_column"}
    def records(file):
        with file.open(newline="") as f:
            return Counter(tuple((k, v) for k, v in row.items() if k not in fields)
                           for row in csv.DictReader(f, delimiter="\t"))
    figures = REPLAY / "figures"
    for row in json.loads((figures / "source_figure_comparison.json").read_text())["figures"]:
        name = row["figure"]
        pdf = ROOT / "manuscript/figures" / (name + ".pdf")
        table = name + "_plotted.tsv"
        old, new = (records(figures / folder / table) for folder in ("archived_plot_data", "rebuilt_plot_data"))
        if hashlib.sha256(pdf.read_bytes()).hexdigest() != row["released_pdf_sha256"] or old != new or sum(new.values()) != row["rebuilt_plotted_rows"]:
            raise ValueError(f"Released figure or scientific records changed: {name}")
    print("Verified all seven released figure PDFs and their plotted records.", flush=True)


def environment(work, name, version, requirements):
    uv = shutil.which("uv")
    if not uv:
        raise SystemExit("Install uv to create the pinned analysis environments. See REPRODUCING.md.")
    path = work / "envs" / name
    python = path / "bin/python"
    if not python.is_file():
        path.parent.mkdir(parents=True, exist_ok=True)
        command(uv, "venv", "--python", version, path)
    elif subprocess.check_output([python, "-c", "import platform; print(platform.python_version())"], text=True).strip() != version:
        raise SystemExit(f"{path} contains a different Python version; choose a fresh work directory.")
    command(uv, "pip", "install", "--python", python, "-r", ROOT / requirements)
    if name == "fellowship":
        command(uv, "pip", "install", "--python", python, "--target", path / "phase6", "--no-deps", "gseapy==1.3.2")
    return python


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    stages = ("verify", "paper", "fellowship", "part-a", "part-b", "transfer", "claims", "l1000", "power", "paper-from-run", "follow-up", "all")
    parser.add_argument("stage", choices=stages)
    parser.add_argument("--work", type=Path, help="Data, environments and results directory outside the checkout")
    parser.add_argument("--workers", type=int, default=4, help="Parallel workers for the fellowship, LINCS and power runs")
    args = parser.parse_args()
    if args.stage == "verify":
        verify()
        return
    if not args.work or args.work.expanduser().resolve().is_relative_to(ROOT) or args.workers < 1:
        parser.error("Supply --work outside this repository and at least one worker.")
    if args.stage in ("all", "power") and args.workers > 8:
        parser.error("The power simulation supports 1–8 workers.")
    work = args.work.expanduser().resolve()
    work.mkdir(parents=True, exist_ok=True)
    selected = ("fellowship", "part-b", "part-a", "transfer", "claims", "l1000", "power", "paper-from-run", "verify") if args.stage == "all" else (args.stage,)
    for stage in selected:
        print(f"Starting {stage}", flush=True)
        if stage == "verify":
            verify()
        elif stage in ("paper", "paper-from-run"):
            py = environment(work, "plots", "3.12.11", "manuscript/requirements.txt")
            extra = []
            if stage == "paper-from-run":
                command(py, ROOT / "manuscript/from_run.py", "--work", work)
                extra = ["--inputs", work / "plot-inputs"]
            command(py, ROOT / "manuscript/build.py", "--out", work / stage, *extra)
        elif stage == "claims":
            command(sys.executable, REPLAY / "claim_audit/run.py", "--work", work, "--output", work / "results/claims")
        elif stage == "l1000":
            py = environment(work, stage, "3.12.11", "verification/source_replay/l1000/requirements.txt")
            states = work / "results/A/A_cellline_sex_states.tsv"
            if not states.is_file():
                raise SystemExit("Run part-b first to prepare the bulk RNA states and Broad metadata.")
            path = work / "results/l1000"
            command(py, REPLAY / "l1000/run.py", "--data", work / "sources", "--states", states, "--output", path)
            command(py, REPLAY / "l1000/verify.py", "--states", states, "--output", path)
        elif stage == "fellowship":
            py = environment(work, stage, "3.12.11", "fellowship/requirements.txt")
            command(py, ROOT / "fellowship/run.py", "--data", work / "fellowship-data",
                    "--out", work / "fellowship-results", "--workers", args.workers)
        elif stage in ("part-a", "transfer", "follow-up"):
            py = environment(work, "atlas", "3.12.11", "analysis/tahoe/requirements.txt")
            if stage in ("part-a", "transfer"):
                step = "all" if stage == "part-a" else "transfer-s16"
                command(py, ROOT / "analysis/part_a/run.py", "--work", work, "--stage", step)
            else:
                command(py, ROOT / "analysis/marker_stability/reproduce.py", "--work", work)
        else:
            py = environment(work, "lincs", "3.13.12", "analysis/lincs/requirements.txt")
            script = "analysis/lincs/reproduce.py" if stage == "part-b" else "verification/source_replay/power/run.py"
            command(py, ROOT / script, "--work", work, "--workers", args.workers)


if __name__ == "__main__":
    main()
