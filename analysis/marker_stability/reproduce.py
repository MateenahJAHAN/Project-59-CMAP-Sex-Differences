"""Repeat the treated-well, repeat-plate and donor-label analyses from the published counts."""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def compare(expected, observed, key):
    if isinstance(expected, dict):
        if type(observed) is not dict or expected.keys() != observed.keys():
            raise ValueError(f"Fields changed at {key}")
        for k, v in expected.items():
            compare(v, observed[k], f"{key}.{k}")
    elif isinstance(expected, list):
        if not isinstance(observed, list) or len(expected) != len(observed):
            raise ValueError(f"Length changed at {key}")
        for i, v in enumerate(expected):
            compare(v, observed[i], f"{key}[{i}]")
    elif isinstance(expected, float) and isinstance(observed, (int, float)):
        if not math.isclose(expected, observed, abs_tol=1e-10, rel_tol=1e-10):
            raise ValueError(f"{key}: {expected} vs {observed}")
    elif observed != expected:
        raise ValueError(f"{key}: {expected!r} vs {observed!r}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    dest = args.work.resolve() / "marker_stability"
    if dest.is_relative_to(HERE.parents[1]):
        parser.error("Choose a work directory outside the checkout")
    files = {"": ["analyse.py", "figure.py"],
             "marker_transfer": ["analyze.py", "aggregate.py", "figure.py", "reference_cutoffs.tsv",
                                 "reference_lines.tsv", "sample_key.tsv", "well_counts.tsv.gz"],
             "genomics_extension": ["analyse.py", "figure.py", "check.py", "sources.json",
                                     "counts.npz", "wells.tsv", "lines.tsv", "features.tsv"],
             "donor_response": ["analyze.py"]}
    for name, inputs in files.items():
        (dest / name).mkdir(parents=True, exist_ok=True)
        for filename in inputs:
            shutil.copy2(HERE / name / filename, dest / name / filename)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", MPLBACKEND="Agg",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    jobs = [("analyse.py", "--root", HERE / "inputs"),
            (HERE / "verify.py", "--work", args.work.resolve(), "--out", dest / "results"), ("figure.py",),
            ("marker_transfer/analyze.py",), ("marker_transfer/figure.py",),
            ("genomics_extension/analyse.py",), ("genomics_extension/check.py",),
            ("genomics_extension/figure.py",), ("donor_response/analyze.py",)]
    for name, *args in jobs:
        print(f"Running {name}", flush=True)
        cmd = [sys.executable, "-B", str(dest / name), *map(str, args)]
        if name == "genomics_extension/check.py":
            with (dest / "genomics_extension/independent_check.json").open("w") as report:
                subprocess.run(cmd, cwd=dest, env=env, check=True, stdout=report)
        else:
            subprocess.run(cmd, cwd=dest, env=env, check=True)
    for name in ("results/summary.json", "results/independent_verification.json",
                 "marker_transfer/results.json", "genomics_extension/results.json",
                 "genomics_extension/independent_check.json", "donor_response/results.json"):
        compare(json.loads((HERE / name).read_text()), json.loads((dest / name).read_text()), name)
        print(f"Verified {name}")


if __name__ == "__main__":
    main()
