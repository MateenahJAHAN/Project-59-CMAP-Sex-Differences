"""Check Note S17's single-plant permutation claim against archived and rerun records."""

import argparse
import csv
import gzip
import json
import math
from collections import Counter
from pathlib import Path


TERMS = ("factY", "factX")
SINGLE = ("C1", "G1")


def single_rows(rows):
    return [row for row in rows if row["scenario"] in SINGLE]


def read_aggregates(folder):
    records = []
    for name in ("power_runs.tsv.gz", "power_planted.tsv.gz"):
        with gzip.open(folder / name, "rt") as src:
            records.append(single_rows(csv.DictReader(src, delimiter="\t")))
    return records


def read_checkpoints(folder):
    import numpy as np

    calls, plants = [], []
    paths = sorted(folder.glob("rep_*.npz"))
    if not paths:
        raise ValueError(f"No completed checkpoints in {folder}")
    manifest = json.loads((folder / "run.json").read_text())
    if manifest["R"] != 1000 or manifest["B"] != 2000:
        raise ValueError("Expected the 1,000-assignment, 2,000-permutation run")
    for path in paths:
        with np.load(path, allow_pickle=False) as data:
            rows, planted = json.loads(data["records"].item())
        rep = int(path.stem.split("_")[-1])
        if len(rows) != 26 or len(planted) != 120 or any(row["rep"] != rep for row in rows + planted):
            raise ValueError(f"Incomplete or mislabeled assignment in {path}")
        calls.extend(single_rows(rows))
        plants.extend(single_rows(planted))
    return calls, plants


def key(row):
    return int(row["rep"]), row["term"], row["scenario"]


def checked_calls(calls, plants):
    runs = {key(row): row for row in calls}
    result = {}
    for row in plants:
        assignment = key(row)
        run = runs[assignment]
        p, q = float(row["p"]), float(row["q"])
        permutations = int(float(run["B"]))
        exceedances = round(p * (permutations + 1)) - 1
        total = int(float(run["n_q_lt_01"]))
        detected = int(float(run["n_planted_q_lt_01"]))
        assert permutations == 2000 and 0 <= exceedances <= permutations
        assert math.isclose(p, (1 + exceedances) / (1 + permutations), abs_tol=1e-14, rel_tol=0)
        assert int(float(run["n_planted"])) == 1 and detected == int(q < 0.1)
        assert total == detected + int(float(run["n_other_q_lt_01"]))
        assert int(float(row["n_ranked"])) == int(float(run["n_tested"]))
        assert assignment not in result
        result[assignment] = dict(compound=row["compound"], p=p, q=q, exceedances=exceedances,
                                  T_statistic_rank=int(float(row["rank"])),
                                  compounds_tested=int(float(row["n_ranked"])), total_q_lt_01=total)
    expected = {(rep, term, scenario) for rep in {k[0] for k in runs}
                for term in TERMS for scenario in SINGLE}
    assert len(calls) == len(runs) and set(runs) == set(result) == expected
    return result


def summarize(records):
    passed = Counter()
    exceptions = []
    ids = {rep for rep, _, _ in records}
    for (rep, term, scenario), data in sorted(records.items()):
        if data["q"] < 0.1:
            passed[f"{term}/{scenario}"] += 1
            if data["exceedances"] > 0:
                exceptions.append(dict(assignment=rep, term=term, scenario=scenario, **data))
    return dict(completed_assignments=len(ids), all_1000_assignments=ids == set(range(1000)),
                single_plant_records=len(records), significant_single_plants=sum(passed.values()),
                significant_by_term_scenario=dict(sorted(passed.items())),
                significant_despite_nonzero_exceedances=exceptions)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path, help="Directory with the archived power TSVs")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--checkpoints", type=Path, help="Directory containing run.json and rep_*.npz")
    group.add_argument("--final", type=Path, help="Directory with the completed rerun's power TSVs")
    args = parser.parse_args()
    archived = checked_calls(*read_aggregates(args.archive))
    assert {rep for rep, _, _ in archived} == set(range(1000))
    fresh = checked_calls(*(read_checkpoints(args.checkpoints) if args.checkpoints else read_aggregates(args.final)))
    if not set(fresh) <= set(archived):
        raise ValueError("Unexpected rerun assignment, term or scenario")
    if args.final and set(fresh) != set(archived):
        raise ValueError("The final TSVs do not contain all 1,000 assignments")
    deltas = [dict(assignment=rep, term=term, scenario=scenario, rerun=data,
                   archive=archived[rep, term, scenario])
              for (rep, term, scenario), data in sorted(fresh.items()) if data != archived[rep, term, scenario]]
    print(json.dumps(dict(archive=summarize(archived), rerun=summarize(fresh),
                          records_differing_on_matched_assignments=deltas), indent=2))
    if deltas:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
