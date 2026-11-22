#!/usr/bin/env python3
"""Check the released Project 59 result tables with the Python standard library."""

import csv
import gzip
import hashlib
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path
import sys

from verify_phases import verify_fellowship_phases


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
RESULT_FILES = {
    "part_a.tsv", "part_b_h1.tsv", "part_b_h3.tsv", "part_b_compounds.tsv",
    "part_b_controls.tsv", "phase6_counts.tsv", "phase6_comparison_labels.tsv.gz",
}
TERMS = {"Y", "XIST"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_table(name, columns):
    opener = gzip.open if name.endswith(".gz") else open
    with opener(RESULTS / name, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames == columns.split(), f"Unexpected columns in {name}")
        for row in reader:
            require(None not in row and all(value is not None for value in row.values()),
                    f"Malformed row in {name}")
            yield row


def verify_checksums():
    expected = {f"results/{name}" for name in RESULT_FILES}
    verified = set()
    for line in (HERE / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        checksum, path = line.split("  ", 1)
        require(path in expected and path not in verified, f"Unexpected checksum path: {path}")
        result = HERE / path
        require(not result.is_symlink(), f"Result must be a regular file: {path}")
        require(hashlib.sha256(result.read_bytes()).hexdigest() == checksum,
                f"Checksum mismatch: {path}")
        verified.add(path)
    require(verified == expected, "Manifest does not cover every result table")


def verify_part_a():
    metrics = {}
    for row in read_table("part_a.tsv", "metric count total unit"):
        metric = row["metric"]
        require(metric not in metrics, f"Duplicate Part A metric: {metric}")
        count, total = int(row["count"]), int(row["total"])
        require(row["unit"] in ("cells", "lines") and 0 <= count <= total and total > 0,
                f"Invalid Part A count: {metric}")
        metrics[metric] = count, total

    def pair(metric, count, total):
        require(metrics[metric] == (count, total), f"Unexpected Part A count: {metric}")

    pair("vehicle_cells_read", 2_205_786, 2_205_786)
    pair("vehicle_cells_retained", 2_175_916, 2_205_786)
    pair("tahoe_cell_lines", 50, 50)
    states = {"XY": 18, "XX_like": 5, "Y_lost": 14, "XIST_lost": 12, "discordant": 1}
    for state, count in states.items():
        pair(f"state_{state}", count, 50)
    require(sum(states.values()) == metrics["tahoe_cell_lines"][0], "State counts do not sum to 50")
    pair("annotation_Male", 33, 50)
    pair("annotation_Female", 16, 50)
    pairs = {
        "donor_annotation_pattern_mismatch": (27, 49),
        "external_RNA_seq_expression_pattern_agreement": (47, 49),
        "external_RNA_seq_annotation_dependent_state_agreement": (46, 49),
        "EmeraldBay_agreement_confident_in_both": (40, 41),
        "EmeraldBay_agreement_all_shared": (42, 47),
        "Amelogenin_Y_agreement_post_hoc": (48, 50),
    }
    for metric, counts in pairs.items():
        pair(metric, *counts)
    require(len(metrics) == 16, "Unexpected number of Part A metrics")
    print("Part A: state counts and comparison denominators verified.")


def indexed_by_term(table, columns):
    rows = list(read_table(table, columns))
    indexed = {row["term"]: row for row in rows}
    require(len(rows) == len(indexed) == 2 and set(indexed) == TERMS,
            f"Unexpected terms in {table}")
    return indexed


def probability(value):
    number = Decimal(value)
    require(number.is_finite() and 0 <= number <= 1, f"Invalid probability: {value}")
    return number


def bh_adjusted(ordered_p):
    running = Decimal(1)
    count = len(ordered_p)
    q_values = []
    for rank in range(count, 0, -1):
        running = min(running, ordered_p[rank - 1] * count / rank)
        q_values.append(running)
    return reversed(q_values)


def verify_part_b():
    panels = indexed_by_term(
        "part_b_h1.tsv",
        "term n_lines n_with_indicator n_without_indicator n_lineages n_permutations "
        "n_tested_compounds panel_statistic panel_p n_nominal_p_lt_0_05 n_bh_q_lt_0_10 "
        "vehicle_p vehicle_rank_if_inserted",
    )
    sizes = indexed_by_term(
        "part_b_h3.tsv",
        "term n_lines residual_df estimate_pct lo90_pct hi90_pct lo95_pct hi95_pct "
        "probability_interval_inside_margin",
    )
    groups = {term: [] for term in TERMS}
    seen = set()
    for row in read_table("part_b_compounds.tsv", "term compound tested in_ranking p bh_q finding"):
        term, compound = row["term"], row["compound"]
        require(term in TERMS and compound and (term, compound) not in seen,
                f"Unrecognized or repeated compound: {term} {compound}")
        require(all(row[field] in ("True", "False") for field in ("tested", "in_ranking", "finding")),
                f"Invalid flags: {term} {compound}")
        seen.add((term, compound))
        groups[term].append(row)

    for term, rows in groups.items():
        h1, h3 = panels[term], sizes[term]
        require(len(rows) == 116, f"Expected 116 compound and control records for {term}")
        ranked = [row for row in rows if row["in_ranking"] == "True"]
        require(len(ranked) == int(h1["n_tested_compounds"]) == 114, f"Wrong compound count: {term}")
        require(all(row["tested"] == "True" for row in ranked), f"Untested compound ranked: {term}")
        ranked.sort(key=lambda row: probability(row["p"]))
        p_values = [probability(row["p"]) for row in ranked]
        q_values = [probability(row["bh_q"]) for row in ranked]
        for observed, calculated in zip(q_values, bh_adjusted(p_values)):
            require(abs(observed - calculated) < Decimal("0.000000000002"),
                    f"BH correction does not reproduce: {term}")
        require(sum(p < Decimal(".05") for p in p_values) == int(h1["n_nominal_p_lt_0_05"]),
                f"Nominal count does not reproduce: {term}")
        require(sum(q < Decimal(".10") for q in q_values) == int(h1["n_bh_q_lt_0_10"]),
                f"Adjusted count does not reproduce: {term}")
        vehicle = next(row for row in rows if row["compound"] == "DMSO (pseudo-compound)")
        excluded = next(row for row in rows if row["compound"] == "AZD-8055")
        require(vehicle["in_ranking"] == "False" and vehicle["tested"] == "True" and
                vehicle["bh_q"] == "" and probability(vehicle["p"]) == probability(h1["vehicle_p"]),
                f"Unexpected vehicle record: {term}")
        require(excluded["in_ranking"] == excluded["tested"] == "False" and
                excluded["p"] == excluded["bh_q"] == "", f"Unexpected untested record: {term}")
        require(1 <= Decimal(h1["vehicle_rank_if_inserted"]) <= len(ranked) + 1,
                f"Unexpected vehicle rank: {term}")
        require(all(row["finding"] == "False" for row in rows), f"Unexpected finding: {term}")
        panel_p = probability(h1["panel_p"])
        require(int(h1["n_permutations"]) == 10_000 and panel_p >= Decimal(".025"),
                f"Unexpected panel result: {term}")
        require(int(h1["n_with_indicator"]) + int(h1["n_without_indicator"]) ==
                int(h1["n_lines"]) == int(h3["n_lines"]) == 41 and int(h1["n_lineages"]) == 11,
                f"Line counts disagree: {term}")
        interval = [Decimal(h3[key]) for key in ("lo95_pct", "lo90_pct", "estimate_pct", "hi90_pct", "hi95_pct")]
        require(all(value.is_finite() for value in interval) and sorted(interval) == interval,
                f"Invalid nested intervals: {term}")
        require(probability(h3["probability_interval_inside_margin"]) < Decimal(".8"),
                f"Unexpected equivalence diagnostic: {term}")

    controls = list(read_table("part_b_controls.tsv", "control response_scale compound n_lines p rank n_ranked passes"))
    require(len(controls) == 2 and {row["control"] for row in controls} == {"PC1", "K4"},
            "Unexpected Part B control rows")
    for row in controls:
        require(row["compound"] == "AMG-232" and Decimal(row["rank"]) == 1 and
                Decimal(row["n_ranked"]) == 114 and row["passes"] == "True" and
                probability(row["p"]) <= Decimal(".05"), "TP53 control record does not match")
    print("Part B: panel summaries, compound multiplicity correction, intervals and controls verified.")


def verify_phase6():
    labels = {}
    outcomes_by_drug = {}
    for row in read_table("phase6_comparison_labels.tsv.gz", "drug lineage tissue category skip_reason"):
        key = row["drug"], row["lineage"], row["tissue"]
        category = row["category"]
        require(all(key) and key not in labels and category in "ABCD" and len(category) == 1,
                f"Malformed or repeated Phase 6 comparison: {key}")
        require(bool(row["skip_reason"]) == (category == "D"), f"Inconsistent skip reason: {key}")
        labels[key] = category
        outcomes_by_drug.setdefault(row["drug"], set()).add(category)
    count = Counter(labels.values())
    saved = list(read_table("phase6_counts.tsv", "category count total"))
    require(len(saved) == 4 and all(int(row["total"]) == len(labels) for row in saved),
            "Wrong number of Phase 6 count rows")
    require({row["category"]: int(row["count"]) for row in saved} == dict(count)
            == {"A": 29, "B": 718, "C": 5375, "D": 33874}, "Phase 6 category counts do not match")
    require(len(labels) == 39_996 and len(outcomes_by_drug) == 50 and
            len({key[:2] for key in labels}) == 909 and len({key[2] for key in labels}) == 44,
            "Phase 6 comparison coverage does not match")
    require(sum(bool(categories & {"A", "B"}) for categories in outcomes_by_drug.values()) == 17 and
            sum("A" in categories for categories in outcomes_by_drug.values()) == 6,
            "Phase 6 compound counts do not match")
    transitions = Counter()
    changed_keys = set()
    path = HERE.parent / "Phase6_NoSampleFilter/results/replay_category_changes.csv"
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            key = row["drug"], row["lineage"], row["tissue"]
            before, after = row["archived_category"], row["replay_category"]
            require(key not in changed_keys and labels.get(key) == after, f"Invalid Phase 6 change: {key}")
            changed_keys.add(key)
            transitions[before, after] += 1
    require(transitions == {("C", "B"): 5, ("B", "C"): 4}, "Phase 6 archive-to-replay transitions")
    print("Earlier Phase 6 replay: 39,996 labels, nine archive differences, and compound counts verified.")


def main():
    try:
        verify_checksums()
        verify_part_a()
        verify_part_b()
        verify_phase6()
        verify_fellowship_phases()
    except (OSError, ValueError, KeyError, InvalidOperation, StopIteration) as error:
        print(f"Verification failed: {error}", file=sys.stderr)
        return 1
    print("All released result-table checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
