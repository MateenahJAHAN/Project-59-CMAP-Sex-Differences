"""Recount the Tahoe donor-label comparison from the public 50-line expression table."""
from collections import Counter
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
with (ROOT / "marker_transfer/reference_cutoffs.tsv").open() as source:
    cutoffs = {r["score"]: float(r["cutoff"]) for r in csv.DictReader(source, delimiter="\t")}
with (ROOT.parent / "data/tahoe-emeraldbay-expression-states.tsv").open() as source:
    rows = [r for r in csv.DictReader(source, delimiter="\t") if r["source"] == "Tahoe-100M"]
assert len(rows) == len({r["cellosaurus_id"] for r in rows}) == 50
counts = Counter()
for row in rows:
    y = float(row["y_score"]) >= cutoffs["y_score"]
    xist = float(row["xist"]) >= cutoffs["XIST"]
    counts[row["sex_label"], y, xist] += 1
    if row["cell_name"] == "A498":
        assert row["sex_label"] == "Male" and xist and not y
for sex in ("Female", "Male", "Unknown"):
    print(sex, {"Y only": counts[sex, True, False], "XIST only": counts[sex, False, True],
                "both low": counts[sex, False, False], "both high": counts[sex, True, True]})
assert counts["Male", False, False] == 14
assert counts["Male", False, False] + counts["Male", False, True] == 15
assert counts["Female", False, False] == 12
