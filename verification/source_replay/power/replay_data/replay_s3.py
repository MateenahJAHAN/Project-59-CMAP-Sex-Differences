"""Recalculate Figure S3 from the public packet with Python's standard library."""

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
from statistics import fmean, stdev

R, B = 1000, 2000
TERMS = ("factY", "factX")
POOLED = ("N0", "M10", "M20", "M30", "M50")
SPARSE = ("C1", "C6", "C23", "G1", "G6", "G23")
PANEL = {
    "A": "A: pooled magnitude",
    "B": "B: 50% larger response in k compounds",
    "C": "C: +1 SD in a 49-gene module in k compounds",
}
E3 = ["E3 on level 3: 95% CI excludes 0", "E3 on level 5: 95% CI excludes 0"]
H1 = "panel-level test: p < 0.025"
SHARE = "planted compounds at q < 0.1 (share of planted)"
ANY = "at least one planted compound at q < 0.1"
RUN_FIELDS = [
    "rep",
    "term",
    "scenario",
    "tested_drugs",
    "panel_exceedances",
    "planted_drugs",
    "planted_hits",
    "l3_est",
    "l3_lo95",
    "l3_hi95",
    "l5_est",
    "l5_lo95",
    "l5_hi95",
]
SINGLE_FIELDS = [
    "rep",
    "term",
    "scenario",
    "compound",
    "p",
    "q",
    "bh_hit",
    "other_bh_hits",
]
FIG_FIELDS = [
    "figure",
    "panel",
    "series",
    "term",
    "scenario",
    "x",
    "x_label",
    "planted_e3_pct",
    "y",
    "y_lo95",
    "y_hi95",
    "y_mc_se",
    "n_rep",
    "B",
]
TEXT = {"figure", "panel", "series", "term", "scenario", "x_label"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(file):
    with file.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def read_rows(file, fields):
    with gzip.open(file, "rt", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        require(reader.fieldnames == fields, f"unexpected columns in {file.name}")
        return list(reader)


def as_num(string):
    if string == "":
        return None
    value = float(string)
    require(math.isfinite(value), "nonfinite field must be encoded as empty")
    return value


def mean_ci(values, share):
    if not values:
        return (None,) * 4
    p = fmean(values)
    n = len(values)
    if share:
        if n < 2:
            return p, None, None, None
        se = stdev(values) / math.sqrt(n)
        half = 1.959964 * se
        return p, max(0.0, p - half), min(1.0, p + half), se
    require(
        all(x in (0, 1) for x in values), "Wilson interval requires binary observations"
    )
    z = 1.959963984540054
    se = math.sqrt(p * (1 - p) / n)
    fac = n + z * z
    c = (p * n + 0.5 * z * z) / fac
    half = z * math.sqrt(n * p * (1 - p) + z * z / 4) / fac
    return p, max(0.0, c - half), min(1.0, c + half), se


def excludes_zero(row, prefix):
    lower = as_num(row[prefix + "_lo95"])
    upper = as_num(row[prefix + "_hi95"])
    require((lower is None) == (upper is None), "incomplete CI")
    if lower is None:
        return None
    require(lower <= upper, "reversed CI")
    return lower > 0 or upper < 0


def grid(rows, scenarios):
    indexed = {(r["term"], r["scenario"], int(r["rep"])): r for r in rows}
    require(
        len(rows) == len(indexed) == R * len(TERMS) * len(scenarios),
        "wrong or duplicate assignment keys",
    )
    require(
        set(indexed)
        == {(term, sc, i) for term in TERMS for sc in scenarios for i in range(R)},
        "missing assignment keys",
    )
    return indexed


def coordinates(rows):
    data = grid(rows, POOLED + SPARSE)
    results = []
    denoms = []
    for term in TERMS:
        for sc in POOLED + SPARSE:
            label = "A" if sc in POOLED else "B" if sc.startswith("C") else "C"
            x = 0.0 if sc == "N0" else int(sc[1:]) / (100 if label == "A" else 1)
            for metric in E3 + [H1] if sc in POOLED else [SHARE, H1, ANY]:
                prefix = "l5" if metric == E3[1] else "l3"
                vals = []
                shifts = []
                for i in range(R):
                    row = data[term, sc, i]
                    if metric in E3:
                        val = excludes_zero(row, prefix)
                    elif metric == H1:
                        exceeds = row["panel_exceedances"]
                        if exceeds == "":
                            val = None
                        else:
                            b = int(exceeds)
                            require(0 <= b <= B, "invalid panel permutation count")
                            val = (1 + b) / (B + 1) < 0.025
                    else:
                        n = int(row["planted_drugs"])
                        hits = int(row["planted_hits"])
                        require(0 <= hits <= n, "impossible planted-hit count")
                        require(
                            n == int(sc[1:]),
                            "sparse scenario differs from its planned denominator",
                        )
                        val = (
                            (hits / n if metric == SHARE else hits > 0)
                            if n > 0
                            else None
                        )
                    if val is not None:
                        vals.append(val)
                    if sc in POOLED:
                        est = as_num(row[prefix + "_est"])
                        null = as_num(data[term, "N0", i][prefix + "_est"])
                        if est is not None and null is not None:
                            shifts.append(est - null)
                y, lo, hi, se = mean_ci(vals, metric == SHARE)
                results.append(
                    dict(
                        figure="S3",
                        panel=PANEL[label],
                        series=metric,
                        term=term,
                        scenario=sc,
                        x=x,
                        x_label="planted pooled response-size difference (fraction)"
                        if label == "A"
                        else "number of compounds with the planted effect (k)",
                        planted_e3_pct=fmean(shifts) if shifts else None,
                        y=y,
                        y_lo95=lo,
                        y_hi95=hi,
                        y_mc_se=se,
                        n_rep=R,
                        B=B,
                    )
                )
                denoms.append(
                    dict(
                        term=term,
                        scenario=sc,
                        metric=metric,
                        planned=R,
                        usable=len(vals),
                        interval="replicate-share normal"
                        if metric == SHARE
                        else "Wilson",
                    )
                )
    return data, results, denoms


def single_plants(rows, data):
    records = grid(rows, ("C1", "G1"))
    by = {}
    examples = []
    for (term, sc, i), row in sorted(records.items()):
        ar = data[term, sc, i]
        m = int(ar["tested_drugs"])
        p = float(row["p"])
        q = float(row["q"])
        require(
            1 <= m <= 115
            and 0 < p <= 1
            and p - 1e-14 <= q <= 1
            and row["compound"] != "",
            "invalid single-plant statistics",
        )
        k = round(p * (B + 1)) - 1
        require(
            0 <= k <= B and abs(p - (1 + k) / (1 + B)) < 2e-16,
            "p inconsistent with permutation grid",
        )
        hit = q < 0.1
        others = int(row["other_bh_hits"])
        require(
            row["bh_hit"] in {"0", "1"}
            and hit == bool(int(row["bh_hit"]))
            and int(ar["planted_drugs"]) == 1
            and hit == int(ar["planted_hits"])
            and 0 <= others < m,
            "S17 case inconsistent",
        )
        stats = by.setdefault(
            (term, sc),
            dict(cases=0, hits=0, hits_zero_exceedances=0, hits_positive_exceedances=0),
        )
        stats["cases"] += 1
        stats["hits"] += int(hit)
        if hit:
            stats[
                "hits_zero_exceedances" if k == 0 else "hits_positive_exceedances"
            ] += 1
        if hit and k > 0:
            # Any BH hit at rank j requires p <= q_threshold*j/m at some step-up cutoff.
            require(
                p < 0.1 * (others + 1) / m,
                "q hit inconsistent with total rejection count",
            )
            examples.append(
                dict(
                    term=term,
                    scenario=sc,
                    rep=i,
                    compound=row["compound"],
                    permutation_exceedances=k,
                    p=p,
                    q=q,
                    tested_drugs=m,
                    other_compounds_q_below_0p1=others,
                )
            )
    return [dict(term=t, scenario=s, **v) for (t, s), v in sorted(by.items())], examples


def compare(results, ref):
    def keys(rows):
        return {(r["term"], r["scenario"], r["series"]): r for r in rows}

    actual, expected = keys(results), keys(ref)
    require(
        len(ref) == len(actual) == 66 and actual.keys() == expected.keys(),
        "S3 reference grid differs",
    )
    maxerr = 0.0
    for key, row in actual.items():
        for col in FIG_FIELDS:
            other = expected[key][col]
            val = row[col]
            if col in TEXT:
                require(other == val, f"S3 label mismatch: {key} {col}")
            else:
                v = as_num(other)
                if v is None or val is None:
                    require(v == val, f"S3 missingness mismatch: {key} {col}")
                else:
                    require(
                        math.isclose(val, v, rel_tol=8e-15, abs_tol=8e-15),
                        f"S3 number mismatch: {key} {col}",
                    )
                    maxerr = max(abs(val - v), maxerr)
    return maxerr


def write_tsv(file, fields, rows):
    with file.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def read_packet(packet):
    rows = read_rows(packet / "assignments.tsv.gz", RUN_FIELDS)
    for r in rows:
        if r["scenario"] in POOLED:
            require(
                r["planted_drugs"] == r["planted_hits"] == "",
                "pooled or null assignment has a sparse endpoint",
            )
        else:
            require(
                all(r[k] == "" for k in RUN_FIELDS[7:]),
                "sparse assignment includes unused pooled fields",
            )
    data, coords, den = coordinates(rows)
    single = read_rows(packet / "single_plants.tsv.gz", SINGLE_FIELDS)
    counts, cases = single_plants(single, data)
    delta = compare(
        coords, read_rows(packet / "reference_coordinates.tsv.gz", FIG_FIELDS)
    )
    cnt = [
        [
            t,
            sum(
                r["usable"]
                for r in den
                if r["term"] == t and r["metric"] == H1 and r["scenario"] == "N0"
            ),
            sum(
                (1 + int(row["panel_exceedances"])) / (B + 1) < 0.025
                for (term, sc, _), row in data.items()
                if term == t and sc == "N0" and row["panel_exceedances"] != ""
            ),
        ]
        for t in TERMS
    ]
    receipt = dict(
        figure_points=len(coords),
        planned_curves=18,
        points_with_all_assignments=sum(r["usable"] == R for r in den),
        null_panel_diagnostics=[
            dict(term=t, usable_assignments=n, rejections=h) for t, n, h in cnt
        ],
        maximum_absolute_coordinate_difference=delta,
        scenario_assignment_rows=len(rows),
        single_plant_rows=len(single),
        single_plant_counts=counts,
        positive_exceedance_q_hits=cases,
    )
    return coords, den, receipt


def verify(packet, out=None):
    public_files = {
        "README.md",
        "assignments.tsv.gz",
        "single_plants.tsv.gz",
        "reference_coordinates.tsv.gz",
        "replay_s3.py",
        "verification.json",
        "provenance.json",
    }
    checks = {}
    for line in (packet / "SHA256SUMS").read_text().splitlines():
        h, name = line.split("  ")
        require(
            name in public_files and name not in checks, "unexpected checksum entry"
        )
        require(digest(packet / name) == h, f"checksum failure in {name}")
        checks[name] = h
    require(set(checks) == public_files, "missing checksum entry")
    coords, den, receipt = read_packet(packet)
    saved = json.loads((packet / "verification.json").read_text())
    require(saved.keys() == receipt.keys(), "unexpected verification summary")
    # Floating summation (Python 3.11 vs 3.12+) can vary slightly.
    for k in saved:
        if k == "maximum_absolute_coordinate_difference":
            continue
        require(saved[k] == receipt[k], f"independent replay changed {k}")
    if out:
        out.mkdir(exist_ok=True, parents=True)
        write_tsv(out / "figS3_recalculated.tsv", FIG_FIELDS, coords)
        write_tsv(out / "figS3_denominators.tsv", list(den[0]), den)
        (out / "s17_checks.json").write_text(json.dumps(receipt, indent=2) + "\n")
    return receipt


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--packet", type=Path, default=Path(__file__).resolve().parent)
    p.add_argument(
        "--out",
        type=Path,
        help="Optional directory for recomputed coordinates, denominators and S17 cases.",
    )
    a = p.parse_args()
    print(json.dumps(verify(a.packet, a.out), indent=2))
