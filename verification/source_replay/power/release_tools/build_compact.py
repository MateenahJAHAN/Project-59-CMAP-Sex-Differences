"""Export a public, standard-library replay of Figure S3 and Note S17."""

import argparse
import csv
import gzip
import io
import json
import math
import os
from pathlib import Path
import shutil
import tempfile
from package import checkpoints
from replay_s3 import (
    B,
    R,
    POOLED,
    SPARSE,
    FIG_FIELDS,
    RUN_FIELDS,
    SINGLE_FIELDS,
    require,
    digest,
    read_packet,
    verify,
)


def rows(path):
    with gzip.open(path, "rt", newline="", encoding="utf8") as src:
        yield from csv.DictReader(src, delimiter="\t")


def tsv_gz(path, fields, items):
    with (
        path.open("wb") as output,
        gzip.GzipFile(
            filename="", mode="wb", fileobj=output, mtime=0, compresslevel=9
        ) as compressed,
        io.TextIOWrapper(compressed, encoding="utf8", newline="") as text,
    ):
        writer = csv.DictWriter(
            text, fieldnames=fields, lineterminator="\n", delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(items)


def number(value):
    if not value:
        return ""
    v = float(value)
    require(math.isfinite(v) and v.is_integer(), "invalid integer field")
    return int(v)


def assignments(root):
    for r in rows(root / "raw/power_runs.tsv.gz"):
        sc = r["scenario"]
        if sc not in POOLED + SPARSE:
            continue
        p = float(r["panel_p"]) if r["panel_p"] else None
        b = round(p * (B + 1)) - 1 if p is not None else ""
        require(
            int(r["B"]) == B and (p is None or abs(p - (b + 1) / (B + 1)) < 2e-16),
            "panel p not on permutation grid",
        )
        out = dict(
            rep=int(r["rep"]),
            term=r["term"],
            scenario=sc,
            tested_drugs=number(r["n_tested"]),
            panel_exceedances=b,
            planted_drugs=number(r["n_planted"]) if sc in SPARSE else "",
            planted_hits=number(r["n_planted_q_lt_01"]) if sc in SPARSE else "",
        )
        for scale in ("l3", "l5"):
            for short, field in [
                ("est", "estimate"),
                ("lo95", "lo95"),
                ("hi95", "hi95"),
            ]:
                out[scale + "_" + short] = (
                    r["e3" + scale + "_" + field] if sc in POOLED else ""
                )
        yield out


def single_plants(root):
    others = {
        (r["rep"], r["term"], r["scenario"]): number(r["n_other_q_lt_01"])
        for r in rows(root / "raw/power_runs.tsv.gz")
        if r["scenario"] in ("C1", "G1")
    }
    for r in rows(root / "raw/power_planted.tsv.gz"):
        if r["scenario"] not in ("C1", "G1"):
            continue
        yield {k: r[k] for k in ("rep", "term", "scenario", "compound", "p", "q")} | {
            "bh_hit": int(float(r["q"]) < 0.1),
            "other_bh_hits": others[(r["rep"], r["term"], r["scenario"])],
        }


def prepare(args, out):
    final = args.final
    if args.archive_fixture:
        require(
            final is None
            and args.checkpoints is None
            and "fixture" in args.out.name.lower(),
            "archive fixture requires a separate, labelled output",
        )
        prov = dict(
            mode="ARCHIVAL FIXTURE ONLY", complete_fresh_assignment_receipt=False
        )
        reference = args.run / "results/figS3_power_values.tsv"
    else:
        require(
            final is not None and args.checkpoints is not None,
            "complete release and matching checkpoints required",
        )
        checks = {
            name: h
            for h, name in (
                line.split("  ")
                for line in (final / "SHA256SUMS").read_text().splitlines()
            )
        }
        reference = final / "figS3_power_values.tsv"
        receipt = final / "receipt.json"
        require(
            all(
                digest(final / f) == checks[f]
                for f in ["receipt.json", "figS3_power_values.tsv"]
            ),
            "invalid fresh packet certificate",
        )
        rec = json.loads(receipt.read_text())
        cp = args.checkpoints.resolve()
        require(
            rec["source_mode"] == "fresh study with all 1000 checkpoints verified"
            and rec["matched_checkpoints"] == R
            and rec["planned_assignments"] == R
            and rec["permutations"] == B
            and rec["checkpoint_set"] == cp.name,
            "requires 1000 matched fresh assignments",
        )
        verified = checkpoints(args.run.resolve(), cp, rec["manifest_sha256"])
        require(
            all(rec[k] == value for k, value in verified.items()),
            "records differ from certified fresh packet",
        )
        prov = dict(
            mode="fresh complete simulation",
            complete_fresh_assignment_receipt=True,
            matched_checkpoints=R,
            fresh_receipt_sha256=digest(receipt),
            checkpoint_manifest_sha256=rec["manifest_sha256"],
        )
    tsv_gz(out / "assignments.tsv.gz", RUN_FIELDS, assignments(args.run))
    tsv_gz(out / "single_plants.tsv.gz", SINGLE_FIELDS, single_plants(args.run))
    with reference.open(newline="", encoding="utf8") as stream:
        ref = csv.DictReader(stream, delimiter="\t")
        require(ref.fieldnames == FIG_FIELDS, "unexpected figure reference schema")
        tsv_gz(out / "reference_coordinates.tsv.gz", FIG_FIELDS, ref)
    shutil.copyfile(Path(__file__).with_name("replay_s3.py"), out / "replay_s3.py")
    _, _, report = read_packet(out)
    (out / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    prov["reference_coordinate_sha256"] = digest(reference)
    (out / "provenance.json").write_text(json.dumps(prov, indent=2) + "\n")
    observed = "\n".join(
        f"| {e['term']} | {e['scenario']} | {e['rep']} | {e['compound']} | {e['permutation_exceedances']} | {e['p']:.6f} | {e['q']:.6f} | {e['other_compounds_q_below_0p1']} |"
        for e in report["positive_exceedance_q_hits"]
    )
    head = (
        "ARCHIVAL FIXTURE ONLY"
        if args.archive_fixture
        else "Replay Figure S3 and check Note S17"
    )
    hits = sum(s["hits"] for s in report["single_plant_counts"])
    n = len(report["positive_exceedance_q_hits"])
    evidence = (
        f"Of {hits:,} planted compounds with q < 0.1, {n} had at least one permutation statistic reach the observed statistic. "
        f"The other {hits - n} had zero exceedances."
    )
    text = f"""# {head}

Use Python 3.11 or later to reconstruct all 66 coordinates and their Monte Carlo intervals. This command runs offline with the standard library. The optional output directory receives the recalculated coordinates, the number of usable assignments for every point and the individual cases behind Note S17.

```sh
python3 replay_s3.py --out replay
```

The packet contains 22,000 assignment records across the two terms and eleven plotted scenarios. The input for the panel test is the permutation exceedance count b. Both terms use p = (b + 1) / 2,001 and p < 0.025, equivalent to b ≤ 49. The pooled endpoint compares each 95% interval with zero. Its paired shift is the planted estimate minus the null estimate from the same assignment. Detection of a planted compound requires BH q < 0.1.

The orange series averages the proportion of planted compounds detected *within each assignment*. The green series counts assignments with one or more detected planted compounds. The orange intervals use the standard error of those assignment-level fractions. Other intervals use Wilson's formula for event counts. A missing measurement is excluded from that endpoint's denominator. The null has no planted compounds, so a planted-hit share is undefined. Null panel rejections and pooled interval exclusions have their own valid denominators.

## Single-compound simulations

The file `single_plants.tsv.gz` contains all 4,000 single-plant cases: the 50% response increase (C1) and the 49-gene shift (G1) for Y and XIST. {evidence} The BH step-up calculation adjusts p values jointly across the eligible drugs in an assignment. A single *planted* compound can pass with positive exceedances when other compounds also pass. These paired records establish what happened at 2,000 permutations. A rate at 10,000 permutations would require running that comparison.

| Term | Scenario | Assignment | Planted compound | Exceedances | p | BH q | Other rejections |
|---|---|---:|---|---:|---:|---:|---:|
{observed}

`assignments.tsv.gz` contains the tested-drug count, panel exceedances, sparse-scenario planted counts and hits, and the pooled estimates and intervals used by the figure. Blank fields indicate an inapplicable or unavailable endpoint. `reference_coordinates.tsv.gz` supplies the figure values used to check the recalculation. `verification.json` records the full-precision comparisons and denominators. `provenance.json` identifies the source of this packet.
"""
    (out / "README.md").write_text(text)
    (out / "SHA256SUMS").write_text(
        "".join(f"{digest(f)}  {f.name}\n" for f in sorted(out.iterdir()))
    )
    verify(out)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ["run", "out", "final", "checkpoints"]:
        p.add_argument("--" + arg, type=Path, required=arg in ("run", "out"))
    p.add_argument("--archive-fixture", action="store_true")
    a = p.parse_args()
    a.out = a.out.resolve()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=a.out.parent, prefix=".compact-s3-") as tmp:
        packet = Path(tmp) / "packet"
        packet.mkdir()
        prepare(a, packet)
        if a.out.exists() and any(a.out.iterdir()):
            require(
                {f.name: digest(f) for f in a.out.iterdir()}
                == {f.name: digest(f) for f in packet.iterdir()},
                "existing packet differs; select a new output directory",
            )
        else:
            if a.out.exists():
                a.out.rmdir()
            os.rename(packet, a.out)
    print(f"{sum(f.stat().st_size for f in a.out.iterdir()):,} bytes in {a.out}")


if __name__ == "__main__":
    main()
