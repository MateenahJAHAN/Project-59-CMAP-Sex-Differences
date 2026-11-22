"""Build a deterministic, aggregate-only S3 release from complete fresh outputs."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
import numpy as np
import pandas as pd
from audit import B, R, inspect_outputs, read, require

PUBLIC = (
    "figS3_power_values",
    "power_summary",
    "attenuation_e3",
    "coverage_e2",
    "coverage_e3",
    "e3_power_interpolated",
    "size_null_scenario",
    "size_diagnostic",
    "size_panel_test",
)
OTHER = (
    "detectable_effects",
    "e3_all_panels_scales",
    "e3_weak_templates",
    "misclass_e3",
    "misclass_flags",
    "tp53_compounds_q_lt_01",
    "tp53_yardstick",
    "harness_check5_power",
    "harness_checks",
    "harness_inputs_sha256",
)
RAW = ("power_runs", "power_planted", "size_diagnostic_rows")
PRIVATE_PATH_COLUMNS = {
    "size_panel_test": {"source_file"},
    "misclass_e3": {"flag_source"},
    "misclass_flags": {"flag_source"},
    "e3_all_panels_scales": {"source_file"},
    "harness_inputs_sha256": {"path"},
}


def sha(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def dump(p, obj):
    p.write_text(json.dumps(obj, sort_keys=True, indent=2, allow_nan=False) + "\n")


def numbers(a, b):
    a = np.asarray(a)
    b = np.asarray(b)
    if a.dtype.kind in "bifuc" and b.dtype.kind in "bifuc":
        nan = np.isnan(a) & np.isnan(b)
        changed = (a != b) & ~nan
        both = np.isfinite(a) & np.isfinite(b)
        delta = (
            float(np.max(np.abs(a[both].astype(float) - b[both].astype(float))))
            if both.any()
            else 0.0
        )
        return changed, delta
    sa = pd.Series(a.ravel(), dtype=object)
    sb = pd.Series(b.ravel(), dtype=object)
    changed = ~(sa.eq(sb) | (sa.isna() & sb.isna()))
    return changed.to_numpy().reshape(a.shape), None


def table_comparison(left, right, study_table):
    require(
        list(left) == list(right) and left.shape == right.shape,
        f"schema/row count mismatch: {study_table}",
    )
    any_changed = np.zeros(len(left), dtype=bool)
    changed_fields = {}
    numerical = 0
    for k in left:
        if k == "time_s" and study_table == "power_runs":
            continue
        a, b = left[k], right[k]
        if k in PRIVATE_PATH_COLUMNS.get(study_table, set()):
            # A relocated absolute filename describes the runtime only; retain the source basename.
            a, b = [
                c.map(lambda s: Path(s).name if isinstance(s, str) else s)
                for c in (a, b)
            ]
        changed, d = numbers(a, b)
        if d is not None:
            numerical += 1
        if changed.any():
            any_changed |= changed
            changed_fields[k] = {
                "entries": int(changed.sum()),
                "max_abs_finite_delta": d,
            }
    return dict(
        rows=len(left),
        scientific_columns=len(left.columns) - int(study_table == "power_runs"),
        numeric_columns=numerical,
        changed_rows=int(any_changed.sum()),
        changed_fields=changed_fields,
    )


def compare(root, archive):
    tables = {}
    for kind, files in [("results", PUBLIC + OTHER), ("raw", RAW)]:
        for fn in files:
            name = fn + (".tsv.gz" if kind == "raw" else ".tsv")
            tables[name] = table_comparison(
                read(root / kind / name), read(archive / kind / name), fn
            )
    arrays = {}
    with (
        np.load(root / "raw/power_e2.npz", allow_pickle=False) as a,
        np.load(archive / "raw/power_e2.npz", allow_pickle=False) as b,
    ):
        require(set(a.files) == set(b.files), "E2 fields differ")
        for k in sorted(a.files):
            x, y = a[k], b[k]
            require(
                x.shape == y.shape and x.dtype == y.dtype,
                f"E2 array schema differs: {k}",
            )
            dif, maximum = numbers(x, y)
            arrays[k] = {
                "shape": list(x.shape),
                "elements": int(x.size),
                "changed": int(dif.sum()),
                "maximum_finite_delta": maximum,
            }
    return {
        "tables": tables,
        "arrays": arrays,
        "all_scientific_fields_identical": not any(
            t["changed_rows"] for t in tables.values()
        )
        and not any(a["changed"] for a in arrays.values()),
    }


def checkpoints(root, cp, expected):

    require(
        cp.is_relative_to((root / "raw/checkpoints").resolve()),
        "checkpoint directory must belong to the fresh run",
    )
    meta = json.loads((cp / "run.json").read_text())
    require(
        (meta["R"], meta["B"], meta["tag"]) == (R, B, "power"), "wrong checkpoint run"
    )
    require(
        hashlib.sha256(json.dumps(meta, sort_keys=True).encode()).hexdigest()
        == cp.name,
        "checkpoint signature differs",
    )
    require(
        {f.name for f in cp.glob("rep_*.npz")}
        == {f"rep_{r:04d}.npz" for r in range(R)},
        "requires all 1000 completed assignments",
    )
    require(expected == sha(cp / "run.json"), "checkpoint manifest changed")
    unique, nruns, nplants, nelements = compare_saved(cp, root / "raw", range(R))
    return dict(
        checkpoint_set=cp.name,
        manifest_sha256=expected,
        matched_checkpoints=R,
        unique_assignment_fingerprints=unique,
        scenario_rows_from_checkpoints=nruns,
        planted_rows_from_checkpoints=nplants,
        interval_elements_from_checkpoints=nelements,
    )


def compare_saved(cp, source, ids):
    """Compare the supplied actual checkpoints with selected rows of an assembled output."""
    import io

    runs, plants = (
        read(source / (f + ".tsv.gz")) for f in ["power_runs", "power_planted"]
    )
    runs = runs[runs.rep.isin(ids)].reset_index(drop=True)
    plants = plants[plants.rep.isin(ids)].reset_index(drop=True)
    fresh_calls = []
    fresh_plants = []
    fingerprints = set()
    with np.load(source / "power_e2.npz", allow_pickle=False) as final:
        cube = final["E2"]
    for idx in ids:
        with np.load(cp / f"rep_{idx:04d}.npz", allow_pickle=False) as row:
            fp = str(row["fingerprint"].item())
            # Identical random assignments, if any, are legitimate.
            require(
                len(fp) == 64 and all(c in "0123456789abcdef" for c in fp),
                "malformed assignment fingerprint",
            )
            fingerprints.add(fp)
            records, planted = json.loads(row["records"].item())
            require(
                len(records) == 26
                and len(planted) == 120
                and all(r["rep"] == idx for r in records + planted),
                "incomplete checkpoint",
            )
            require(
                np.array_equal(row["intervals"], cube[idx], equal_nan=True),
                "saved E2 differs from the fresh assembled array",
            )
            fresh_calls.extend(records)
            fresh_plants.extend(planted)
    for name, records, tab in [
        ("power_runs", fresh_calls, runs),
        ("power_planted", fresh_plants, plants),
    ]:
        chunks = []
        for lo in range(0, len(records), 3000):
            columns = list(tab)
            # The study serializes dicts to TSV with pandas object/NA conventions.
            parts = pd.DataFrame(records[lo : lo + 3000]).reindex(columns=columns)
            chunks.append(
                pd.read_csv(
                    io.StringIO(parts.to_csv(sep="\t", index=False)),
                    sep="\t",
                    float_precision="round_trip",
                )
            )
        cmp = table_comparison(pd.concat(chunks, ignore_index=True), tab, name)
        require(
            cmp["changed_rows"] == 0,
            f"fresh output is inconsistent with saved {name} records: {cmp['changed_fields']}",
        )
    return (
        len(fingerprints),
        len(runs),
        len(plants),
        int(len(ids) * np.prod(cube.shape[1:])),
    )


def write_outputs(args, tmp):
    root, ref, cal = args.run, args.reference, args.calibration
    require(
        (root == ref) == args.archive_fixture,
        "supply a distinct fresh run; use --archive-fixture only for the archive",
    )
    if args.archive_fixture:
        require(
            "fixture" in args.out.name.lower(),
            "an archival test must target a directory explicitly named fixture",
        )
        provenance = {
            "source_mode": "archival fixture only; no fresh assignment verified"
        }
    else:
        require(
            args.checkpoints is not None,
            "final output requires the fresh checkpoint directory",
        )
        provenance = checkpoints(
            root, args.checkpoints, sha(args.checkpoints / "run.json")
        )
        provenance["source_mode"] = "fresh study with all 1000 checkpoints verified"
    for fn in ("harness_check5_power", "harness_checks"):
        harness = read(root / "results" / (fn + ".tsv"))
        require(harness.passes.eq(True).all(), f"{fn} reports a failed check")
    curves, findings = inspect_outputs(root, cal)
    comparisons = compare(root, ref)
    for f in PUBLIC:
        df = read(root / "results" / (f + ".tsv"))
        if f == "size_panel_test":
            df = df.drop(columns="source_file")
        require(
            not df.select_dtypes(exclude=["number", "bool"])
            .apply(
                lambda s: s.str.contains(
                    r"(?:/Users/|/private/|/home/|[A-Za-z]:\\)", na=False
                )
            )
            .any()
            .any(),
            "local path in exported table",
        )
        df.to_csv(tmp / (f + ".tsv"), sep="\t", index=False)
    curves.to_csv(tmp / "figure_denominators.tsv", sep="\t", index=False)
    dump(tmp / "comparison.json", comparisons)
    receipt = {
        **provenance,
        "planned_assignments": R,
        "permutations": B,
        "planned_figure_curves": 18,
        "planned_figure_points": 66,
        "complete_figure_curves": int(
            curves.assign(ok=curves.usable_assignments == R)
            .groupby(["panel", "series", "term"])
            .ok.all()
            .sum()
        ),
        "fully_observed_points": int((curves.usable_assignments == R).sum()),
        "figure_interval_method_counts": {
            k: int(v) for k, v in curves.interval.value_counts().items()
        },
        "independent_recounts": findings,
        "g_calibration_source_sha256": sha(cal / "G_calibration.tsv"),
        "reference_output_identity": {
            f: sha(ref / "results" / (f + ".tsv")) for f in PUBLIC
        },
        "verifier_sha256": {
            f: sha(Path(__file__).parent / f) for f in ["audit.py", "package.py"]
        },
    }
    dump(tmp / "receipt.json", receipt)
    ci = "| Term | Primary Level 3 detection at a 20% change | Null panel rejections, p < 0.025 | 95% exact interval |\n|---|---:|---:|---:|\n"
    inter = read(tmp / "e3_power_interpolated.tsv")
    size = read(tmp / "size_null_scenario.tsv")
    for term, label in [("factY", "Y marker"), ("factX", "XIST")]:
        row = size.loc[size.term == term].iloc[0]
        it = inter.loc[
            (inter.term == term) & inter.measure.str.contains("level 3")
        ].iloc[0]
        ci += f"| {label} | {it.interpolated_rate:.1%} | {row.n_panel_p_lt_0025}/{row.n_rep} | {row.cp95_lo_0025:.2%}–{row.cp95_hi_0025:.2%} |\n"
    heading = (
        "ARCHIVAL FIXTURE ONLY"
        if args.archive_fixture
        else "Figure S3 power and calibration"
    )
    text = f"""# {heading}

On the primary drug panel, effects were planted into measured expression profiles and Y/XIST labels reassigned within lineage. Each of {R:,} assignments used {B:,} permutations. These post hoc detection probabilities describe the 41-line design and the particular effect planted. They do not estimate the prevalence of biological effects.

{ci}
The 20% figures are interpolations between measured simulation points. The null counts measure the empirical rejection rate when the reassigned markers carry no effect. The [paired design diagnostic](size_diagnostic.tsv) uses the first 200 assignments. The [configured calibration](size_panel_test.tsv) summarizes a separate set of 1,200 assignment–contrast records.

[Figure S3](figS3_power_values.tsv) contains {receipt["fully_observed_points"]} of 66 coordinates with all assignments usable, across {receipt["complete_figure_curves"]} of 18 fully observed curves. The [denominators](figure_denominators.tsv) identify the usable assignments and eligible compound range for each coordinate. Sparse effects were planted only among eligible compounds. Fifty-four binary detection proportions use Wilson intervals. The twelve orange points average the fraction of planted compounds detected per assignment and use normal intervals from the variation among those fractions. The [coverage tables for pooled effects](coverage_e3.tsv) and [per-compound effects](coverage_e2.tsv) separately evaluate the 90% and 95% study intervals.

The [comparison](comparison.json) checks the raw simulation records, every E2 interval and all scientific fields of the revision's tables against the archived study. The [independent recount](receipt.json) derives figure coordinates, intervals, coverage and calibration from replicate records. Source tables and per-line expression values remain in the local study workspace.
"""
    (tmp / "README.md").write_text(text)
    files = sorted(p for p in tmp.iterdir() if p.is_file())
    (tmp / "SHA256SUMS").write_text("".join(f"{sha(p)}  {p.name}\n" for p in files))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ["run", "reference", "calibration", "out"]:
        p.add_argument("--" + key, required=True, type=Path)
    p.add_argument("--checkpoints", type=Path)
    p.add_argument("--archive-fixture", action="store_true")
    args = p.parse_args()
    for key in ["run", "reference", "calibration", "out", "checkpoints"]:
        value = getattr(args, key)
        if value is not None:
            setattr(args, key, value.resolve())
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=args.out.parent, prefix=".s3-build-") as temp:
        target = Path(temp) / "output"
        target.mkdir()
        write_outputs(args, target)
        if args.out.exists() and any(args.out.iterdir()):
            require(
                {p.name: sha(p) for p in target.iterdir()}
                == {p.name: sha(p) for p in args.out.iterdir()},
                "existing output differs; specify a new output directory",
            )
        else:
            if args.out.exists():
                args.out.rmdir()
            os.rename(target, args.out)
    print(
        "ARCHIVAL FIXTURE ONLY"
        if args.archive_fixture
        else "Complete fresh figure packet",
        args.out,
    )


if __name__ == "__main__":
    main()
