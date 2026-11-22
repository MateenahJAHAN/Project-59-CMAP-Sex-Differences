"""Independent checks from the full simulation records. Does not import the study summarizer."""

import math
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import beta, norm

R, B = 1000, 2000
TERMS = ("factY", "factX")
SCENARIOS = (
    "N0",
    "M10",
    "M20",
    "M30",
    "M50",
    "C1",
    "C6",
    "C23",
    "G1",
    "G6",
    "G23",
    "A1",
    "A2",
)
Z = norm.ppf(0.975)
Z_MEAN = 1.959964  # precision in the original mean-share interval
E3 = ("e3l3", "e3l5", "e3l3w", "e3l3s", "e3l5w", "e3l5s")


def require(test, message):
    if not test:
        raise ValueError(message)


def read(path):
    return pd.read_csv(path, sep="\t", float_precision="round_trip")


class Checks:
    def __init__(self):
        self.count = 0
        self.max_abs = 0.0
        self.fields = Counter()

    def same(self, field, actual, expected):
        if isinstance(expected, (str, bool, np.bool_)):
            require(actual == expected, f"{field}: {actual!r} versus {expected!r}")
        elif isinstance(expected, (int, np.integer)):
            require(float(actual) == expected, f"{field}: {actual} versus {expected}")
        else:
            a, b = float(actual), float(expected)
            if not (np.isnan(a) and np.isnan(b)):
                # Identical inputs can differ by a few ulps under independent arithmetic ordering.
                require(
                    np.isclose(a, b, rtol=8e-15, atol=8e-15),
                    f"{field}: {a:.17g} versus {b:.17g}",
                )
                self.max_abs = max(abs(a - b), self.max_abs)
        self.fields[field] += 1
        self.count += 1


def cp(x, n):
    return (
        0.0 if x == 0 else beta.ppf(0.025, x, n - x + 1),
        1.0 if x == n else beta.ppf(0.975, x + 1, n - x),
    )


def rate(values, shares=False):
    v = np.asarray(values, dtype=float)
    require(np.isfinite(v).all() and len(v) > 0, "missing or nonfinite numerator")
    n = len(v)
    p = v.mean()
    if shares:
        se = v.std(ddof=1) / math.sqrt(n)
        return p, max(0.0, p - Z_MEAN * se), min(1.0, p + Z_MEAN * se), se
    require(np.isin(v, [0.0, 1.0]).all(), "non-binary Wilson endpoint")
    se = math.sqrt(p * (1 - p) / n)
    denom = n + Z * Z
    centre = (p * n + Z * Z / 2) / denom
    half = Z * math.sqrt(n * p * (1.0 - p) + Z * Z / 4) / denom
    return p, max(0.0, centre - half), min(1.0, centre + half), se


def checked_grid(df, keys, levels):
    require(not df.duplicated(keys).any(), "duplicate simulation keys")
    observed = pd.MultiIndex.from_frame(df[keys])
    expected = pd.MultiIndex.from_product(levels, names=keys)
    require(
        expected.difference(observed).empty and observed.difference(expected).empty,
        "incomplete simulation grid",
    )


def inspect_outputs(root, caldir):
    runs = read(root / "raw/power_runs.tsv.gz")
    plants = read(root / "raw/power_planted.tsv.gz")
    diag = read(root / "raw/size_diagnostic_rows.tsv.gz")

    def result(name):
        return read(root / "results" / (name + ".tsv"))

    require(
        len(runs) == 26000 and runs.B.eq(B).all() and len(plants) == 120000,
        "incorrect simulation size",
    )
    checked_grid(runs, ["rep", "term", "scenario"], [range(R), TERMS, SCENARIOS])
    require(
        not plants.duplicated(["rep", "term", "scenario", "compound"]).any(),
        "duplicate planted record",
    )
    expected_plants = [s for s in SCENARIOS if s[0] in "CG"]
    sizes = plants.groupby(["rep", "term", "scenario"], as_index=False).size()
    checked_grid(sizes, ["rep", "term", "scenario"], [range(R), TERMS, expected_plants])
    require(
        sizes.apply(lambda s: s["size"] == int(s.scenario[1:]), axis=1).all(),
        "wrong number of planted drugs",
    )
    sl = {
        (t, s): g.sort_values("rep") for (t, s), g in runs.groupby(["term", "scenario"])
    }
    cur, checks = [], Checks()
    fig = result("figS3_power_values")
    require(
        len(fig) == 66 and fig.figure.eq("S3").all(), "unexpected number of S3 points"
    )
    series = (
        "E3 on level 3: 95% CI excludes 0",
        "E3 on level 5: 95% CI excludes 0",
        "panel-level test: p < 0.025",
    )
    sparse = (
        "planted compounds at q < 0.1 (share of planted)",
        "at least one planted compound at q < 0.1",
        series[2],
    )
    required = {
        (t, s, z)
        for t in TERMS
        for s in SCENARIOS
        if s[0] != "A"
        for z in (series if s[0] in "NM" else sparse)
    }
    require(
        set(zip(fig.term, fig.scenario, fig.series)) == required,
        "unexpected S3 estimand grid",
    )
    compound_hits = (
        plants.assign(hit=plants.q < 0.1).groupby(["term", "scenario", "rep"]).hit
    )
    planted_count = compound_hits.sum()
    for row in fig.itertuples():
        term, sc, endpoint = row.term, row.scenario, row.series
        group = sl[term, sc]
        shares = False
        family = "A" if sc[0] in "NM" else "B" if sc.startswith("C") else "C"
        require(
            row.panel.startswith(family + ":") and row.n_rep == R and row.B == B,
            "unexpected S3 metadata",
        )
        checks.same(
            "x",
            row.x,
            0.0 if sc == "N0" else int(sc[1:]) / (100 if family == "A" else 1),
        )
        if endpoint.startswith("E3"):
            est = "e3l5" if "level 5" in endpoint else "e3l3"
            lo = group[est + "_lo95"]
            hi = group[est + "_hi95"]
            good = lo.notna() & hi.notna()
            v = (lo[good] > 0) | (hi[good] < 0)
            require(
                np.array_equal(
                    v.to_numpy(),
                    group.loc[good, est + "_difference_declared"].astype(bool),
                ),
                "E3 decision inconsistent with interval",
            )
        elif endpoint == series[2]:
            v = group.panel_p.dropna() < 0.025
        else:
            hits = planted_count.loc[term, sc].reindex(group.rep)
            require(
                np.array_equal(hits, group.n_planted_q_lt_01)
                and group.n_planted.eq(int(sc[1:])).all(),
                "planted discoveries do not recount",
            )
            shares = endpoint.startswith("planted compounds")
            v = hits / int(sc[1:]) if shares else hits > 0
        if family == "A":
            # The released coordinate file supplies the Level 3 shift on the panel-test row as well.
            est = "e3l5" if "level 5" in endpoint else "e3l3"
            null = sl[term, "N0"]
            shift = np.mean(
                group[est + "_estimate"].to_numpy() - null[est + "_estimate"].to_numpy()
            )
            checks.same("planted_e3_pct", row.planted_e3_pct, shift)
        else:
            checks.same("planted_e3_pct", row.planted_e3_pct, np.nan)
        stats = rate(v, shares)
        for name, val in zip(("y", "y_lo95", "y_hi95", "y_mc_se"), stats):
            checks.same(name, getattr(row, name), val)
        cur.append(
            dict(
                term=term,
                panel=family,
                scenario=sc,
                series=endpoint,
                planned_assignments=R,
                usable_assignments=len(v),
                interval="replicate-share normal" if shares else "Wilson",
                eligible_drugs_min=int(group.n_tested.min()),
                eligible_drugs_max=int(group.n_tested.max()),
                analysis="null diagnostic" if sc == "N0" else "conditional detection",
            )
        )
    curves = pd.DataFrame(cur)
    coverage(result, sl, root, checks)
    calibration(result, sl, diag, caldir, checks)
    for row in result("e3_power_interpolated").itertuples():
        key = "level 3" if "level 3" in row.measure else "level 5"
        cc = fig[
            (fig.term == row.term)
            & fig.series.eq("E3 on " + key + ": 95% CI excludes 0")
        ].sort_values("planted_e3_pct")
        for h, nd in [
            ("lo", cc[cc.planted_e3_pct < 20].iloc[-1]),
            ("hi", cc[cc.planted_e3_pct >= 20].iloc[0]),
        ]:
            checks.same(
                "interpolated_scenario_" + h, getattr(row, "scenario_" + h), nd.scenario
            )
            checks.same("bracketing_rate", getattr(row, "rate_" + h), nd.y)
            checks.same(
                "bracketing_shift", getattr(row, "planted_" + h), nd.planted_e3_pct
            )
        checks.same(
            "interpolated_rate",
            row.interpolated_rate,
            np.interp(20.0, cc.planted_e3_pct, cc.y),
        )
    return curves, dict(
        scenario_rows=len(runs),
        planted_drug_rows=len(plants),
        diagnostic_rows=len(diag),
        checked_scalar_values=checks.count,
        checked_values_by_field=dict(sorted(checks.fields.items())),
        independent_formula_max_abs_difference=checks.max_abs,
    )


def coverage(result, sl, root, checks):
    for row in result("coverage_e3").itertuples():
        g, n0, e = sl[row.term, row.scenario], sl[row.term, "N0"], row.estimand
        ehat = g[e + "_estimate"].to_numpy()
        planted = (
            np.nanmean(ehat - n0[e + "_estimate"].to_numpy())
            if row.scenario != "N0"
            else 0.0
        )
        checks.same("E3_coverage_target", row.planted_value, planted)
        for lv in (90, 95):
            name = f"coverage{lv}"
            lo, hi = g[e + f"_lo{lv}"].to_numpy(), g[e + f"_hi{lv}"].to_numpy()
            good = np.isfinite(lo) & np.isfinite(hi)
            stats = rate((lo[good] <= planted) & (planted <= hi[good]))
            for suffix, val in zip(("", "_lo95", "_hi95", "_mc_se"), stats):
                checks.same("E3_" + name + suffix, getattr(row, name + suffix), val)
            checks.same(
                "E3_interval_halfwidth",
                getattr(row, "median_halfwidth" + str(lv)),
                np.nanmedian((hi - lo) / 2),
            )
    with np.load(root / "raw/power_e2.npz", allow_pickle=False) as raw:
        cube = raw["E2"]
        terms = list(raw["terms"])
        scen = list(raw["scenarios"])
        require(
            cube.shape == (1000, 2, 5, 115, 5)
            and terms == list(TERMS)
            and scen == list(SCENARIOS[:5])
            and list(raw["values"]) == ["est", "lo90", "hi90", "lo95", "hi95"],
            "unexpected E2 cube",
        )
        for row in result("coverage_e2").itertuples():
            ti, si = terms.index(row.term), scen.index(row.scenario)
            idx = (
                np.arange(115)
                if row.stratum == "all compounds"
                else raw["weak5"]
                if row.stratum.startswith("weak")
                else raw["strong5"]
            )
            checks.same("E2_stratum", row.n_compounds, len(idx))
            est = cube[:, ti, si][:, idx]
            valid = np.isfinite(est[:, :, 0])
            diff = est[:, :, 0] - cube[:, ti, 0][:, idx, 0]
            target = (
                np.divide(
                    np.nansum(diff, axis=0),
                    valid.sum(axis=0),
                    out=np.full(len(idx), np.nan),
                    where=valid.sum(axis=0) > 0,
                )
                if si
                else np.zeros(len(idx))
            )
            for lv, il, ih in [(90, 1, 2), (95, 3, 4)]:
                key = f"coverage{lv}"
                lower, upper = est[:, :, il], est[:, :, ih]
                good = np.isfinite(lower) & np.isfinite(upper)
                within = (lower <= target) & (upper >= target) & good
                count = good.sum(1)
                checks.same(
                    "E2_eligible_pairs",
                    getattr(row, "n_pairs" + str(lv)),
                    int(count.sum()),
                )
                rates = within.sum(1)[count > 0] / count[count > 0]
                se = rates.std(ddof=1) / np.sqrt(len(rates))
                p = within.sum() / count.sum()
                for suffix, val in [
                    ("", p),
                    ("_mc_se", se),
                    ("_lo95", p - Z_MEAN * se),
                    ("_hi95", p + Z_MEAN * se),
                ]:
                    checks.same("E2_" + key + suffix, getattr(row, key + suffix), val)


def calibration(result, sl, diag, caldir, checks):
    baseline = read(caldir / "G_calibration.tsv")
    baseline = baseline[baseline.run.eq("configured")]
    require(
        len(baseline) == 1200
        and not baseline.duplicated(["panel", "term", "pseudo"]).any(),
        "incomplete configured calibration",
    )
    variants = diag.variant.unique()
    require(
        len(variants) == 2 and len(diag) == 800 and diag.B.eq(B).all(),
        "incomplete paired diagnostic",
    )
    checked_grid(diag, ["rep", "term", "variant"], [range(200), TERMS, variants])
    for table, kind in [
        ("size_null_scenario", "N0"),
        ("size_panel_test", "configured"),
        ("size_diagnostic", "paired"),
    ]:
        for row in result(table).itertuples():
            null = sl[row.term, "N0"] if row.term in TERMS else None
            if kind == "N0":
                g = null
            elif kind == "configured":
                g = baseline[
                    (baseline.panel == row.panel) & (baseline.term == row.term)
                ].sort_values("pseudo")
            else:
                g = (
                    null.iloc[:200]
                    if row.variant.startswith("a:")
                    else diag[
                        (diag.term == row.term) & (diag.variant == row.variant)
                    ].sort_values("rep")
                )
                dif = null.panel_p.to_numpy()[:200] - g.panel_p.to_numpy()
                checks.same(
                    "paired_calibration_delta", row.paired_diff_vs_frozen, dif.mean()
                )
                checks.same(
                    "paired_calibration_se",
                    row.paired_diff_se,
                    dif.std(ddof=1) / np.sqrt(200),
                )
            p = g.panel_p.dropna().to_numpy()
            n = len(p)
            checks.same(
                "calibration_n", int(row.n_rep if kind == "N0" else row.n_pseudo), n
            )
            checks.same("mean_panel_p", row.mean_panel_p, p.mean())
            checks.same("median_panel_p", row.median_panel_p, np.median(p))
            for al, code in [(0.05, "005"), (0.025, "0025")]:
                x = int((p < al).sum())
                lo, hi = cp(x, n)
                prefix = "panel_" if kind == "N0" else ""
                checks.same(
                    "calibration_hits", getattr(row, "n_" + prefix + "p_lt_" + code), x
                )
                checks.same(
                    "calibration_rate",
                    getattr(row, "share_" + prefix + "p_lt_" + code),
                    x / n,
                )
                checks.same("Clopper_Pearson_lo", getattr(row, "cp95_lo_" + code), lo)
                checks.same("Clopper_Pearson_hi", getattr(row, "cp95_hi_" + code), hi)
                checks.same(
                    "calibration_reading",
                    getattr(row, "reading_" + code),
                    "conservative"
                    if hi < al
                    else "liberal"
                    if lo > al
                    else "consistent with nominal",
                )
                if kind == "configured":
                    checks.same(
                        "near_permutation_threshold",
                        getattr(row, "n_within_2se_perm_" + code),
                        int(np.sum(np.abs(p - al) <= 2 * math.sqrt(al * (1 - al) / B))),
                    )


def verify_packet(root, calibration, packet, allow_fixture=False):
    import hashlib
    import json

    names = {
        "figS3_power_values",
        "power_summary",
        "attenuation_e3",
        "coverage_e2",
        "coverage_e3",
        "e3_power_interpolated",
        "size_null_scenario",
        "size_diagnostic",
        "size_panel_test",
    }
    expected = {n + ".tsv" for n in names} | {
        "README.md",
        "comparison.json",
        "figure_denominators.tsv",
        "receipt.json",
    }
    hashes = {}
    for row in (packet / "SHA256SUMS").read_text().splitlines():
        checksum, filename = row.split("  ")
        require(
            filename in expected and filename not in hashes,
            "unexpected release filename",
        )
        with (packet / filename).open("rb") as stream:
            actual = hashlib.file_digest(stream, "sha256").hexdigest()
        require(checksum == actual, f"checksum failed: {filename}")
        hashes[filename] = actual
    require(
        set(hashes) == expected
        and {p.name for p in packet.iterdir()} == expected | {"SHA256SUMS"},
        "incomplete or unexpected packet files",
    )
    receipt = json.loads((packet / "receipt.json").read_text())
    if not allow_fixture:
        require(
            receipt.get("matched_checkpoints") == R
            and receipt["source_mode"].startswith("fresh study"),
            "requires a complete, fresh study",
        )
    curves, results = inspect_outputs(root, calibration)
    pd.testing.assert_frame_equal(
        read(packet / "figure_denominators.tsv"), curves, check_exact=True
    )
    require(
        results == receipt["independent_recounts"],
        "receipt disagrees with independent recount",
    )
    for name in names:
        src = read(root / "results" / (name + ".tsv"))
        if name == "size_panel_test":
            src = src.drop(columns="source_file")
        pd.testing.assert_frame_equal(
            src, read(packet / (name + ".tsv")), check_exact=True
        )
    return hashes


if __name__ == "__main__":
    import argparse
    from pathlib import Path

    p = argparse.ArgumentParser(
        description="Verify the published aggregate packet against the complete simulation output."
    )
    for arg in ["run", "calibration", "packet"]:
        p.add_argument("--" + arg, type=Path, required=True)
    p.add_argument("--allow-archive-fixture", action="store_true")
    a = p.parse_args()
    hashes = verify_packet(a.run, a.calibration, a.packet, a.allow_archive_fixture)
    print(f"{len(hashes)} public payloads verified against the source records.")
