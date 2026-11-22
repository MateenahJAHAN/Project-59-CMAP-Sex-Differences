"""Audit the power calculation's input matrices and random designs. No hypothesis tests are run."""

import argparse
import ast
import hashlib
import json
import os
import sys
import types
from pathlib import Path

for name in [
    "OPENBLAS_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
]:
    os.environ[name] = "1"
import numpy as np
import pandas as pd


def digest(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def comparison(ref, obs, full_digest=False):
    a, b = np.asarray(ref), np.asarray(obs)
    if a.shape != b.shape:
        raise ValueError(f"Shape differs: {a.shape} versus {b.shape}")
    if a.dtype.kind in "USO" or b.dtype.kind in "USO":
        diff = int(np.sum(a != b))
        rec = dict(shape=list(a.shape), changed=diff)
        if full_digest:
            rec["ordered_sha256"] = hashlib.sha256(
                "\n".join(map(str, a.ravel())).encode()
            ).hexdigest()
    else:
        diff = missing = present = 0
        maximum = 0.0
        for x, y in zip(
            np.array_split(a, min(len(a), 32)), np.array_split(b, min(len(b), 32))
        ):
            fx, fy = np.isfinite(x), np.isfinite(y)
            missing += int((fx != fy).sum())
            finite = fx & fy
            present += int(finite.sum())
            if finite.any():
                diff += int((x[finite] != y[finite]).sum())
                maximum = max(
                    maximum,
                    float(
                        np.max(
                            np.abs(x[finite].astype(float) - y[finite].astype(float))
                        )
                    ),
                )
        rec = dict(
            shape=list(a.shape),
            finite_entries=present,
            changed_finite=diff,
            changed_availability=missing,
            max_abs_delta=maximum,
        )
        if full_digest:
            rec["ordered_sha256"] = hashlib.sha256(
                np.ascontiguousarray(a).tobytes()
            ).hexdigest()
        diff += missing
    if diff:
        raise ValueError(f"Input mismatch: {rec}")
    return rec


def load_power(path, engine, input_root):
    # Preserve the released scientific code. Relocate only its two installation paths.
    module = ast.parse(path.read_text())
    assignments = {"REPO": engine, "RUN_DIR": input_root}
    for node in module.body:
        if (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
        ):
            key = node.targets[0].id
            if key in assignments:
                node.value = ast.copy_location(
                    ast.Constant(str(assignments.pop(key))), node.value
                )
    assert not assignments
    ctx = {"__file__": str(path), "__name__": "power_audit_common"}
    exec(compile(module, str(path), "exec"), ctx)
    ctx["OVERRIDES"] = {
        k: str(input_root / rel)
        for k, rel in [
            ("inputs", "inputs"),
            ("tensors_D", "tensors/D"),
            ("tensors_E", "tensors/E"),
        ]
    }
    return types.SimpleNamespace(**ctx)


def source_context(K, root, old):
    p = root / "part_b/run"
    dd, ee = p / "D", p / "E"
    labs = p / "D_labels/results/D"
    dten = dd / "data/D/tensors"
    axes = dd / "results/D/tensors"
    eroot = ee / "data/E/tensors"
    lines = pd.read_csv(axes / "P1_lines.tsv", sep="\t").cell_iname.tolist()
    comps = pd.read_csv(axes / "P1_compounds.tsv", sep="\t").cmap_name.tolist()
    gene = pd.read_csv(dd / "results/D/D_landmarks.tsv", sep="\t").gene_symbol.tolist()
    panel = K.C.Panel(
        "P1",
        dten / "P1_resp.npy",
        axes / "P1_lines.tsv",
        axes / "P1_compounds.tsv",
        gene,
        K.C.read_covariates(str(labs / "D_line_covariates.tsv")),
        old.cfg.run("key_seed_base") + old.cfg.const("panel_numbers")["P1"],
        old.cfg.run("B_full"),
        dten / "P1_DMSO_resp.npy",
    )
    el = (eroot / "P1_d111_lines.txt").read_text().splitlines()
    ec = (eroot / "P1_d111_compounds.txt").read_text().splitlines()
    eg = pd.read_csv(ee / "results/E/E_landmarks.tsv", sep="\t").gene_symbol.tolist()
    assert all(
        len(a) == len(set(a)) and set(a) == set(b)
        for a, b in [(el, lines), (ec, comps), (eg, gene)]
    )
    indexes = [
        [orig.index(v) for v in dst]
        for orig, dst in [(el, lines), (ec, comps), (eg, gene)]
    ]
    r3 = np.load(eroot / "R3_P1_d111.npy", mmap_mode="r")[np.ix_(*indexes)]
    ctx = types.SimpleNamespace(
        P=panel,
        R3=r3,
        avail5=panel.avail,
        avail3=np.isfinite(r3[:, :, 0]),
        cfg=old.cfg,
        tpX=K.GS.tp53_columns(panel.cov, K.GS.tp53_code(old.cfg))[0],
        rule=K.GS.rule_for(old.cfg),
    )
    return ctx, indexes


def well_reconstruction(root, D, indexes):
    e = root / "part_b/run/E"
    saved = e / "data/E/tensors"
    lm = lambda s: np.load(saved / s, mmap_mode="r")
    xt, vc, stored = map(lm, ["l3_trt_P1.npy", "l3_ctl_P1.npy", "r3_wells_P1.npy"])
    tmeta = pd.read_csv(saved / "l3_trt_P1_rows.tsv", sep="\t")
    cmeta = pd.read_csv(saved / "l3_ctl_P1_rows.tsv", sep="\t")
    wells = pd.read_csv(e / "results/E/E_wells.tsv.gz", sep="\t", low_memory=False)

    def align(rows, role):
        selected = wells[(wells.panel == "P1") & (wells.role == role)].set_index(
            "sample_id"
        )
        assert selected.index.is_unique and rows.sample_id.is_unique
        assert set(rows.sample_id) == set(selected.index)
        out = selected.loc[rows.sample_id].reset_index()
        assert np.array_equal(out.row_index, rows.row_index)
        return out

    wt, wc = align(tmeta, "treated"), align(cmeta, "dmso")
    groups = wc.groupby("det_plate").indices
    permitted = {
        plate for plate in set(wt.det_plate) if len(groups.get(plate, ())) >= 2
    }
    mask = wt.det_plate.isin(permitted)
    wm = wt[mask].reset_index(drop=True)
    stored_meta = pd.read_csv(saved / "r3_wells_P1_rows.tsv", sep="\t")
    assert wm[stored_meta.columns].equals(stored_meta)
    plate_pos = {p: i for i, p in enumerate(sorted(permitted))}
    means = np.stack(
        [np.asarray(vc[groups[p]], dtype=float).mean(axis=0) for p in sorted(permitted)]
    )
    plate_ids = wm.det_plate.map(plate_pos).to_numpy()
    selected = np.flatnonzero(mask.to_numpy())
    selected_dose = wm.pert_idose.eq("1.11 uM").to_numpy()
    lines = (saved / "P1_d111_lines.txt").read_text().splitlines()
    comps = (saved / "P1_d111_compounds.txt").read_text().splitlines()
    unit = (
        wm.cell_iname.map({s: i for i, s in enumerate(lines)}).to_numpy() * len(comps)
        + wm.cmap_name.map({s: i for i, s in enumerate(comps)}).to_numpy()
    )
    sums = np.zeros((len(lines) * len(comps), len(D.P.genes)))
    n = np.zeros(len(sums), dtype="int64")
    maxdelta = 0.0
    for lo in range(0, len(selected), 128):
        ix = selected[lo : lo + 128]
        size = len(ix)
        computed = np.asarray(xt[ix], dtype=float) - means[plate_ids[lo : lo + size]]
        other = stored[lo : lo + size]
        maximum = float(np.max(np.abs(computed - other)))
        assert maximum == 0
        maxdelta = max(maxdelta, maximum)
        use = selected_dose[lo : lo + size]
        np.add.at(sums, unit[lo : lo + size][use], computed[use])
        np.add.at(n, unit[lo : lo + size][use], 1)
    means = np.divide(
        sums, n[:, None], out=np.full_like(sums, np.nan), where=n[:, None] > 0
    ).reshape(len(lines), len(comps), -1)
    r3 = comparison(D.R3, means[np.ix_(*indexes)])
    counts = comparison(lm("N3_P1_d111.npy"), n.reshape(len(lines), len(comps)))
    # Re-average the complete source signature extract indexed by its matching metadata.
    signatures = pd.read_csv(e / "results/E/E_panel_signatures.tsv.gz", sep="\t")
    signatures = signatures.loc[signatures.panel.eq("P1")].set_index("sig_id")
    order = pd.read_csv(saved / "l5_P1_rows.tsv", sep="\t")
    records = lm("l5_P1.npy")
    assert signatures.index.is_unique and order.sig_id.is_unique
    assert set(signatures.index) == set(order.sig_id) and len(order) == len(records)
    sm = signatures.loc[order.sig_id]
    use = sm.pert_idose.eq("1.11 uM").to_numpy()
    assert sm.pert_time.eq(24).all()
    groups = (
        sm.cell_iname.map({v: i for i, v in enumerate(lines)}) * len(comps)
        + sm.cmap_name.map({v: i for i, v in enumerate(comps)})
    ).to_numpy()
    counts5 = np.bincount(groups[use], minlength=len(lines) * len(comps))
    sums.fill(0)
    for start in range(0, len(order), 256):
        section = slice(start, start + 256)
        selected_rows = use[section]
        np.add.at(sums, groups[section][selected_rows], records[section][selected_rows])
    avg5 = np.divide(
        sums,
        counts5[:, None],
        out=np.full_like(sums, np.nan),
        where=counts5[:, None] > 0,
    )
    avg5 = avg5.reshape(len(lines), len(comps), -1)
    means5 = comparison(lm("Z5_P1_d111.npy"), avg5)
    cells5 = comparison(lm("N5_P1_d111.npy"), counts5.reshape(len(lines), len(comps)))
    e5 = avg5[np.ix_(*indexes)]
    rounded = comparison(D.P.R, e5.astype(np.float32))
    same = np.isfinite(e5)
    rawdelta = np.abs(e5[same] - D.P.R[same])
    return dict(
        all_P1_treated_wells=len(wt),
        all_P1_DMSO_wells=len(wc),
        all_treatment_plates=wt.det_plate.nunique(),
        plates_with_two_or_more_controls=len(permitted),
        retained_treatment_wells=len(selected),
        selected_1p11_uM_wells=int(n.sum()),
        all_retained_well_marker_vectors_equal=True,
        max_abs_vehicle_centered_well_delta=maxdelta,
        numerical_level3_after_reaggregation=r3,
        selected_well_counts_per_cell=counts,
        all_P1_level5_signatures=len(records),
        selected_1p11_uM_signatures=int(counts5.sum()),
        level5_counts_per_cell=cells5,
        level5_signature_reaggregation=means5,
        independent_level5_after_float32_storage=rounded,
        original_E_float64_mean_vs_stored_D_float32_max=float(np.max(rawdelta)),
        original_E_float64_mean_vs_stored_D_float32_nonexact=int((rawdelta != 0).sum()),
    )


def table_differences(before, after):
    x = pd.read_csv(before, sep="\t", float_precision="round_trip").set_index(
        "cell_iname"
    )
    y = pd.read_csv(after, sep="\t", float_precision="round_trip").set_index(
        "cell_iname"
    )
    assert x.index.equals(y.index) and list(x) == list(y)
    diff = {}
    for k in x:
        same = x[k].eq(y[k]) | (x[k].isna() & y[k].isna())
        if not same.all():
            xv = x.loc[~same, k]
            yv = y.loc[~same, k]
            diff[k] = dict(
                rows=int((~same).sum()), max_abs_delta=float((xv - yv).abs().max())
            )
    return dict(total_rows=len(x), total_columns=len(x.columns) + 1, differences=diff)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--engine", required=True, type=Path)
    parser.add_argument(
        "--frozen-g",
        required=True,
        type=Path,
        help="Root containing inputs, tensors/D and tensors/E",
    )
    parser.add_argument(
        "--checkpoint-manifest", type=Path, help="The current power checkpoint run.json"
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root, engine, frozen = (
        p.resolve() for p in [args.run_root, args.engine, args.frozen_g]
    )
    K = load_power(root / "power/scripts/common.py", engine, frozen)
    old = K.Data()
    new, ix = source_context(K, root, old)
    inputs = root / "part_b/run"
    dd = inputs / "D"
    ee = inputs / "E"
    labels = inputs / "D_labels/results/D"
    hashes = {}

    def checkfile(key, left, right):
        a, b = digest(left), digest(right)
        hashes[key] = {
            "frozen_sha256": a,
            "source_rebuilt_sha256": b,
            "byte_identical": a == b,
        }
        return a == b

    for fname in ["P1_resp.npy", "P1_DMSO_resp.npy"]:
        assert checkfile(
            "D/" + fname, frozen / "tensors/D" / fname, dd / "data/D/tensors" / fname
        )
    for fname in [
        "R3_P1_d111.npy",
        "Z5_P1_d111.npy",
        "N3_P1_d111.npy",
        "N5_P1_d111.npy",
        "P1_d111_lines.txt",
        "P1_d111_compounds.txt",
        "l3_trt_P1.npy",
        "l3_ctl_P1.npy",
        "r3_wells_P1.npy",
        "l3_trt_P1_rows.tsv",
        "l3_ctl_P1_rows.tsv",
        "r3_wells_P1_rows.tsv",
        "l5_P1.npy",
        "l5_P1_rows.tsv",
    ]:
        assert checkfile(
            "E/" + fname, frozen / "tensors/E" / fname, ee / "data/E/tensors" / fname
        )
    for lhs, rhs in [
        ("taskD/tensors/P1_lines.tsv", dd / "results/D/tensors/P1_lines.tsv"),
        ("taskD/tensors/P1_compounds.tsv", dd / "results/D/tensors/P1_compounds.tsv"),
        ("taskD/D_landmarks.tsv", dd / "results/D/D_landmarks.tsv"),
        ("taskE/E_landmarks.tsv", ee / "results/E/E_landmarks.tsv"),
    ]:
        assert checkfile(lhs, frozen / "inputs" / lhs, rhs)
    diffs = {}
    for filename in ["D_line_covariates.tsv", "D_line_table.tsv"]:
        path = frozen / "inputs/taskD" / filename
        checkfile("taskD/" + filename, path, labels / filename)
        diffs[filename] = table_differences(path, labels / filename)
        assert not (
            set(diffs[filename]["differences"])
            - {"ploidy", "ploidy_wes", "ploidy_snp6"}
        )
    for f in ["G_sets_counts.tsv", "G_sets_lines.tsv"]:
        assert checkfile(f, engine / "results/G" / f, root / "part_b/results/model" / f)
    arrays = {}
    for key, a, b in [
        ("line_order", old.P.lines, new.P.lines),
        ("compound_order", old.P.comps, new.P.comps),
        ("gene_order", old.P.genes, new.P.genes),
        ("P1_lineages", old.P.lin, new.P.lin),
        ("level5_treatment", old.P.R, new.P.R),
        ("level5_vehicle_loaded", old.P.dmso, new.P.dmso),
        ("level3_response", old.R3, new.R3),
        ("level3_availability", old.avail3, new.avail3),
        ("level5_availability", old.avail5, new.avail5),
        ("TP53_design_columns", old.tpX, new.tpX),
        ("frozen_model_permutation_keys", old.P.keys, new.P.keys),
    ]:
        arrays[key] = comparison(
            a,
            b,
            full_digest=key.endswith("_order")
            or key == "frozen_model_permutation_keys",
        )
    for key, x, y, ma, mb in [
        ("level5", old.P.R, new.P.R, old.avail5, new.avail5),
        ("level3", old.R3, new.R3, old.avail3, new.avail3),
    ]:
        arrays[key + "_unplanted_response_sizes"] = comparison(
            K.C.response_size(x, ma), K.C.response_size(y, mb)
        )
        template1 = K.template_classes(x)
        template2 = K.template_classes(y)
        for name, a, b in zip(
            ["template_norms", "weak_compound_indices", "strong_compound_indices"],
            template1,
            template2,
        ):
            arrays[key + "_" + name] = comparison(a, b)
    real = pd.read_csv(labels / "D_line_table.tsv", sep="\t")
    aligned = K.GS.Aligned(new.P, real)
    fixed_settings = K.GS.primary_settings(old.cfg)
    frozen_labels = pd.read_csv(old.I.real_line_table_path(), sep="\t")
    frozen_aligned = K.GS.Aligned(old.P, frozen_labels)
    for term in K.TERMS:
        frozen_design = frozen_aligned.design(term, fixed_settings)
        rebuilt_design = aligned.design(term, fixed_settings)
        for attr in ["Xo", "base", "mask", "lin"]:
            arrays["observed_" + term + "_" + attr] = comparison(
                frozen_design[attr], rebuilt_design[attr]
            )
        assert (
            frozen_design["k"] == rebuilt_design["k"]
            and frozen_design["col_names"] == rebuilt_design["col_names"]
        )
    line_names, cnt, design = aligned.set_and_counts(
        "factorial", K.GS.primary_settings(old.cfg)
    )
    lineages = dict(zip(new.P.lines, new.P.lin))
    pairs = [(s, lineages[s]) for s in line_names]
    refpairs, refcnt = K.factorial_set(old.P)
    assert refpairs == pairs
    pd.testing.assert_frame_equal(
        refcnt.reset_index(drop=True), cnt.reset_index(drop=True), check_dtype=False
    )
    assert old.rule == new.rule
    n = 1000
    b = 2000
    a, ma = K.draw_codes(old.P, refpairs, refcnt, n, K.SEED_PSEUDO)
    c, mc = K.draw_codes(new.P, pairs, cnt, n, K.SEED_PSEUDO)
    labs = comparison(a, c, True)
    assert np.array_equal(ma, mc)
    used = 0
    h = hashlib.sha256()
    for rep in range(n):
        for term in K.TERMS:
            d1 = K.term_design(old, a[:, rep], ma, term)
            d2 = K.term_design(new, c[:, rep], mc, term)
            for attr in ["Xo", "mask", "g1"]:
                assert np.array_equal(d1[attr], d2[attr], equal_nan=True)
                h.update(np.ascontiguousarray(d2[attr]).tobytes())
            assert (
                d1["k"] == d2["k"]
                and d1["tcols"] == d2["tcols"]
                and d1["rule"] == d2["rule"]
            )
            used += 1
    pkeys = [(old.P.keys == new.P.keys).all(), (old.P.keys[:b] == new.P.keys[:b]).all()]
    assert all(pkeys)
    manifest = args.checkpoint_manifest
    if manifest is None:
        options = list((root / "power/raw/checkpoints").glob("*/run.json"))
        if len(options) != 1:
            raise ValueError(
                "Supply --checkpoint-manifest for the active power calculation."
            )
        manifest = options[0]
    cp = json.loads(manifest.read_text())
    assert (cp["R"], cp["B"]) == (n, b)
    cfg_files = {
        **{
            f: root / "power/scripts" / f
            for f in ["power_checkpoint.py", "power.py", "common.py"]
        },
        "SHA256SUMS.txt": engine / "src/G/SHA256SUMS.txt",
        "G_defaults.json": engine / "config/G_defaults.json",
        "G_inputs_pinned.tsv": engine / "config/G_inputs_pinned.tsv",
        **{
            f: engine / "results/G" / f
            for f in ["G_sets_lines.tsv", "G_sets_counts.tsv"]
        },
    }
    assert cp["files"] == {f: digest(path) for f, path in cfg_files.items()}
    module_hashes = K.check_frozen_code()
    assert all(r["passes"] for r in module_hashes)
    loader = old.I.hash_table()
    assert loader.status.eq("equal").all()
    metadata_paths = {
        "source_level3_well_metadata": ee / "results/E/E_wells.tsv.gz",
        "source_level5_signature_metadata": ee / "results/E/E_panel_signatures.tsv.gz",
    }
    source_rebuild = well_reconstruction(root, new, ix)
    out = dict(
        scope="Input provenance for the Figure S3 planted-effect simulation. No power models or replicates executed.",
        versions={k: cp[k] for k in ["numpy", "pandas", "scipy"]},
        frozen_engine_verified_files=len(module_hashes),
        frozen_loader_verified_objects=len(loader),
        rebuilt_source_metadata_sha256={
            k: digest(p) for k, p in metadata_paths.items()
        },
        checkpoint_source_hashes=cp["files"],
        materialized_input_hashes=hashes,
        model_input_arrays=arrays,
        covariate_and_label_source_differences=diffs,
        source_response_reconstruction=source_rebuild,
        source_E_to_model_D_axis_alignment={
            name: {
                "length": len(v),
                "positions_changed": sum(i != j for i, j in enumerate(v)),
            }
            for name, v in zip(["lines", "compounds", "genes"], ix)
        },
        availability={
            "observed_level3_line_compound_pairs": int(new.avail3.sum()),
            "observed_level5_line_compound_pairs": int(new.avail5.sum()),
        },
        panel=dict(
            total_lines=len(old.P.lines),
            compounds=len(old.P.comps),
            genes=len(old.P.genes),
            factorial_lines=len(pairs),
            factorial_lineages=len(cnt),
            group_counts=dict(
                zip(
                    ["Y_only", "XIST_only", "neither", "both"],
                    [int(cnt["n_g" + str(k)].sum()) for k in range(1, 5)],
                )
            ),
            TP53_coding=K.GS.tp53_code(old.cfg),
            TP53_columns=K.GS.tp53_columns(new.P.cov, K.GS.tp53_code(old.cfg))[1],
            candidate_rules=dict(old.rule),
            treatment_condition="LINCS ASG, 24 h, exact recorded primary dose 1.11 uM",
        ),
        draws=dict(
            planned_assignments=n,
            permutations=b,
            allocation_seed=K.SEED_PSEUDO,
            scenario_seed=K.SEED_SCEN,
            all_eligible_line_groups_match=True,
            per_lineage_group_counts_match=True,
            full_assignment_arrays=labs,
            complete_term_designs_compared=used,
            term_designs_sha256=h.hexdigest(),
            module_size=K.MODULE_SIZE,
            scenarios=list(K.SCENARIOS),
            scenario_settings=dict(
                magnitude=K.M_DELTA,
                compound_counts=K.C_K,
                gene_module_compound_counts=K.G_K,
                attenuation_unaffected_lines=K.A_U,
                compound_factor=K.C_FACTOR,
                attenuation_delta=K.A_DELTA,
            ),
            engine_uses_first_B_permutation_keys=b,
            number_of_precomputed_permutation_keys=len(old.P.keys),
        ),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2) + "\n")
    print(
        json.dumps(
            {
                "panel": out["panel"],
                "availability": out["availability"],
                "source_reconstruction": source_rebuild,
                "designs": used,
                "covariate_differences": diffs,
            }
        )
    )


if __name__ == "__main__":
    main()
