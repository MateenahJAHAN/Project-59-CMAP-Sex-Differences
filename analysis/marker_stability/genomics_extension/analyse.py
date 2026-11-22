"""Compare autosomal perturbation profiles using whole-line assignments and two replicate plates."""

from itertools import combinations
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parent


def corr_vectors(values):
    centered = values - values.mean(axis=-1, keepdims=True)
    norms = np.linalg.norm(centered, axis=-1, keepdims=True)
    return np.divide(centered, norms, out=np.zeros_like(centered), where=norms > 0)


def expr(count, umi):
    return np.log2(
        1
        + np.divide(
            count,
            umi[..., None],
            out=np.zeros(count.shape, dtype=float),
            where=umi[..., None] > 0,
        )
        * 1e6
    )


def contrast(matrix, labels):
    n = len(labels)
    same0, same1, opposite = [], [], []
    for i, j in combinations(range(n), 2):
        if labels[i] != labels[j]:
            opposite.append(matrix[i, j])
        elif labels[i]:
            same1.append(matrix[i, j])
        else:
            same0.append(matrix[i, j])
    assert min(map(len, [same0, same1, opposite])) > 0
    return (np.mean(same0) + np.mean(same1)) / 2 - np.mean(opposite)


def main():
    sources = json.loads((P / "sources.json").read_text())
    assert (
        hashlib.sha256((P / "counts.npz").read_bytes()).hexdigest()
        == sources["counts_sha256"]
    )
    npz = np.load(P / "counts.npz")
    counts, n, umi = [npz[k] for k in ["gene_counts", "cells", "umi"]]
    wells = pd.read_csv(P / "wells.tsv", sep="\t")
    lines = pd.read_csv(P / "lines.tsv", sep="\t")
    genes = pd.read_csv(P / "features.tsv", sep="\t")
    meta = [wells[wells.plate.eq(plate)] for plate in ["plate6", "plate14"]]
    treated = [
        m[~m.drug.eq("DMSO_TF")].set_index("condition", drop=False) for m in meta
    ]
    assert all(d.index.is_unique for d in treated) and set(treated[0].index) == set(
        treated[1].index
    )
    conditions = sorted(treated[0].index)
    assert len(conditions) == 94
    data = []
    controls = []
    cell_counts = []
    vehicle_n = []
    for m, tx in zip(meta, treated):
        wind = (
            tx.loc[conditions]["sample"]
            .map(dict(zip(wells["sample"], wells.index)))
            .to_numpy()
        )
        vi = m.index[m.drug.eq("DMSO_TF")].to_numpy()
        assert len(vi) == 2 and list(wells.loc[vi, "sample"]) == sorted(
            wells.loc[vi, "sample"]
        )
        vehicle = expr(counts[vi].sum(axis=0), umi[vi].sum(axis=0))
        data.append(expr(counts[wind], umi[wind]) - vehicle[None, :, :])
        controls.append(
            expr(counts[vi[1]], umi[vi[1]]) - expr(counts[vi[0]], umi[vi[0]])
        )
        cell_counts.append(n[wind])
        vehicle_n.append(n[vi].min(axis=0))
    response = np.stack(data)
    cell_counts = np.stack(cell_counts)
    vmin = np.min(vehicle_n, axis=0)
    land = genes.landmark.to_numpy()
    xn = corr_vectors(response[:, :, :, land])
    vn = corr_vectors(np.stack(controls)[:, :, land])
    score = np.einsum("dig,djg->dij", xn[0], xn[1], optimize=True)
    symmetric = (score + score.transpose(0, 2, 1)) / 2
    vehicle_cross = vn[0] @ vn[1].T
    vehicle_cross = (vehicle_cross + vehicle_cross.T) / 2

    def calculate(
        min_cells=200, marker="marker", drop_colo=False, complete=False, save=False
    ):
        ok = (np.min(cell_counts, axis=0) >= min_cells) & (vmin >= 200)[None, :]
        sub = lines.primary & (vmin >= 200)
        if drop_colo:
            sub &= ~lines.cell_name.eq("COLO 205")
        portions = {}
        perms = []
        line_n = 0
        pair_row = []
        eligibility = {}
        vehicle_by = {}
        for organ, ld in lines[sub].groupby("organ", sort=True):
            ix = ld.index.to_numpy()
            lab = ld[marker].to_numpy(dtype=bool)
            if min(lab.sum(), (~lab).sum()) < 2:
                continue
            eligible = ok[:, ix]
            common = eligible.all(axis=1)
            nmatched = eligible[:, :, None] & eligible[:, None, :]
            if complete:
                nmatched &= common[:, None, None]
            coverage = nmatched.sum(axis=0)
            ii, jj = np.triu_indices(len(ix), 1)
            eligibility[organ] = {
                "lines": list(map(str, ld.cell_name)),
                "same_conditions_all_lines": int(common.sum()),
                "range_pair_conditions": [
                    int(coverage[ii, jj].min()),
                    int(coverage[ii, jj].max()),
                ],
            }
            if np.min(coverage[ii, jj]) < 20:
                continue
            selected = symmetric[:, ix][:, :, ix]
            matrix = np.divide(
                (selected * nmatched).sum(axis=0),
                coverage,
                out=np.zeros(coverage.shape),
                where=coverage > 0,
            )
            portions[organ] = float(contrast(matrix, lab))
            line_n += len(ix)
            vehicle_by[organ] = float(contrast(vehicle_cross[np.ix_(ix, ix)], lab))
            for i, j in zip(ii, jj):
                pair_row.append(
                    dict(
                        organ=organ,
                        line1=ld.iloc[i].cell_name,
                        line2=ld.iloc[j].cell_name,
                        label1=bool(lab[i]),
                        label2=bool(lab[j]),
                        n=int(coverage[i, j]),
                        r=float(matrix[i, j]),
                        vehicle=float(vehicle_cross[ix[i], ix[j]]),
                    )
                )
            if save:
                these = []
                for subset in combinations(range(len(ix)), int(lab.sum())):
                    alt = np.zeros(len(ix), dtype=bool)
                    alt[list(subset)] = True
                    these.append(contrast(matrix, alt))
                perms.append(np.array(these))
        if save:
            pd.DataFrame(pair_row).to_csv(
                P / "distinct_line_pairs.tsv", sep="\t", index=False
            )
            if perms:
                null = perms[0]
                for g in perms[1:]:
                    null = np.add.outer(null, g).ravel()
                null /= len(portions)
                tol = 1e-12
                est = np.mean(list(portions.values()))
                p_exact = float(np.mean(np.abs(null) >= abs(est) - tol))
                pd.DataFrame(
                    {
                        "q": [0, 0.025, 0.5, 0.975, 1],
                        "mean_contrast": np.quantile(null, [0, 0.025, 0.5, 0.975, 1]),
                    }
                ).to_csv(P / "permutation_quantiles.tsv", sep="\t", index=False)
        result = {
            "estimate": float(np.mean(list(portions.values()))) if portions else None,
            "organs": portions,
            "n_lines": line_n,
            "coverage": eligibility,
            "n_pairs": len(pair_row),
            "vehicle": float(np.mean(list(vehicle_by.values()))) if portions else None,
        }
        if save:
            result.update(
                exact_p=p_exact,
                labelings=len(null),
                reference_95=np.quantile(null, [0.025, 0.975]).tolist(),
                organ_omissions={
                    k: float(np.mean([v for o, v in portions.items() if o != k]))
                    for k in portions
                },
            )
        return result

    result = {"primary": calculate(save=True)}
    result["sensitivity"] = {
        "alternate_baseline": calculate(marker="alternate"),
        "exclude_COLO205": calculate(drop_colo=True),
        "min50": calculate(50),
        "min500": calculate(500),
        "common_conditions_within_organ": calculate(complete=True),
    }
    elig = (cell_counts.min(axis=0) >= 200) & (vmin >= 200)[None, :]
    reps = []
    for li in lines.index[lines.primary & (vmin >= 200)]:
        inds = np.flatnonzero(elig[:, li])
        assert len(inds) >= 2
        cross = xn[0, inds, li] @ xn[1, inds, li].T
        same = float(np.diag(cross).mean())
        other = float((cross.sum() - np.trace(cross)) / (len(inds) * (len(inds) - 1)))
        reps.append(
            {
                "name": lines.loc[li, "cell_name"],
                "organ": lines.loc[li, "organ"],
                "n": len(inds),
                "matched": same,
                "different": other,
            }
        )
    rep = pd.DataFrame(reps)
    rep.to_csv(P / "same_line_reproducibility.tsv", sep="\t", index=False)
    result["replication"] = {
        "matched": float(rep.matched.mean()),
        "other": float(rep.different.mean()),
        "delta": float((rep.matched - rep.different).mean()),
        "lines": len(rep),
        "positive_delta_lines": int((rep.matched > rep.different).sum()),
    }
    ci = next(i for i, c in enumerate(conditions) if c.startswith("Triamcinolone |"))
    ids = np.flatnonzero(elig[ci])
    controls = ["FKBP5", "DDIT4", "SGK1", "PER1"]
    genes_pos = [genes.index[genes.gene_symbol.eq(x)].item() for x in controls]
    ctrl_expr = np.take(response[:, ci, ids, :], genes_pos, axis=-1)
    records = [
        (plate, lines.cell_name.iloc[li], *ctrl_expr[ip, k, :])
        for ip, plate in enumerate(["plate6", "plate14"])
        for k, li in enumerate(ids)
    ]
    pd.DataFrame(records, columns=["plate", "cell_line"] + controls).to_csv(
        P / "triamcinolone_by_line.tsv", sep="\t", index=False
    )
    result["triamcinolone"] = {
        "lines": len(ids),
        "genes": controls,
        "median_by_gene": np.median(ctrl_expr, axis=1).tolist(),
        "median_line_module_per_plate": np.median(
            ctrl_expr.mean(axis=2), axis=1
        ).tolist(),
        "paired_median": float(np.median(ctrl_expr.mean(axis=(0, 2)))),
    }
    result["n_landmark"] = int(land.sum())
    result["n_lines_raw"] = len(lines)
    result["shared_conditions"] = len(conditions)
    result["primary_compared_conditions"] = int(
        elig[:, lines.primary.to_numpy()].sum(axis=1).astype(bool).sum()
    )
    result["cells_two_plates"] = int(n.sum())
    result["eligible_line_conditions_all_lines"] = int(elig.sum())
    result["eligible_line_conditions_primary"] = int(
        elig[:, lines.primary.to_numpy()].sum()
    )
    pd.DataFrame(
        {
            "condition": conditions,
            "eligible_all_lines": elig.sum(axis=1),
            "eligible_primary_lines": elig[:, lines.primary.to_numpy()].sum(axis=1),
        }
    ).to_csv(P / "coverage_by_condition.tsv", sep="\t", index=False)
    (P / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "sensitivity"}, indent=2))


if __name__ == "__main__":
    main()
