"""Calculate donor-sex group differences in drug-specific transcription on both plates."""

import hashlib
from itertools import combinations
import json
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

P = Path(__file__).resolve().parent
SOURCE = P.parent / "genomics_extension"


def compare(null, observed):
    assert np.min(abs(null - observed)) < 1e-10
    return {
        "permutations": len(null),
        "reference_mean": float(null.mean()),
        "observed_excess": float(observed - null.mean()),
        "one_sided_p": float(np.mean(null >= observed - 1e-12)),
        "central_95_reference": np.quantile(null, [0.025, 0.975]).tolist(),
    }


def main():
    source = json.loads((SOURCE / "sources.json").read_text())
    assert (
        hashlib.sha256((SOURCE / "counts.npz").read_bytes()).hexdigest()
        == source["counts_sha256"]
    )
    arrays = np.load(SOURCE / "counts.npz")
    wells = pd.read_csv(SOURCE / "wells.tsv", sep="\t")
    features = pd.read_csv(SOURCE / "features.tsv", sep="\t")
    lines = pd.read_csv(SOURCE / "lines.tsv", sep="\t")
    assert list(lines.index) == list(range(50))
    landscape = features.landmark.to_numpy()
    chosen = json.loads((SOURCE / "results.json").read_text())["primary"]["coverage"]
    results = {}
    distributions = []
    conditional = []
    pooled = []
    logs = []
    for organ, member in chosen.items():
        candidates = lines[
            lines.cell_name.isin(member["lines"])
            & lines.annotation.isin(["Male", "Female"])
        ]
        male = candidates.annotation.eq("Male").to_numpy()
        if min(male.sum(), (~male).sum()) < 2:
            continue
        marker = candidates.marker.to_numpy(bool)
        treated = []
        controls = []
        for plate in ["plate6", "plate14"]:
            same = wells[wells.plate == plate]
            treated.append(same[same.drug != "DMSO_TF"].sort_values("condition"))
            controls.append(same[same.drug == "DMSO_TF"])
        assert treated[0].condition.to_list() == treated[1].condition.to_list()
        li = candidates.index.to_numpy()
        ok = np.ones(len(treated[0]), dtype=bool)
        for tx, veh in zip(treated, controls):
            assert (arrays["cells"][veh.index][:, li] >= 200).all()
            ok &= (arrays["cells"][tx.index][:, li] >= 200).all(axis=1)
        selected = []
        independent = []
        for tx, veh in zip(treated, controls):
            at = tx.index[ok]
            counts = arrays["gene_counts"][at][:, li][:, :, landscape]
            total = arrays["umi"][at][:, li, None]
            expr = np.log2(1 + counts.astype(float) * 1e6 / total)
            v_counts = arrays["gene_counts"][veh.index][:, li][:, :, landscape].sum(
                axis=0
            )
            v_umi = arrays["umi"][veh.index][:, li].sum(axis=0)
            dmso = np.log2(1 + 1e6 * v_counts / v_umi[:, None])
            responses = expr - dmso
            centered = responses - responses.mean(axis=0)
            selected.append(np.moveaxis(centered, 1, 0).reshape(len(li), -1))
            independent.append(expr)
        gram = selected[0] @ selected[1].T / selected[0].shape[1]
        allocations = []
        allowed = []
        labels = []
        for ids in combinations(range(len(li)), male.sum()):
            mask = np.zeros(len(li), bool)
            mask[list(ids)] = True
            weights = mask / mask.sum() - (~mask) / (~mask).sum()
            fast = float(weights @ gram @ weights)
            # Direct group means on drug expression. Centering removes the vehicle term algebraically.
            direct = []
            for expr in independent:
                field = expr[:, mask].mean(axis=1) - expr[:, ~mask].mean(axis=1)
                direct.append(field - field.mean(axis=0))
            slow = float((direct[0] * direct[1]).mean())
            assert abs(fast - slow) < 1e-10
            allocations.append(fast)
            allowed.append((mask & marker).sum() == (male & marker).sum())
            labels.append(np.array_equal(mask, male))
        dist = np.array(allocations)
        restricted = dist[np.array(allowed)]
        observed = dist[np.array(labels)].item()
        distributions.append(dist)
        conditional.append(restricted)
        pooled.append(observed)
        results[organ] = {
            "eligible_drug_doses": int(ok.sum()),
            "female": int((~male).sum()),
            "male": int(male.sum()),
            "observed": observed,
            **compare(dist, observed),
            "within_marker": compare(restricted, observed),
        }
        for name in treated[0].condition[ok]:
            logs.append({"organ": organ, "condition": name})
    assert set(results) == {"Bowel", "Lung"}
    observed = float(np.mean(pooled))
    null = np.add.outer(*distributions).ravel() / 2
    stricter = np.add.outer(*conditional).ravel() / 2
    report = {
        "genes": int(landscape.sum()),
        "organs": results,
        "observed": observed,
        "unrestricted": compare(null, observed),
        "within_marker": compare(stricter, observed),
        "direct_and_matrix_implementations_agreed_for_every_organ_assignment": True,
    }
    (P / "results.json").write_text(json.dumps(report, indent=2) + "\n")
    pd.DataFrame(logs).to_csv(P / "eligible_conditions.tsv", sep="\t", index=False)
    fig, axs = plt.subplots(1, 2, figsize=(8.1, 3), layout="constrained", sharex=True)
    for ax, distribution, title in zip(
        axs, [null, stricter], ["A   Within organ", "B   Within organ and marker group"]
    ):
        ax.hist(
            distribution, bins=30, color="#71959e", edgecolor="white", linewidth=0.4
        )
        ax.axvline(observed, color="#963f21", lw=1.8, label="Observed donor labels")
        ax.set_title(title, fontsize=10, loc="left")
        ax.set_xlabel("Cross-plate product of group differences", fontsize=9)
        ax.tick_params(labelsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.legend(fontsize=8, frameon=False)
    axs[0].set_ylabel("Label assignments")
    fig.savefig(P / "donor_response.png", dpi=230, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
