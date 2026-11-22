"""Independent replay of main pair correlations and exact line-label distribution."""

from pathlib import Path
import itertools, json
import numpy as np, pandas as pd

P = Path(__file__).resolve().parent
n = np.load(P / "counts.npz")
y = n["gene_counts"]
cells = n["cells"]
umi = n["umi"]
w = pd.read_csv(P / "wells.tsv", sep="\t")
lines = pd.read_csv(P / "lines.tsv", sep="\t")
ft = pd.read_csv(P / "features.tsv", sep="\t")
lm = np.flatnonzero(ft.landmark)
res = json.loads((P / "results.json").read_text())
matrix = pd.read_csv(P / "distinct_line_pairs.tsv", sep="\t")
ids = {s: i for i, s in enumerate(lines.cell_name)}
eligible = lines[lines.primary & (cells[w.drug.eq("DMSO_TF")].min(axis=0) >= 200)]
meta = w[~w.drug.eq("DMSO_TF")].copy()
meta["i"] = meta.index
common = meta.pivot(index="condition", columns="plate", values="i")
assert common.shape == (94, 2)
response = {}
for plate in ["plate6", "plate14"]:
    ctrls = w.index[(w.plate == plate) & w.drug.eq("DMSO_TF")]
    rows = common[plate].to_numpy()
    base = np.log2(1 + 1e6 * y[ctrls].sum(0)[:, lm] / umi[ctrls].sum(0)[:, None])
    for name in eligible.cell_name:
        i = ids[name]
        a = (
            np.log2(1 + 1e6 * y[rows, i][:, lm] / np.maximum(umi[rows, i], 1)[:, None])
            - base[i]
        )
        response[plate, name] = a
results = []
for organ, g in eligible.groupby("organ"):
    names = g.cell_name.tolist()
    ng = len(g)
    inds = np.array([ids[name] for name in names])
    labels = g.marker.to_numpy()
    K = np.diag(np.zeros(ng))
    for a, b in itertools.combinations(range(ng), 2):
        i, j = names[a], names[b]
        mask = np.logical_and.reduce(
            [
                cells[common[p], ids[z]] >= 200
                for z in [i, j]
                for p in ["plate6", "plate14"]
            ]
        )
        vv = []
        for c in np.flatnonzero(mask):
            vv.append(
                (
                    np.corrcoef(response["plate6", i][c], response["plate14", j][c])[
                        0, 1
                    ]
                    + np.corrcoef(response["plate14", i][c], response["plate6", j][c])[
                        0, 1
                    ]
                )
                / 2
            )
        K[a, b] = K[b, a] = np.mean(vv)
        published = matrix[
            ((matrix.line1 == i) & (matrix.line2 == j))
            | ((matrix.line1 == j) & (matrix.line2 == i))
        ].iloc[0]
        assert len(vv) == published.n and abs(K[a, b] - published.r) < 1e-12

    def calc(x):
        pos = np.flatnonzero(x)
        neg = np.flatnonzero(~x)
        rpos = K[np.ix_(pos, pos)].sum() / (len(pos) * (len(pos) - 1))
        rneg = K[np.ix_(neg, neg)].sum() / (len(neg) * (len(neg) - 1))
        return 0.5 * (rpos + rneg) - K[np.ix_(pos, neg)].mean()

    emp = calc(labels)
    assert abs(emp - res["primary"]["organs"][organ]) < 1e-12
    null = [
        calc(np.isin(np.arange(ng), cb))
        for cb in itertools.combinations(range(ng), labels.sum())
    ]
    results.append((emp, np.asarray(null)))
observed = np.mean([r[0] for r in results])
distribution = np.mean(
    np.meshgrid(*[r[1] for r in results], indexing="ij"), axis=0
).ravel()
assert (
    abs(observed - res["primary"]["estimate"]) < 1e-12
    and len(distribution) == res["primary"]["labelings"]
)
p = float(np.mean(np.abs(distribution) >= abs(observed) - 1e-12))
assert abs(p - res["primary"]["exact_p"]) < 1e-12
print(
    json.dumps(
        {
            "input_genes": len(lm),
            "pair_correlations_independently_recomputed": len(matrix),
            "exact_labelings": len(distribution),
            "estimate": observed,
            "exact_p": p,
            "all_counts_estimates_match": True,
        },
        indent=2,
    )
)
