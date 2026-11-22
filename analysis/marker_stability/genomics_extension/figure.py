from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

P = Path(__file__).resolve().parent
stats = json.loads((P / "results.json").read_text())
q = stats["primary"]
rep = pd.read_csv(P / "same_line_reproducibility.tsv", sep="\t")
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.spines.right": False,
        "axes.spines.top": False,
        "axes.labelcolor": "#27333b",
        "text.color": "#20292e",
    }
)
fig, axs = plt.subplots(
    1,
    3,
    figsize=(10.35, 3.35),
    gridspec_kw=dict(width_ratios=[0.97, 1.18, 1.0]),
    layout="constrained",
)
a = axs[0]
a.plot(
    [0, 1],
    rep[["different", "matched"]].to_numpy().T,
    color="#9fb4bb",
    linewidth=0.8,
    alpha=0.88,
)
a.scatter([0] * len(rep), rep.different, s=13, color="#758a91")
a.scatter([1] * len(rep), rep.matched, s=13, color="#006f77")
a.plot(
    [0, 1],
    [rep.different.mean(), rep.matched.mean()],
    color="#8c401c",
    lw=2.2,
    marker="D",
    ms=4,
)
a.set(
    xticks=[0, 1],
    xticklabels=["Different drug", "Same drug"],
    xlim=(-0.28, 1.28),
    ylabel="Cross-plate correlation (r)",
)
a.set_title("A   Same cell line", loc="left", fontweight="bold", fontsize=11)
a.grid(axis="y", alpha=0.16)
a = axs[1]
contrast = list(q["organs"].values()) + [q["estimate"]]
loc = [3, 2, 1, 0]
a.axvline(0, color="#a0aaae", lw=0.8, ls="--")
a.barh(loc, contrast, color=["#508c91"] * 3 + ["#a6532b"], height=0.32)
for y, x in zip(loc, contrast):
    a.text(max(x, 0) + 0.001, y, f"{x:+.3f}", va="center", fontsize=8)
a.set(
    xlim=(-0.008, 0.055),
    yticks=loc,
    yticklabels=["Bowel (9)", "Lung (11)", "Skin (4)", "All, equal organ weight"],
    xlabel="Same minus opposite marker group (Δr)",
    xticks=[0, 0.02, 0.04],
)
a.set_title(
    "B   Different lines, same organ", loc="left", fontweight="bold", fontsize=11
)
a = axs[2]
v = np.array(stats["triamcinolone"]["median_by_gene"])
g = stats["triamcinolone"]["genes"]
xs = np.arange(len(g))
a.bar(xs - 0.18, v[0], width=0.36, color="#17888a", label="Plate 6")
a.bar(xs + 0.18, v[1], width=0.36, color="#9cb7c3", label="Plate 14")
a.set(xticks=xs, xticklabels=g, ylabel="Median change in log₂(CPM + 1)")
a.tick_params(axis="x", labelrotation=28)
a.legend(frameon=False, fontsize=8)
a.set_title("C   Glucocorticoid control", loc="left", fontweight="bold", fontsize=11)
fig.savefig(P / "landmark_responses.png", dpi=230, bbox_inches="tight")
fig.savefig(P / "landmark_responses.svg", bbox_inches="tight")
plt.close(fig)
