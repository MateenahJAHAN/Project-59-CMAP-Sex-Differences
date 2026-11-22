"""Methods design and a descriptive decomposition of the COLO 205 boundary case."""

from pathlib import Path
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "marker_transfer"
GENES = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y", "NLGN4Y"]
BLUE, RUST, GREY = "#25556a", "#994622", "#505a5d"
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def methods():
    fig, ax = plt.subplots(figsize=(8.3, 3.8), layout="constrained")
    ax.set(xlim=(0, 1), ylim=(0, 1))
    ax.axis("off")
    ax.text(0.01, 0.95, "Tahoe-100M", fontweight="bold", fontsize=15, color=BLUE)
    ax.text(
        0.01,
        0.85,
        "50 cancer cell lines pooled across 14 plates  ·  28 vehicle wells  ·  1,316 drug wells",
        fontsize=10,
    )
    tracks = [
        (
            "RNA markers",
            "Each treated line–well",
            "Vehicle cells of the same line\non the other 13 plates",
            "Joint Y / XIST\ncall agreement",
        ),
        (
            "Autosomal RNA by marker group",
            "Plates 6 and 14, same drug\nand dose at the same position",
            "Drug minus same-plate vehicle\nacross autosomal landmarks",
            "Cross-plate correlation\nbetween lines",
        ),
        (
            "Autosomal RNA by donor label",
            "Bowel and lung, same drug and dose",
            "Male−female difference per gene,\ncentered over drugs",
            "Agreement in group\ndifferences across plates",
        ),
        (
            "Cell recovery",
            "Drug at 5 and 0.05 µM",
            "Each line’s share of cells relative\nto vehicle on its own plate",
            "High–low difference\nwithin each organ",
        ),
    ]
    for i, (label, unit, comparator, endpoint) in enumerate(tracks):
        y = 0.72 - i * 0.208
        ax.axhline(y + 0.052, xmin=0.01, xmax=0.99, color="#bfc8ca", lw=0.7)
        ax.text(0.01, y, label, fontsize=10.5, color=BLUE, fontweight="bold", va="top")
        ax.text(0.01, y - 0.056, unit, fontsize=8.8, va="top")
        ax.text(0.39, y - 0.02, comparator, fontsize=9.4, va="top")
        ax.annotate(
            "",
            xy=(0.724, y - 0.055),
            xytext=(0.687, y - 0.055),
            arrowprops={"arrowstyle": "->", "color": GREY},
        )
        ax.text(0.735, y - 0.02, endpoint, fontsize=9.4, va="top", color=BLUE)
    fig.savefig(HERE / "design.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def boundary():
    frame = pd.read_csv(SOURCE / "well_counts.tsv.gz", sep="\t")
    frame = frame[frame.cell_line_id == "CVCL_0218"].copy()
    vehicle = frame.drug == "DMSO_TF"
    reference = frame[vehicle].groupby("plate")[["total_umi", *GENES]].sum()
    avg = np.log2(1 + reference[GENES].mul(1e6 / reference.total_umi, axis=0))
    frame = frame[frame.cells >= 200].copy()
    scores = np.log2(1 + frame[GENES].mul(1e6 / frame.total_umi, axis=0))
    treated = frame.drug != "DMSO_TF"
    diff = scores.loc[treated] - avg.loc[frame.loc[treated, "plate"]].set_axis(
        scores.index[treated]
    )
    quant = diff.quantile([0.25, 0.5, 0.75]).T
    quant.columns = ["q25", "median", "q75"]
    quant.to_csv(
        HERE / "colo205_gene_shifts.tsv",
        sep="\t",
        index_label="gene",
        float_format="%.8f",
    )
    s = scores.mean(axis=1)
    assert (int(treated.sum()), int((~treated).sum())) == (1281, 28)
    cuts = json.loads((SOURCE / "independent_check_thresholds.json").read_text())
    fig, (a, b) = plt.subplots(
        1,
        2,
        figsize=(8.6, 3.3),
        layout="constrained",
        gridspec_kw={"width_ratios": [1.15, 1]},
    )
    for mask, label, color in [
        (treated, "Drug (1,281)", BLUE),
        (~treated, "Vehicle (28)", RUST),
    ]:
        sorted_values = np.sort(s[mask])
        a.step(
            sorted_values,
            np.arange(1, len(sorted_values) + 1) / len(sorted_values),
            where="post",
            lw=1.8,
            label=label,
            color=color,
        )
    for k, style in [("lineplate", "-"), ("line", "--")]:
        cut = cuts[k]["cutoffs"]["Y"]
        a.axvline(cut, lw=1, ls=style, color=GREY)
        a.text(
            cut,
            0.03,
            f"{cut:.2f}",
            rotation=90,
            ha="right",
            va="bottom",
            fontsize=8,
            color=GREY,
        )
    a.set(
        xlabel="Mean Y-gene score, log₂(CPM + 1)",
        ylabel="Fraction of wells at or below score",
        ylim=(0, 1.03),
    )
    a.legend(frameon=False, loc="upper left", fontsize=8)
    a.set_title("A   COLO 205 across wells", loc="left", fontsize=11)
    ix = np.arange(len(quant))
    b.axvline(0, color="#8f9999", lw=0.8)
    b.hlines(ix, quant.q25, quant.q75, color=BLUE, lw=2)
    b.scatter(quant["median"], ix, color=RUST, s=18, zorder=3)
    b.set(yticks=ix, yticklabels=quant.index, xlabel="Drug − vehicle, log₂(CPM + 1)")
    b.invert_yaxis()
    b.set_title("B   Shift from same-plate vehicle", loc="left", fontsize=11)
    fig.savefig(HERE / "colo205_continuous.png", dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(quant, s[treated].median(), s[~treated].median())


if __name__ == "__main__":
    methods()
    boundary()
