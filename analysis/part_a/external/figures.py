"""Step 5: figures from results/F1 tables (plan section 11)."""
import sys, os, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "verification/source_replay/part_a/scripts"))
from common import Y_GENES, YP_GENES, YQ_GENES

ROOT = sys.argv[1]
RES = os.path.join(ROOT, "results/F1"); FIG = os.path.join(ROOT, "figures/F1"); os.makedirs(FIG, exist_ok=True)
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
ANN = {"Male": ("#2a78d6", "o"), "Female": ("#eb6834", "s"), "other": ("#8a8984", "^")}
GOOD, CRIT, NEUTRAL = "#0ca30c", "#d03b3b", "#8a8984"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                     "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF,
                     "axes.facecolor": SURF, "savefig.facecolor": SURF})

em = pd.read_csv(os.path.join(RES, "F1_emerald_states.tsv"), sep="\t")
sp = pd.read_csv(os.path.join(RES, "F1_sciplex_states.tsv"), sep="\t")
cut = pd.read_csv(os.path.join(RES, "F1_cutoffs.tsv"), sep="\t")
ec = cut[cut.dataset == "EmeraldBay"].set_index("score").cutoff
sc = cut[cut.dataset != "EmeraldBay"].set_index("score").cutoff


def annkey(a):
    return a if a in ("Male", "Female") else "other"


def scatter(ax, df, xcol, ycol, label_col, cx, cy, title, label_all=True):
    for k, (col, mk) in ANN.items():
        d = df[df.annotation.map(annkey) == k]
        lab = {"Male": "Male (Cellosaurus)", "Female": "Female (Cellosaurus)", "other": "Unspecified or unknown"}[k]
        ax.scatter(d[xcol], d[ycol], s=36, c=col, marker=mk, edgecolors=SURF, linewidths=1.5, label=f"{lab} (n={len(d)})", zorder=3)
    ax.axvline(cx, color=INK2, lw=1, ls="--", zorder=1); ax.axhline(cy, color=INK2, lw=1, ls="--", zorder=1)
    if label_all:
        for _, r in df.iterrows():
            ax.annotate(str(r[label_col]), (r[xcol], r[ycol]), fontsize=6, color=INK2, xytext=(3, 2), textcoords="offset points")
    ax.set_title(title, fontsize=10, color=INK, loc="left")


# ---- Figure 1: XIST vs Y score ----
fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1.6, 1]})
scatter(axes[0], em, "XIST", "y_score", "cell_name", ec["XIST"], ec["y_score"],
        f"EmeraldBay vehicle (DMSO_TF), 52 lines. Dashed: KDE cut-offs (XIST {ec['XIST']:.2f}, Y {ec['y_score']:.2f})")
axes[0].set_xlabel("XIST, log2(CPM+1), line pseudobulk"); axes[0].set_ylabel("Y score, mean log2(CPM+1) of 8 Y genes")
axes[0].legend(frameon=False, fontsize=8, loc="upper right")
sp2 = sp.assign(lab=sp["sample"] + " " + sp.cell_name)
scatter(axes[1], sp2, "XIST", "y_score", "lab", sc["XIST"], sc["y_score"],
        f"sci-Plex vehicle. Dashed: Task A CCLE cut-offs (XIST {sc['XIST']:.2f}, Y {sc['y_score']:.2f})")
axes[1].set_xlabel("XIST, log2(CPM+1)"); axes[1].set_ylabel("Y score")
axes[1].legend(frameon=False, fontsize=8, loc="upper right")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "F1_fig1_xist_vs_y.png"), dpi=160); plt.close(fig)

# ---- Figure 2: arm-level scores ----
fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), gridspec_kw={"width_ratios": [1.6, 1]})
scatter(axes[0], em, "yp_score", "yq_score", "cell_name", ec["yp_score"], ec["yq_score"],
        f"EmeraldBay: Yp (RPS4Y1, ZFY) vs Yq (6 genes). Cut-offs Yp {ec['yp_score']:.2f}, Yq {ec['yq_score']:.2f}")
axes[0].set_xlabel("Yp score, log2(CPM+1)"); axes[0].set_ylabel("Yq score, log2(CPM+1)")
axes[0].legend(frameon=False, fontsize=8, loc="lower right")
scatter(axes[1], sp2, "yp_score", "yq_score", "lab", sc["yp_score"], sc["yq_score"],
        f"sci-Plex: cut-offs {sc['yp_score']:.2f} (Task A Y cut-off for both arms)")
axes[1].set_xlabel("Yp score"); axes[1].set_ylabel("Yq score"); axes[1].legend(frameon=False, fontsize=8, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "F1_fig2_arm_scores.png"), dpi=160); plt.close(fig)

# ---- Figure 3: per-cell Y UMIs, A549 vs a homogeneous XY comparator ----
cells = pd.read_csv(os.path.join(RES, "F1_emerald_vehicle_cells.tsv.gz"), sep="\t")
cells = cells[cells.drug == "DMSO_TF"]
xy = em[(em.state == "XY") & (em.confidence == "high") & (em.cellosaurus_id != "CVCL_0023")]
comp = xy.sort_values("n_vehicle_cells", ascending=False).iloc[0]
json.dump({"comparator_cellosaurus_id": comp.cellosaurus_id, "comparator_name": comp.cell_name,
           "comparator_n_vehicle_cells": int(comp.n_vehicle_cells),
           "rule": "EmeraldBay XY line with high confidence and the most vehicle cells, excluding A549 (plan 11)"},
          open(os.path.join(RES, "F1_fig3_comparator.json"), "w"), indent=1)
mos = pd.read_csv(os.path.join(RES, "F1_emerald_mosaic.tsv"), sep="\t")
spv = []
for tag in ("sciPlex2", "sciPlex3", "sciPlex4"):
    v = pd.read_parquet(os.path.join(ROOT, f"data/F1/sciplex_proc/{tag}_vehicle_cells.parquet"))
    spv.append(v[v.line == "A549"].assign(sample=tag))
spv = pd.concat(spv)
fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
pairs = [("A549", "CVCL_0023", "#2a78d6"), (comp.cell_name, comp.cellosaurus_id, "#eb6834")]
for nm, cid, col in pairs:
    d = cells[cells.line == cid]; y = d[Y_GENES].sum(axis=1); rate = 1000 * y / d.depth
    axes[0].hist(rate.clip(upper=15), bins=np.linspace(0, 15, 61), histtype="step", lw=2, color=col, density=True,
                 label=f"{nm} (n={len(d)}; {100*np.mean(y==0):.1f}% cells with 0 Y UMIs)")
    dec = pd.qcut(d.depth, 10, labels=False, duplicates="drop")
    g = pd.DataFrame({"dec": dec, "zero": (y == 0), "depth": d.depth, "y": y}).groupby("dec")
    lam = y.sum() / d.depth.sum()
    obs = g.zero.mean(); md = g.depth.median()
    exp1 = g.apply(lambda s: np.mean(np.exp(-lam * s.depth)))
    m = mos[(mos.line == cid) & (mos.marker == "Y")]
    axes[1].plot(md, obs, "-o", color=col, lw=2, ms=6, label=f"{nm} observed")
    axes[1].plot(md, exp1, ":", color=col, lw=2, label=f"{nm} expected, one Poisson rate")
    if len(m):
        pi, l1 = m.pi_neg.iloc[0], m.lam_per_1k.iloc[0] / 1000
        b = m.background_per_1k.iloc[0] / 1000
        expm = g.apply(lambda s: np.mean(pi * np.exp(-b * s.depth) + (1 - pi) * np.exp(-l1 * s.depth)))
        axes[1].plot(md, expm, "--", color=col, lw=1.5, label=f"{nm} expected, fitted mixture (pi={pi:.2f})")
axes[0].set_xlabel("Y UMIs per 1,000 UMIs per cell (clipped at 15)"); axes[0].set_ylabel("density")
axes[0].set_title("EmeraldBay vehicle cells: per-cell Y-UMI rate", loc="left", fontsize=10)
axes[0].legend(frameon=False, fontsize=7.5)
axes[1].set_xscale("log"); axes[1].set_xlabel("median depth of decile (UMIs, log scale)"); axes[1].set_ylabel("share of cells with 0 Y UMIs")
axes[1].set_title("Zero-Y share by depth decile", loc="left", fontsize=10); axes[1].legend(frameon=False, fontsize=6.5)
spm = pd.read_csv(os.path.join(RES, "F1_sciplex_mosaic.tsv"), sep="\t")
for tag, col in (("sciPlex2", "#2a78d6"), ("sciPlex3", "#1baf7a"), ("sciPlex4", "#6250d6")):
    d = spv[spv["sample"] == tag]; y = d[Y_GENES].sum(axis=1)
    dec = pd.qcut(d.depth, 10, labels=False, duplicates="drop")
    g = pd.DataFrame({"dec": dec, "zero": (y == 0), "depth": d.depth}).groupby("dec")
    lam = y.sum() / d.depth.sum()
    axes[2].plot(g.depth.median(), g.zero.mean(), "-o", color=col, lw=2, ms=6, label=f"{tag} A549 observed (n={len(d)})")
    axes[2].plot(g.depth.median(), g.apply(lambda s: np.mean(np.exp(-lam * s.depth))), ":", color=col, lw=2, label=f"{tag} expected, one rate")
axes[2].set_xscale("log"); axes[2].set_xlabel("median depth of decile (UMIs, log scale)"); axes[2].set_ylabel("share of cells with 0 Y UMIs")
axes[2].set_title("sci-Plex A549 vehicle cells (no other XY line in sci-Plex)", loc="left", fontsize=10)
axes[2].legend(frameon=False, fontsize=6.5)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "F1_fig3_a549_y_umis.png"), dpi=160); plt.close(fig)

# ---- Figure 4: agreement ----
ag = pd.read_csv(os.path.join(RES, "F1_agreement.tsv"), sep="\t")
st = ag[ag.pct_agree.notna()].reset_index(drop=True)
sh = ag[ag.value.notna() | (ag.result == "not run")].reset_index(drop=True)
fig, axes = plt.subplots(1, 2, figsize=(14, 0.42 * max(len(st), len(sh)) + 1.6), gridspec_kw={"width_ratios": [1.3, 1]})
ax = axes[0]
for i, r in st.iterrows():
    col = GOOD if r.result == "pass" else CRIT if r.result == "fail" else NEUTRAL
    ax.barh(i, 100 * r.pct_agree, color=col, height=0.6)
    ax.text(100 * r.pct_agree + 1, i, f"{int(r.n_agree)}/{int(r.n_compared)} ({100*r.pct_agree:.0f}%) {r.result}", va="center", fontsize=7.5, color=INK)
ax.set_yticks(range(len(st))); ax.set_yticklabels([f"{c}\n[{s[:60]}]" for c, s in zip(st.comparison, st.subset.fillna(""))], fontsize=7)
ax.invert_yaxis(); ax.set_xlim(0, 135); ax.axvline(90, color=INK2, ls="--", lw=1)
ax.set_xlabel("% of lines with the same state (dashed: 90% pre-registered threshold, first row only)")
ax.set_title("State agreement (green pass, red fail, grey descriptive)", loc="left", fontsize=10)
ax = axes[1]
for i, r in sh.iterrows():
    if r.result == "not run":
        ax.text(0.02, i, "not run", va="center", fontsize=8); continue
    col = GOOD if r.result == "pass" else CRIT
    ax.axhspan(i - 0.3, i + 0.3, xmin=0, xmax=0, color=SURF)
    ax.plot([r.reference - 0.15, r.reference + 0.15], [i, i], color=GRID, lw=8, solid_capstyle="butt")
    ax.plot(r.reference, i, "|", color=INK2, ms=14, mew=2)
    ax.plot(r.value, i, "o", color=col, ms=9, mec=SURF, mew=1.5)
    ax.text(1.02, i, f"{r.value:.3f} vs {r.reference:.3f}: {r.result}", va="center", fontsize=7.5, color=INK, transform=ax.get_yaxis_transform())
ax.set_yticks(range(len(sh))); ax.set_yticklabels(sh.comparison, fontsize=7); ax.invert_yaxis(); ax.set_xlim(-0.02, 1.02)
ax.set_xlabel("marker-negative share (EM mixture). Grey band: pilot value +/- 0.15")
ax.set_title("Subpopulation shares vs Tahoe pilot", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "F1_fig4_agreement.png"), dpi=160); plt.close(fig)
print("figures written")
