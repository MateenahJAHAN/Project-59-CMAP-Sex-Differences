"""Figures for task F2 (all from results/F2 tables)."""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R, F = "results/F2", "figures/F2"
os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e4e3dd", "grid.linewidth": 0.6, "axes.edgecolor": "#8a8980", "savefig.dpi": 200})
STATES = ["XY", "XX-like", "Y-lost", "XIST-lost", "discordant", "ambiguous"]
COL = dict(zip(STATES, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#6250d6"]))
MRK = dict(zip(STATES, ["o", "s", "^", "v", "D", "X"]))
INK, INK2 = "#2b2b28", "#6b6a62"
st = pd.read_csv(f"{R}/F2_states.tsv", sep="\t")
lp = pd.read_csv(f"{R}/F2_line_plate_values.tsv", sep="\t")
cut = pd.read_csv(f"{R}/F2_cutoffs.tsv", sep="\t")
c = cut.set_index(["cutoff_set", "score"]).cutoff
lp = lp.merge(st[["cell_name", "state"]], left_on="line", right_on="cell_name", validate="m:1")

# 1. XIST against Y score
fig, ax = plt.subplots(1, 2, figsize=(11, 5), sharex=True, sharey=True)
for s in STATES:
    q = lp[lp.state == s]
    ax[0].scatter(q.y_score, q.xist, s=10, c=COL[s], marker=MRK[s], alpha=0.6, linewidths=0, label=f"{s} (line call)")
ax[0].axvline(c["final_line_x_plate", "y_score"], color=INK2, lw=1, ls="--")
ax[0].axhline(c["final_line_x_plate", "xist"], color=INK2, lw=1, ls="--")
ax[0].set_title("Line x plate pseudobulks (dashed: line x plate cut-offs)")
for s in STATES:
    q = st[st.state_linelevel == s]
    ax[1].scatter(q.y_score, q.xist, s=36, c=COL[s], marker=MRK[s], edgecolors="white", linewidths=0.8, label=s)
dense = (st.y_score < 2.0) & (st.xist < 2.5)  # crowded cluster: not labelled (all listed in F2_states.tsv)
for _, r in st[~dense].iterrows():
    ax[1].annotate(r.cell_name, (r.y_score, r.xist), fontsize=6, color=INK, xytext=(3, 2), textcoords="offset points")
ax[1].annotate(f"{int(dense.sum())} lines with Y score < 2 and XIST < 2.5\n(unlabelled; see F2_states.tsv)", (1.0, 3.2), fontsize=6.5, color=INK2)
ax[1].axvline(c["final_line_level", "y_score"], color=INK2, lw=1, ls=":")
ax[1].axhline(c["final_line_level", "xist"], color=INK2, lw=1, ls=":")
ax[1].axvline(c["final_line_x_plate", "y_score"], color=INK2, lw=1, ls="--")
ax[1].axhline(c["final_line_x_plate", "xist"], color=INK2, lw=1, ls="--")
ax[1].set_title("Line level (colour: call with line-level cut-offs, dotted)")
for x in ax:
    x.set_xlabel("Y score, mean log2(CPM+1) of 8 Y genes"); x.set_ylabel("XIST, log2(CPM+1)")
ax[0].legend(fontsize=7, frameon=False, loc="upper right"); ax[1].legend(fontsize=7, frameon=False, loc="upper right")
fig.tight_layout(); fig.savefig(f"{F}/F2_fig1_xist_vs_y.png"); plt.close(fig)

# 2. Yp against Yq
fig, ax = plt.subplots(figsize=(6.5, 5.5))
male = st.cellosaurus_sex == "Male"
for s in STATES:
    q = st[(st.state == s)]
    ax.scatter(q.yq_score, q.yp_score, s=36, c=COL[s], marker=MRK[s], edgecolors="white", linewidths=0.8, label=s)
for _, r in st[male].iterrows():
    ax.annotate(r.cell_name, (r.yq_score, r.yp_score), fontsize=6, color=INK, xytext=(3, 2), textcoords="offset points")
ax.axvline(c["final_line_x_plate", "yq_score"], color=INK2, lw=1, ls="--")
ax.axhline(c["final_line_x_plate", "yp_score"], color=INK2, lw=1, ls="--")
ax.set_xlabel("Yq score (DDX3Y, USP9Y, UTY, NLGN4Y, KDM5D, EIF1AY)"); ax.set_ylabel("Yp score (RPS4Y1, ZFY)")
ax.set_title("Line-level arm scores (labels: Cellosaurus Male lines; dashed: cut-offs)")
ax.legend(fontsize=7, frameon=False); fig.tight_layout(); fig.savefig(f"{F}/F2_fig2_yp_vs_yq.png"); plt.close(fig)

# 3. zero-Y fraction by depth
dc = pd.read_csv(f"{R}/F2_depth_curve.tsv", sep="\t")
fig, ax = plt.subplots(figsize=(7, 4.8))
LC = {"A549": "#2a78d6", "NCI-H460": "#eb6834", "LoVo": "#1baf7a"}
for line, col in LC.items():
    q = dc[(dc.cell_name == line) & (dc.marker == "Y")]
    if len(q):
        ax.plot(q.median_umis, q.zero_frac, "-o", color=col, lw=2, ms=5, label=f"{line}: observed")
        ax.plot(q.median_umis, q.expected_zero_frac_homogeneous, ":", color=col, lw=1.5, label=f"{line}: expected if homogeneous")
ax.set_xscale("log"); ax.set_xlabel("UMIs per cell (bin median)"); ax.set_ylabel("share of cells with 0 Y UMIs")
ax.set_title("Zero-Y fraction by depth, vehicle cells"); ax.legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{F}/F2_fig3_zero_y_by_depth.png"); plt.close(fig)

# 4. mosaic shares by estimator and depth
M = pd.read_csv(f"{R}/F2_mosaic.tsv", sep="\t")
show = M[M.cell_name.isin(["A549", "SW48", "NCI-H460", "LoVo"]) | M.flagged_b_d_e]
est = [("a_zero_frac_ge3000", "(a) >=3,000"), ("a_zero_frac_ge8000", "(a) >=8,000"), ("b_share_1000", "(b) 1,000"),
       ("b_share_2000", "(b) 2,000"), ("b_share_4000", "(b) 4,000"), ("b_share_8000", "(b) 8,000"), ("c_share", "(c) Poisson")]
fig, ax = plt.subplots(figsize=(max(7, 0.9 * len(show) + 3), 4.8))
ecol = ["#2a78d6", "#3987e5", "#eb6834", "#d95926", "#c94e1d", "#a8411a", "#1baf7a"]
for j, (k, lab) in enumerate(est):
    ax.scatter(np.arange(len(show)) + (j - 3) * 0.1, show[k], s=24, color=ecol[j], label=lab, zorder=3)
ax.axhspan(0.15, 0.85, color="#f0efec", zorder=0)
ax.set_xticks(np.arange(len(show))); ax.set_xticklabels([f"{r.cell_name}\n{r.marker}" for _, r in show.iterrows()], fontsize=7)
ax.set_ylabel("marker-negative share"); ax.set_ylim(-0.02, 1.02)
ax.set_title("Mosaic share by estimator (shaded: 0.15 to 0.85 band of the rule for (b) at 2,000)")
ax.legend(fontsize=7, frameon=False, ncol=4, loc="upper left"); fig.tight_layout(); fig.savefig(f"{F}/F2_fig4_mosaic_shares.png"); plt.close(fig)

# 5. co-occurrence odds ratios
my = M[M.marker == "Y"]
fig, ax = plt.subplots(figsize=(6, 5.5))
ax.scatter(my.e_or_xhalves_control, my.e_or_yq0_yp0, s=30, color="#2a78d6", edgecolors="white", linewidths=0.8)
for _, r in my.iterrows():
    ax.annotate(r.cell_name, (r.e_or_xhalves_control, r.e_or_yq0_yp0), fontsize=6, color=INK, xytext=(3, 2), textcoords="offset points")
lo, hi = np.nanmin(np.r_[my.e_or_xhalves_control, my.e_or_yq0_yp0]) * 0.8, np.nanmax(np.r_[my.e_or_xhalves_control, my.e_or_yq0_yp0]) * 1.25
xx = np.geomspace(lo, hi, 50)
ax.plot(xx, xx, color=INK2, lw=1, ls=":", label="equal"); ax.plot(xx, 2 * xx, color=INK2, lw=1, ls="--", label="2 x control (rule)")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel("MH odds ratio, X-gametolog thinned halves (control)"); ax.set_ylabel("MH odds ratio, Yq = 0 with Yp = 0")
ax.set_title("Co-occurrence of zero Yq and zero Yp (Y-positive lines)"); ax.legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{F}/F2_fig5_cooccurrence_or.png"); plt.close(fig)

# 6. ambient per plate
amb = pd.read_csv(f"{R}/F2_ambient.tsv", sep="\t")
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].bar(amb.plate_n, amb.y_per_1000_umis, color="#2a78d6", width=0.7)
ax[0].set_title("Y UMIs per 1,000 UMIs, Y-negative lines (foreign cells removed)")
ax[1].bar(amb.plate_n, amb.xist_per_1000_umis, color="#eb6834", width=0.7)
ax[1].set_title("XIST UMIs per 1,000 UMIs, XIST-negative lines (foreign cells removed)")
for x in ax:
    x.set_xlabel("plate"); x.set_xticks(amb.plate_n)
fig.tight_layout(); fig.savefig(f"{F}/F2_fig6_ambient_per_plate.png"); plt.close(fig)

# 7. foreign cells against abundance
fr = pd.read_csv(f"{R}/F2_foreign_cells.tsv", sep="\t")
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
for x, mk, col in [(ax[0], "Y", "#2a78d6"), (ax[1], "XIST", "#eb6834")]:
    q = fr[fr[f"foreign_{mk}_share"].notna()]
    x.scatter(q.n_vehicle_cells, q[f"foreign_{mk}_share"], s=28, color=col, edgecolors="white", linewidths=0.8, label="observed")
    x.scatter(q.n_vehicle_cells, q[f"expected_foreign_{mk}_cells_from_background"] / q.n_vehicle_cells, s=14, color=INK2, marker="x", label="expected from plate background")
    for _, r in q.iterrows():
        x.annotate(r.cell_name, (r.n_vehicle_cells, r[f"foreign_{mk}_share"]), fontsize=6, color=INK, xytext=(3, 2), textcoords="offset points")
    x.set_xscale("log"); x.set_xlabel("vehicle cells of the line"); x.set_ylabel(f"share of cells with >=3 {mk} UMIs")
    x.set_title(f"Foreign-marker cells, {mk}-negative lines"); x.legend(fontsize=7, frameon=False)
fig.tight_layout(); fig.savefig(f"{F}/F2_fig7_foreign_vs_abundance.png"); plt.close(fig)

# 8. agreement table as a plot
ab = pd.read_csv(f"{R}/F2_agreement_by_line.tsv", sep="\t")
ab = ab[ab.cutoffs == "line_x_plate_cutoffs"]
comps = ["TaskC_expressed_state", "TaskA_state", "pilot_tahoe_state"]
lines = st.sort_values(["state", "cell_name"]).cell_name.tolist()
fig, ax = plt.subplots(figsize=(6.5, 0.2 * len(lines) + 1.5))
for i, line in enumerate(lines):
    for j, cp in enumerate(comps):
        q = ab[(ab.cell_name == line) & (ab.comparison == cp)]
        if len(q):
            ok = bool(q.agree.iloc[0])
            ax.scatter(j, i, s=60, marker="o" if ok else "X", color="#1baf7a" if ok else "#e34948")
            if not ok:
                ax.annotate(f"{q.ref_state.iloc[0]}", (j, i), fontsize=6, color=INK, xytext=(7, -2), textcoords="offset points")
        else:
            ax.scatter(j, i, s=12, marker=".", color="#c3c2b7")
ax.set_yticks(range(len(lines))); ax.set_yticklabels([f"{l} ({st.set_index('cell_name').state[l]})" for l in lines], fontsize=6)
ax.set_xticks(range(len(comps))); ax.set_xticklabels(["Task C", "Task A", "pilot"]); ax.set_xlim(-0.5, len(comps) - 0.2)
ax.invert_yaxis(); ax.grid(False)
ax.set_title("F2 call against reference (o agree, X disagree with reference state shown, . no reference)", fontsize=8)
fig.tight_layout(); fig.savefig(f"{F}/F2_fig8_agreement.png"); plt.close(fig)
print("figures written")
