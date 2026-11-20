"""Figures of task G (task draft, Deliverables): written to <out>/figures/."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def qq(perm_by_term, h2, out):
    fig, axes = plt.subplots(1, len(perm_by_term), figsize=(5 * len(perm_by_term), 4.5), squeeze=False)
    for ax, (term, (permT, cols)) in zip(axes[0], perm_by_term.items()):
        t = h2[(h2.term == term) & h2.in_ranking]
        obs = np.sort(t.T_obs.values)[::-1]
        null = np.sort(permT[:, cols], axis=1)[:, ::-1]
        exp, lo, hi = null.mean(0), np.quantile(null, 0.05, axis=0), np.quantile(null, 0.95, axis=0)
        ax.fill_between(exp, lo, hi, color="0.85", label="5th to 95th percentile of the permutation null")
        ax.plot(exp, obs, "o", ms=3, color="C0", label="observed")
        m = max(exp.max(), obs.max())
        ax.plot([0, m], [0, m], "k:", lw=0.8)
        ax.set_xlabel("expected T_c (permutation null, ordered)"); ax.set_ylabel("observed T_c, ordered")
        ax.set_title(f"P1, {term}"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(out, "G_fig1_qq_Tc_P1.png"), dpi=120); plt.close(fig)


def panel_null(null_by_term, h1, out):
    fig, axes = plt.subplots(1, len(null_by_term), figsize=(5 * len(null_by_term), 4), squeeze=False)
    for ax, (term, null) in zip(axes[0], null_by_term.items()):
        r = h1[(h1.term == term) & (h1.panel == "P1") & (h1.analysis == "primary")].iloc[0]
        ax.hist(null, bins=60, color="0.7")
        ax.axvline(r.panel_stat, color="C3", label=f"observed (p = {r.panel_p:.4f})")
        ax.set_xlabel("panel statistic (mean T_c)"); ax.set_title(f"P1, {term}"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(out, "G_fig2_panel_statistic_null.png"), dpi=120); plt.close(fig)


def e3(h3, margin, out):
    t = h3[np.isfinite(h3.get("estimate", pd.Series(dtype=float)))].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7, 0.45 * max(len(t), 4) + 1))
    y = np.arange(len(t))
    ax.axvspan(-margin, margin, color="0.9", label=f"+-{margin:g}% margin")
    ax.hlines(y, t.lo95, t.hi95, color="0.5", lw=1, label="95% CI")
    ax.hlines(y, t.lo90, t.hi90, color="C0", lw=3, label="90% CI")
    ax.plot(t.estimate, y, "ko", ms=4)
    ax.axvline(0, color="k", lw=0.6)
    ax.set_yticks(y); ax.set_yticklabels([f"{a} {b} {c} ({d})" for a, b, c, d in zip(t.panel, t.term, t.scale, t.role)], fontsize=7)
    ax.set_xlabel("E3, percent of mean response size"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(out, "G_fig3_E3.png"), dpi=120); plt.close(fig)


def p2_shape(agree, out):
    fig, ax = plt.subplots(figsize=(6, 4))
    for i, (term, g) in enumerate(agree.groupby("term")):
        ax.scatter(g.spearman_P1_P2, -np.log10(g.p_perm_one_sided), s=14, label=term, color=f"C{i}")
    ax.axhline(-np.log10(0.05), color="k", lw=0.6, ls=":")
    ax.set_xlabel("Spearman correlation of P1 and P2 coefficients (978 genes)"); ax.set_ylabel("-log10 one-sided permutation p")
    ax.set_title("P1 vs P2 shape agreement, shared compounds"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(out, "G_fig4_P1_P2_shape.png"), dpi=120); plt.close(fig)


def viability(vt, out):
    fig, ax = plt.subplots(figsize=(5, 5))
    for i, (term, g) in enumerate(vt.groupby("term")):
        ax.scatter(g.E2_same_lines, g.E2_auc_adjusted, s=12, label=term, color=f"C{i}")
    lim = np.nanmax(np.abs(vt[["E2_same_lines", "E2_auc_adjusted"]].values)) if len(vt) else 1
    ax.plot([-lim, lim], [-lim, lim], "k:", lw=0.8)
    ax.set_xlabel("E2 (same lines, no viability covariate)"); ax.set_ylabel("E2 with PRISM AUC covariate")
    ax.set_title("E2 with and without the viability adjustment, P1"); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(out, "G_fig5_E2_viability.png"), dpi=120); plt.close(fig)
