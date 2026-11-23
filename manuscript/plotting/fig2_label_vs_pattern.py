"""Donor labels and chromosome-marker expression."""


import p59style as S
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
from matplotlib.lines import Line2D      # noqa: E402
from matplotlib.transforms import offset_copy  # noqa: E402

import layout as V                     # noqa: E402

ORDER = ["female", "male", "no label"]   # drawing order, bottom to top (see figures_v2)
MARK = {"male": dict(marker="o", color=S.BLUE, s=16),
        "female": dict(marker="s", color=S.ORANGE, s=14),
        "no label": dict(marker="^", color=S.GRAY, s=18)}
LABEL_WORD = {"male": "labelled male", "female": "labelled female", "no label": "no label"}
EDGE = 0.4
LABEL_PT = 7
LOG2 = V.LOG2 + "(CPM+1)"
PANELS = {
    "A": dict(x="y_score", y="XIST", title="Marker pattern and label",
              quad={"tl": "XIST only", "tr": "both", "bl": "neither", "br": "Y only"},
              xlim=(0.0, 7.2), ylim=(0.0, 10.6), xticks=[0, 2, 4, 6], yticks=[0, 2, 4, 6, 8, 10]),
    "B": dict(x="yp_score", y="yq_score", title="Y-arm scores",
              quad={"tl": "Yq high, Yp low", "tr": "both high", "bl": "both low", "br": "Yp high, Yq low"},
              xlim=(0.0, 7.7), ylim=(0.0, 7.45), xticks=[0, 1, 2, 3, 4, 5, 6, 7], yticks=[0, 1, 2, 3, 4, 5, 6, 7]),
}
CUT_ROWS = {"A": ("y_score", "XIST"), "B": ("yp_score", "yq_score")}
CUT_NAMES = {"A": ("Y cut-off", "XIST cut-off"), "B": ("Yp cut-off", "Yq cut-off")}
LABELLED_A = ["RKO", "A498"]             # the two lines the Results text names next to this panel
FIG_H = 3.12
GAP = 0.56


def load(inputs):
    atlas = inputs / "part_a/tables"
    d = V.read_tsv(atlas / "partA_local_states.tsv")
    gene_rows = V.read_tsv(inputs / "figures/rebuilt_plot_data/fig2_label_vs_pattern_plotted.tsv")
    genes = {r["name"]: r["text"].split(", ") for _, r in gene_rows[gene_rows.element == "gene_list"].iterrows()}
    d["label"] = d["annotation"].map({"Male": "male", "Female": "female"}).fillna("no label")
    cuts = V.read_tsv(atlas / "partA_local_cutoffs.tsv")
    cut = {p: tuple(V.one_row(cuts, level="lineplate", score=s) for s in CUT_ROWS[p]) for p in PANELS}
    audit = V.one_row(V.read_tsv(atlas / "partA_local_label_audit_summary.tsv"), level="line x plate cut-offs (primary)")
    for p, cfg in PANELS.items():
        assert d[cfg["x"]].between(*cfg["xlim"]).all() and d[cfg["y"]].between(*cfg["ylim"]).all(), p
    # B: lines with exactly one arm score at or above its cut-off
    yp_c, yq_c = (float(r["cutoff"]) for r in cut["B"])
    one_arm = (d["yp_score"] >= yp_c) != (d["yq_score"] >= yq_c)
    # placed in this order: lines with another label first (they sit in the crowded corner), then by Yq score
    labelled_b = d[one_arm].assign(_o=(d["label"] == "male").astype(int)).sort_values(["_o", "yq_score"],
                                                                                      ascending=[True, False])
    n_male = int((d["label"] == "male").sum())
    n_male_one = int((one_arm & (d["label"] == "male")).sum())
    return dict(data=d, cut=cut, audit=audit, genes=genes, labelled_b=labelled_b,
                n_male=n_male, n_male_one=n_male_one)


def quadrant_names(ax, ycut, quad, panel):
    out = {}
    spec = {"tl": ((0, 1), "axes fraction", (4, -3), "left"), "tr": ((1, 1), "axes fraction", (-4, -3), "right"),
            "bl": ((0, ycut), ("axes fraction", "data"), (4, -3), "left"),
            "br": ((1, ycut), ("axes fraction", "data"), (-4, -3), "right")}
    for key, (xy, coords, off, ha) in spec.items():
        t = ax.annotate(quad[key], xy=xy, xycoords=coords, xytext=off, textcoords="offset points", ha=ha, va="top",
                        fontsize=LABEL_PT, color=S.TEXT2, fontstyle="italic", zorder=5)
        t.set_gid(f"{panel}:quadrant:{key}")
        out[key] = t
    return out


def legend_row(fig, ax_right):
    handles = [Line2D([], [], linestyle="none", marker="")]
    labels = ["Cellosaurus label"]
    for lab in ["male", "female", "no label"]:
        m = MARK[lab]
        handles.append(Line2D([], [], linestyle="none", marker=m["marker"], markersize=float(np.sqrt(m["s"])),
                              markerfacecolor=m["color"], markeredgecolor="white", markeredgewidth=EDGE))
        labels.append(lab)
    tr = offset_copy(ax_right.transAxes, fig=fig, x=0, y=31, units="points")
    leg = fig.legend(handles, labels, loc="lower right", bbox_to_anchor=(1.0, 1.0), bbox_transform=tr, ncol=4,
                     frameon=False, handlelength=0.9, handletextpad=0.35, columnspacing=1.1, borderpad=0,
                     borderaxespad=0, fontsize=7.5, labelcolor=S.TEXT)
    leg.set_gid("legend")
    return leg


def label_text_b(r):
    return r["cell_name"] if r["label"] == "male" else f"{r['cell_name']}\n({LABEL_WORD[r['label']]})"


def build(inputs):
    src = load(inputs)
    d = src["data"]
    fig = plt.figure(figsize=(S.FULL, FIG_H))
    axes, info = {}, {"src": src, "panels": {}}
    for p, cfg in PANELS.items():
        ax = fig.add_axes([0.1, 0.1, 0.3, 0.6])
        axes[p] = ax
        pts = {}
        for lab in ORDER:
            sub = d[d["label"] == lab]
            m = MARK[lab]
            sc = ax.scatter(sub[cfg["x"]], sub[cfg["y"]], marker=m["marker"], s=m["s"], color=m["color"],
                            edgecolors="white", linewidths=EDGE, zorder=3 + ORDER.index(lab))
            sc.set_gid(f"{p}:points:{lab}")
            pts[lab] = sub
        xc, yc = (float(r["cutoff"]) for r in src["cut"][p])
        nx, ny = CUT_NAMES[p]
        S.cutoff(ax, xc, axis="x", label=f"{nx} {xc:.2f}")
        ax.lines[-1].set_gid(f"{p}:cutoff:x")
        ax.texts[-1].set_gid(f"{p}:cutoff_label:x")
        S.cutoff(ax, yc, axis="y", label=f"{ny} {yc:.2f}")
        ax.lines[-1].set_gid(f"{p}:cutoff:y")
        ax.texts[-1].set_gid(f"{p}:cutoff_label:y")
        quads = quadrant_names(ax, yc, cfg["quad"], p)
        ax.set_xlim(*cfg["xlim"])
        ax.set_ylim(*cfg["ylim"])
        ax.set_xticks(cfg["xticks"], [str(v) for v in cfg["xticks"]])
        ax.set_yticks(cfg["yticks"], [str(v) for v in cfg["yticks"]])
        g = src["genes"]
        if p == "A":
            ax.set_xlabel(f"Y score (mean of {len(g['y'])} Y genes), {LOG2}")
            ax.set_ylabel(f"XIST, {LOG2}")
            sub_text = (f"{int(src['audit']['n_without_usual_pattern'])} of {int(src['audit']['n_labelled'])} "
                        f"labelled lines lack their label's usual pattern")
        else:
            ax.set_xlabel(f"Yp score (mean of {' and '.join(g['yp'])}), {LOG2}")
            ax.set_ylabel(f"Yq score (mean of {len(g['yq'])} Yq genes), {LOG2}")
            sub_text = f"{src['n_male_one']} of {src['n_male']} male-labelled lines have one arm score high"
        S.panel(ax, p, cfg["title"], sub=sub_text)
        ax.texts[-2].set_gid(f"{p}:title")
        ax.texts[-1].set_gid(f"{p}:subtitle")
        info["panels"][p] = dict(ax=ax, points=pts, xc=xc, yc=yc, quads=quads, subtitle=sub_text, labels={})
    legend = legend_row(fig, axes["B"])
    W, H = S.FULL, FIG_H

    def place(m):
        left, right, bottom, top = m
        w = (W - left - right - GAP) / 2
        axes["A"].set_position([left / W, bottom / H, w / W, (H - top - bottom) / H])
        axes["B"].set_position([(left + w + GAP) / W, bottom / H, w / W, (H - top - bottom) / H])

    V.fit_page(fig, place, [0.45, 0.05, 0.42, 0.62])
    for p, cfg in PANELS.items():
        ax, pi = axes[p], info["panels"][p]
        groups = [(pi["points"][lab][cfg["x"]].to_numpy(), pi["points"][lab][cfg["y"]].to_numpy(),
                   MARK[lab]["s"], MARK[lab]["marker"]) for lab in ORDER]
        fixed = [t for t in ax.texts if t.get_gid() and (":cutoff_label:" in t.get_gid() or ":quadrant:" in t.get_gid())]
        if p == "A":
            rows = [d[d["cell_name"] == n].iloc[0] for n in LABELLED_A]
            items = [(float(r[cfg["x"]]), float(r[cfg["y"]]), r["cell_name"]) for r in rows]
        else:
            rows = [r for _, r in src["labelled_b"].iterrows()]
            items = [(float(r[cfg["x"]]), float(r[cfg["y"]]), label_text_b(r)) for r in rows]
        placed = V.place_direct_labels(ax, items, groups, fixed, EDGE, LABEL_PT, cut_x=[pi["xc"]], cut_y=[pi["yc"]])
        for r, (t, lead) in zip(rows, placed):
            t.set_gid(f"{p}:label:{r['cell_name']}")
            if lead is not None:
                lead.set_gid(f"{p}:leader:{r['cell_name']}")
            pi["labels"][r["cell_name"]] = t
    V.size_frame(fig)
    info["legend"], info["axes"] = legend, axes
    return fig, info
