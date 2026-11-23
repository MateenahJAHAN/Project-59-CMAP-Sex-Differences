"""Drug-response panel tests and effect estimates."""


import p59style as S
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from matplotlib import transforms        # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

import layout as V                     # noqa: E402

TP53 = "TP53 status (comparison exposure)"
NULLS = [("factY", "A", "Y term", S.BLUE, "main|primary|P1|factY"),
         ("factX", "B", "XIST term", S.ORANGE, "main|primary|P1|factX"),
         ("TP53", "C", "TP53 status (post hoc)", S.AQUA, f"exposure|{TP53}|P1")]
BIN = 0.005
XLIM = (0.40, 1.00)
XTICKS = [0.4, 0.6, 0.8, 1.0]
PANEL_ORDER = ["P1", "P1b", "P2"]
SCALES = [("level3", "primary", "normalised expression (level 3)"), ("level5", "S12", "z-score signature (level 5)")]
D_XLIM, D_XTICKS = (-60, 30), [-60, -40, -20, 0, 20]
FIG_H = 4.55
LABEL_PT = 7
LW_95, LW_90, MS_DOT, MEW_DOT = 1.0, 3.4, 5.0, 1.0
ROW_PT = 7.5


# ------------------------------------------------------------------------------------------------ data
def load(inputs):
    source = inputs / "part_b"
    points = pd.read_csv(source / "fig4_partB_plotted.tsv", sep="\t", float_precision="round_trip")
    draws = pd.read_csv(source / "figure4_nulls.tsv.gz", sep="\t", float_precision="round_trip")
    intervals = pd.read_csv(source / "figure4_intervals.tsv", sep="\t", float_precision="round_trip")
    intervals["source_line"] = np.arange(2, len(intervals) + 2)
    nulls = {}
    for term, letter, *_ in NULLS:
        values = points.loc[(points.panel == letter) & points.value.notna()].set_index("element").value
        observed = float(values["observed"])
        distribution = draws[term].to_numpy()
        split = observed * (1 - 1e-5)
        exceedances = int((distribution >= split).sum())
        assert exceedances == int(values["n_at_or_above"])
        nulls[term] = dict(null=distribution, stat=observed, p=float(values["p"]), split=split,
                           B=len(distribution), n_ge=exceedances, nq=int(values["n_q_lt_01"]),
                           g1=int(values["n_g1"]), g0=int(values["n_g0"]), n_tested=int(values["n_tested"]))
    estimates = {(r.term, r.panel, r.scale): r for r in intervals.itertuples()}
    heads = points[(points.element == "row_label") & (points.name == "head")]
    words = {r.series: r.text.split(" (")[0] for r in heads.itertuples()}
    margin = float(points.loc[points.element == "margin", "value"].iloc[0])
    return dict(nulls=nulls, est=estimates, words=words, margin=margin)


def histogram(d):
    null, split = d["null"], d["split"]
    k0 = int(np.floor((null.min() - split) / BIN)) - 1
    k1 = int(np.floor((null.max() - split) / BIN)) + 2
    k = np.arange(k0, k1 + 1)
    edges = split + BIN * k
    i0 = int(np.flatnonzero(k == 0)[0])
    counts, _ = np.histogram(null, bins=edges)
    assert counts.sum() == d["B"] and int(counts[i0:].sum()) == d["n_ge"]
    return edges, counts, i0


def texts_for(term, d):
    if term == "factY":
        sub = f"{d['g1']} lines with Y genes, {d['g0']} with neither"
    elif term == "factX":
        sub = f"{d['g1']} lines with XIST, {d['g0']} with neither"
    else:
        sub = f"{d['g1']} mutant, {d['g0']} wild-type lines"
    nq = d["nq"]
    comp = "no compound" if nq == 0 else (f"{nq} compound" if nq == 1 else f"{nq} compounds")
    return dict(sub=sub, sub2=f"{comp} at q < 0.1 of {d['n_tested']} tested",
                observed=f"observed {d['stat']:.3f}\np = {d['p']:.4f}")


# ------------------------------------------------------------------------------------------------ drawing
def build(inputs):
    src = load(inputs)
    nulls = src["nulls"]
    hist = {t: histogram(nulls[t]) for t, *_ in NULLS}
    top = max(int(c.max()) for _, c, _ in hist.values())
    ytop = float(np.ceil(top / 0.72 / 50) * 50)          # headroom for the two-line observed labels
    ystep = 200 if ytop <= 1200 else 400
    yticks = list(range(0, int(ytop) + 1, ystep))
    fig = plt.figure(figsize=(S.FULL, FIG_H))
    info = {"src": src, "hist": {}, "rows": []}
    axh = {}
    for term, letter, title, colour, _ in NULLS:
        d, (edges, counts, i0) = nulls[term], hist[term]
        t = texts_for(term, d)
        ax = fig.add_axes([0.1, 0.6, 0.2, 0.25])
        axh[term] = ax
        ax.stairs(counts[:i0], edges[:i0 + 1], fill=True, color=S.LIGHT, linewidth=0, zorder=2, gid=f"{letter}:below")
        ax.stairs(counts[i0:], edges[i0:], fill=True, color=colour, linewidth=0, zorder=2, gid=f"{letter}:tail")
        ax.axvline(d["stat"], color=S.TEXT, linewidth=0.9, zorder=4, gid=f"{letter}:observed_line")
        right = (XLIM[1] - d["stat"]) > (d["stat"] - XLIM[0])        # the label goes on the side with more room
        ax.annotate(t["observed"], xy=(d["stat"], 1.0), xycoords=("data", "axes fraction"),
                    xytext=(3 if right else -3, 0), textcoords="offset points", ha="left" if right else "right",
                    va="top", fontsize=LABEL_PT, color=S.TEXT, linespacing=1.15, zorder=5, gid=f"{letter}:observed")
        if term == NULLS[0][0]:
            # direct label for the coloured bins, in the first panel only (the colour means the same in all three):
            # right of the tail, clear of every bin within 2 bins of its left edge
            anchor_x = d["stat"] + 0.30 * (XLIM[1] - d["stat"])
            tail_h = counts[edges[1:] > anchor_x - 2 * BIN].max()
            ax.annotate("at or above\nobserved", xy=(anchor_x, tail_h), xycoords="data", xytext=(0, 5),
                        textcoords="offset points", ha="left", va="bottom", fontsize=LABEL_PT, color=S.TEXT2,
                        linespacing=1.1, zorder=5, gid=f"{letter}:tail_label")
        ax.set_xlim(*XLIM)
        ax.set_ylim(0, ytop)
        ax.set_xticks(XTICKS, [f"{x:.1f}" for x in XTICKS])
        ax.set_yticks(yticks, [f"{y:,}" if term == "factY" else "" for y in yticks])
        ax.set_xlabel("panel statistic")
        if term == "factY":
            ax.set_ylabel("permutations per bin")
        for a, g in zip(V.panel_lines(ax, letter, title, [t["sub"], t["sub2"]]), ("title", "subtitle", "subtitle2")):
            a.set_gid(f"{letter}:{g}")
        info["hist"][term] = dict(ax=ax, letter=letter, edges=edges, counts=counts, i0=i0, texts=t, title=title)
    # ---------------------------------------------------------------- D: pooled response-size difference
    est, words, margin = src["est"], src["words"], src["margin"]
    ys = {}                                                  # a header row per panel, then one row per scale
    y = 0
    for pn in PANEL_ORDER:
        ys[(pn, "head")] = y
        y -= 1
        for scale, _, _ in SCALES:
            ys[(pn, scale)] = y
            y -= 1
    y_lo, y_hi = y + 0.4, 0.6
    axd = {}
    for term, colour, head in (("factY", S.BLUE, "Y term"), ("factX", S.ORANGE, "XIST term")):
        ax = fig.add_axes([0.3, 0.05, 0.3, 0.3])
        axd[term] = ax
        ax.set_xlim(*D_XLIM)
        ax.set_ylim(y_lo, y_hi)
        ax.add_patch(Rectangle((-margin, y_lo), 2 * margin, y_hi - y_lo, facecolor=S.BAND, edgecolor="none",
                               linewidth=0, zorder=0, gid=f"D:{term}:margin"))
        ax.plot([0, 0], [y_lo, y_hi], color=S.TEXT2, linewidth=0.6, solid_capstyle="butt", zorder=1,
                gid=f"D:{term}:zero")
        for pn in PANEL_ORDER:
            for scale, role, sw in SCALES:
                r = est[(term, pn, scale)]
                yy = ys[(pn, scale)]
                primary = pn == "P1" and scale == "level3"
                c = colour if primary else S.GRAY
                ax.plot([r.lo95, r.hi95], [yy, yy], color=c, linewidth=LW_95, solid_capstyle="butt", zorder=3,
                        gid=f"D:{term}:{pn}:{scale}:ci95")
                ax.plot([r.lo90, r.hi90], [yy, yy], color=c, linewidth=LW_90, solid_capstyle="butt", zorder=4,
                        gid=f"D:{term}:{pn}:{scale}:ci90")
                ax.plot([r.estimate], [yy], marker="o", markersize=MS_DOT, markerfacecolor=c, markeredgecolor="white",
                        markeredgewidth=MEW_DOT, linestyle="none", zorder=5, gid=f"D:{term}:{pn}:{scale}:estimate")
                info["rows"].append(dict(panel="D", term=term, element="interval", series=pn, name=scale,
                                         role=role, y=yy, estimate=float(r.estimate), lo90=float(r.lo90),
                                         hi90=float(r.hi90), lo95=float(r.lo95), hi95=float(r.hi95),
                                         n_lines=int(r.n_lines), colour="term" if primary else "GRAY",
                                         source_line=int(r.source_line)))
        ax.set_xticks(D_XTICKS, [V.fmt_minus(str(v)) for v in D_XTICKS])
        ax.spines["left"].set_visible(False)
        ax.set_yticks([])
        ax.annotate(head, xy=(0.5, 1.0), xycoords="axes fraction", xytext=(0, 3), textcoords="offset points",
                    ha="center", va="bottom", fontsize=ROW_PT, fontweight="bold", color=S.TEXT, gid=f"D:{term}:head")
        # the band's label in the first header row, where no interval is drawn
        ax.text(0, ys[(PANEL_ORDER[0], "head")], f"\u00b1{margin:g}%", ha="center", va="center", fontsize=LABEL_PT,
                color=S.TEXT2, zorder=2, gid=f"D:{term}:margin_label")
    for term in ("factY", "factX"):
        axd[term].set_xlabel("difference, % of mean response size")
    info.update(axh=axh, axd=axd, ys=ys, ytop=ytop, yticks=yticks)
    W, H = S.FULL, FIG_H
    pad = V.content_pad()                      # text anchored at the page edge starts here
    rend = fig.canvas.get_renderer()

    def width(s, **kw):
        t = fig.text(0, 0, s, **kw)
        w = t.get_window_extent(rend).width / fig.dpi
        t.remove()
        return w

    lab_w = max([width(words["P1"] + " (primary panel)", fontsize=ROW_PT, fontweight="bold")] +
                [width(sw, fontsize=ROW_PT) + 0.12 for _, _, sw in SCALES])
    row_labels = []
    for pn in PANEL_ORDER:
        row_labels.append((pn, "head", words[pn] + (" (primary panel)" if pn == "P1" else "")))
        for scale, _, sw in SCALES:
            row_labels.append((pn, scale, sw))
    info["row_labels"] = row_labels
    tax = fig.add_axes([0.0, 0.5, 0.1, 0.001])          # holds D's letter, title and subtitle
    tax.axis("off")
    S.panel(tax, "D", "Pooled response-size difference, every panel and scale",
            sub="each term against lines with neither marker; dot: estimate; thick bar: 90% interval; thin line: 95%",
            x=0)
    tax.texts[-2].set_gid("D:title")
    tax.texts[-1].set_gid("D:subtitle")
    tr_lab = transforms.blended_transform_factory(fig.transFigure, axd["factY"].transData)
    for pn, kind, txt in row_labels:
        bold = kind == "head"
        fig.text((pad + (0 if bold else 0.12)) / W, ys[(pn, kind)], txt, transform=tr_lab, ha="left", va="center",
                 fontsize=ROW_PT, fontweight="bold" if bold else "normal", color=S.TEXT if bold else S.TEXT2,
                 gid=f"D:rowlabel:{pn}:{kind}")

    def place(m):
        left, right, bottom, top = m
        h_top, gap_h = 1.18, 0.22
        wh = (W - left - right - 2 * gap_h) / 3
        y_h = H - top - h_top
        for i, (term, *_rest) in enumerate(NULLS):
            axh[term].set_position([(left + i * (wh + gap_h)) / W, y_h / H, wh / W, h_top / H])
        d_left, gap_d = pad + lab_w + 0.12, 0.22
        wd = (W - d_left - right - gap_d) / 2
        hd = y_h - 1.02 - bottom
        for i, term in enumerate(("factY", "factX")):
            axd[term].set_position([(d_left + i * (wd + gap_d)) / W, bottom / H, wd / W, hd / H])
        tax.set_position([pad / W, (bottom + hd + 0.20) / H, 0.5 / W, 0.001 / H])

    V.fit_page(fig, place, [0.45, 0.05, 0.35, 0.35])
    V.size_frame(fig)
    return fig, info
