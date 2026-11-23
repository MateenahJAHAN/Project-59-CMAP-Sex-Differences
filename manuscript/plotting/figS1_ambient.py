"""Ambient RNA and marker-negative cells."""

import math

import p59style as S
import matplotlib.pyplot as plt           # noqa: E402
import numpy as np                        # noqa: E402
import pandas as pd                       # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D       # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, NullLocator  # noqa: E402

import layout as V                      # noqa: E402

CELLS = "marker < 3 (foreign cells removed)"
SETS = {"A": ("Y (8 genes)", "Y-negative lines with both arms low (partial-Y lines removed)"),
        "B": ("XIST", "XIST-negative lines")}
TOP_XIST = 3
FIG_W, FIG_H = S.FULL, 5.0
LOG_X = (200, 3e5)
C_YLIM = (1e-6, 1.0)
D_YLIM = (10 ** -10.45, 10 ** -2.5)
D_ZERO = 1e-10
D_BREAK = 10 ** -9.5
HALO_PT = 12
MARK = dict(obs=dict(marker="o", s=16, edgecolors="white", linewidths=0.4, zorder=4),
            exp=dict(marker="x", s=13, color=S.GRAY, linewidths=0.7, zorder=3),
            sep=dict(marker="D", s=13, color=S.GRAY, edgecolors="white", linewidths=0.4, zorder=4))
SEG = dict(colors=S.LIGHT, linewidths=0.6, zorder=2)
BAR_W = 0.5
LEG_ROOM = 0.60        # inches between the key's bottom and the bottom of C and D (tick labels, axis label, key)


def sig2(v):
    """Two significant digits, trailing zeros kept (0.0040, not 0.004)."""
    d = max(0, 1 - int(math.floor(math.log10(abs(v)))))
    return f"{v:.{d}f}"


def load(inputs):
    amb = V.read_tsv(inputs / "figures/source_data/partA_local_ambient.tsv")
    plates = {}
    for p, (marker, line_set) in SETS.items():
        q = amb[(amb.marker == marker) & (amb.line_set == line_set) & (amb.cells == CELLS) & (amb.unit_type == "plate")]
        q = q.assign(plate_n=q.unit.str.replace("plate", "", regex=False).astype(int)).sort_values("plate_n")
        assert list(q.plate_n) == list(range(1, 15)) and q.n_lines.nunique() == 1, p
        assert np.allclose(q.per_1000, q.marker_umi / q.umi * 1000, rtol=1e-5), p
        plates[p] = q
    fr = V.read_tsv(inputs / "figures/source_data/partA_local_foreign_cells.tsv")
    assert fr.cell_name.is_unique
    y = fr[fr.y_negative == True]                                  # noqa: E712
    xs = fr[fr.xist_negative == True]                              # noqa: E712
    y_main, y_sep = y[y.partial_y == False], y[y.partial_y == True]  # noqa: E712
    assert len(y_main) == int(plates["A"].n_lines.iloc[0]), "C and A must hold the same lines"
    assert len(xs) == int(plates["B"].n_lines.iloc[0]), "D and B must hold the same lines"
    for q, mk in ((y, "Y"), (xs, "XIST")):
        assert np.allclose(q[f"share_foreign_{mk}"], q[f"n_foreign_{mk}"] / q.n_vehicle_cells_before_filter, rtol=1e-5)
    top = xs.sort_values("share_foreign_XIST", ascending=False)
    assert top.share_foreign_XIST.iloc[TOP_XIST - 1] > top.share_foreign_XIST.iloc[TOP_XIST], "tie at the label cut"
    return dict(plates=plates, y_main=y_main, y_sep=y_sep, xist=xs, top_xist=list(top.cell_name.iloc[:TOP_XIST]))


def pow10(v, _pos=None):
    e = int(round(np.log10(v)))
    return "1" if e == 0 else rf"$\mathdefault{{10^{{{e}}}}}$"


def halo(ax, xs, ys, r_pt):
    px = ax.figure.dpi / 72.0
    centre = ax.transData.transform(np.c_[xs, ys])
    rings = [centre]
    for r in np.arange(1.0, r_pt + 0.01, 1.0):
        th = np.linspace(0, 2 * np.pi, max(8, int(np.ceil(2 * np.pi * r))), endpoint=False)
        ring = centre[:, None, :] + r * px * np.stack([np.cos(th), np.sin(th)], axis=-1)[None, :, :]
        rings.append(ring.reshape(-1, 2))
    return ax.transData.inverted().transform(np.vstack(rings))


def axis_break(ax, ylim):
    lo, hi = np.log10(ylim)
    f = (np.log10(D_BREAK) - lo) / (hi - lo)
    d = 0.016
    ax.plot([0, 0], [f - d, f + d], transform=ax.transAxes, color="white", linewidth=2.4, clip_on=False, zorder=5,
            solid_capstyle="butt", gid="axis_break_gap")
    for yy in (f - d, f + d):
        ax.plot([-0.011, 0.011], [yy - 0.011, yy + 0.011], transform=ax.transAxes, color=S.TEXT2, linewidth=0.6,
                clip_on=False, zorder=6, gid="axis_break")


def plate_panel(ax, q, color, letter, title, ylabel, kind, yticks, ytop):
    x, y = q.plate_n.to_numpy(), q.per_1000.to_numpy()
    assert y.max() < ytop
    ax.bar(x, y, width=BAR_W, color=color, zorder=2, gid=f"{letter}:bars")
    ax.set_xlim(0.4, 14.6)
    ax.set_xticks(x)
    ax.set_ylim(0, ytop)
    ax.yaxis.set_major_locator(FixedLocator(yticks))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:g}"))
    S.hgrid(ax)
    sub = f"{int(q.n_lines.iloc[0])} {kind}; plates {sig2(y.min())} to {sig2(y.max())}"
    S.panel(ax, letter, title, sub)
    ax.texts[-2].set_gid(f"{letter}:title")
    ax.texts[-1].set_gid(f"{letter}:subtitle")
    ax.set_xlabel("plate")
    ax.set_ylabel(ylabel)
    return sub


def share_panel(ax, q, mk, color, letter, title, sub, ylim, yticks, yminor, label_lines, sep=None):
    x = q.n_vehicle_cells_before_filter.to_numpy(float)
    obs = q[f"share_foreign_{mk}"].to_numpy(float)
    exp = q[f"expected_share_{mk}"].to_numpy(float)
    pos = obs > 0
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*LOG_X)
    ax.set_ylim(*ylim)
    ax.add_collection(LineCollection([[(a, e), (a, o)] for a, e, o in zip(x[pos], exp[pos], obs[pos])],
                                     gid=f"{letter}:segment", **SEG))
    ax.scatter(x, exp, gid=f"{letter}:expected", **MARK["exp"])
    ax.scatter(x[pos], obs[pos], color=color, gid=f"{letter}:observed", **MARK["obs"])
    if (~pos).any():
        assert ylim[0] < D_ZERO < D_BREAK
        ax.scatter(x[~pos], np.full((~pos).sum(), D_ZERO), color=color, gid=f"{letter}:observed_zero", **MARK["obs"])
        axis_break(ax, ylim)
    sx = sy = np.array([])
    if sep is not None:
        sx = sep.n_vehicle_cells_before_filter.to_numpy(float)
        sy = sep[f"share_foreign_{mk}"].to_numpy(float)
        assert (sy > 0).all()
        ax.scatter(sx, sy, gid=f"{letter}:separate", **MARK["sep"])
    ax.yaxis.set_major_locator(FixedLocator(yticks))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, p: "0" if v == D_ZERO else pow10(v)))
    ax.yaxis.set_minor_locator(FixedLocator(yminor) if yminor else NullLocator())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.xaxis.set_major_locator(FixedLocator([1e3, 1e4, 1e5]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:,.0f}"))
    ax.xaxis.set_minor_formatter(NullFormatter())
    S.panel(ax, letter, title, sub)
    ax.texts[-2].set_gid(f"{letter}:title")
    ax.texts[-1].set_gid(f"{letter}:subtitle")
    ax.set_xlabel("vehicle cells of the line")
    ax.set_ylabel("share of the line's cells")
    return dict(x=x, obs=obs, exp=exp, pos=pos, sx=sx, sy=sy)


def share_labels(ax, q, mk, letter, label_lines, drawn, sep=None):
    """Labels for the named lines, placed after the layout is final; points that coincide on the page (within 2 pt)
    share one label. Labels keep clear of every mark and stay HALO_PT from every dot they do not name."""
    x, obs, exp, pos, sx, sy = (drawn[k] for k in ("x", "obs", "exp", "pos", "sx", "sy"))
    allq = pd.concat([q, sep]) if sep is not None else q
    lab = allq[allq.cell_name.isin(label_lines)].sort_values(f"share_foreign_{mk}", ascending=False)
    assert len(lab) == len(label_lines)
    fig = ax.figure
    fig.canvas.draw()
    to_pt = 72.0 / fig.dpi
    disp = ax.transData.transform(np.c_[lab.n_vehicle_cells_before_filter, lab[f"share_foreign_{mk}"]]) * to_pt
    groups = []
    for i in range(len(lab)):
        for g in groups:
            if np.hypot(*(disp[i] - disp[g[0]])) < 2.0:
                g.append(i)
                break
        else:
            groups.append([i])
    lx = [lab.n_vehicle_cells_before_filter.iloc[g[0]] for g in groups]
    ly = [lab[f"share_foreign_{mk}"].iloc[g[0]] for g in groups]
    texts = ["\n".join(lab.cell_name.iloc[g]) for g in groups]
    y_drawn = np.where(pos, obs, D_ZERO)
    marks = [(x, exp, MARK["exp"]["s"], "x"), (x, y_drawn, MARK["obs"]["s"], "o")]
    if len(sx):
        marks.append((sx, sy, MARK["sep"]["s"], "D"))
    seg_pts = [np.c_[np.full(60, a), np.geomspace(e, o, 60)] for a, e, o in zip(x[pos], exp[pos], obs[pos])]
    shared = [[(lab.n_vehicle_cells_before_filter.iloc[i], lab[f"share_foreign_{mk}"].iloc[i]) for i in g[1:]]
              for g in groups]
    placed = V.place_direct_labels(ax, list(zip(lx, ly, texts)), marks, [], 0.4, 7,
                                   extra_xy=np.vstack(seg_pts) if seg_pts else None, shared=shared)
    arts = [t for t, _ in placed]
    for (t, lead), g in zip(placed, groups):
        t.set_gid(f"{letter}:label:" + "|".join(lab.cell_name.iloc[g]))
        if lead is not None:
            lead.set_gid(f"{letter}:leader:" + "|".join(lab.cell_name.iloc[g]))
    return {(letter, name): text for g, text in zip(groups, texts) for name in lab.cell_name.iloc[g]}


def build(inputs):
    d = load(inputs)
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    ax = {k: fig.add_axes([0.1, 0.1, 0.3, 0.3]) for k in "ABCD"}
    subs = {}
    subs["A"] = plate_panel(ax["A"], d["plates"]["A"], S.BLUE, "A", "Ambient Y RNA, marker-rich cells removed",
                            "Y UMIs per 1,000 UMIs", "lines with both Y-arm scores low", [0, 0.01, 0.02, 0.03], 0.038)
    subs["B"] = plate_panel(ax["B"], d["plates"]["B"], S.ORANGE, "B", "Ambient XIST RNA, marker-rich cells removed",
                            "XIST UMIs per 1,000 UMIs", "lines with XIST low", [0, 0.001, 0.002, 0.003, 0.004], 0.0044)
    subs["C"] = f"{len(d['y_main'])} lines with both Y-arm scores low"
    subs["D"] = f"{len(d['xist'])} lines with XIST low"
    W, H = FIG_W, FIG_H
    LEFT_GAP = 0.66

    holder = {}

    def place(m):
        left, right, bottom, top = m
        w = (W - left - right - LEFT_GAP) / 2
        h_top, h_bot, gap_v = 1.30, 1.72, 0.92
        if "legend" in holder:                       # the key sits at the bottom margin, under C and D
            holder["legend"].set_bbox_to_anchor((0.5, bottom / H), transform=fig.transFigure)
        bottom = bottom + LEG_ROOM
        ax["C"].set_position([left / W, bottom / H, w / W, h_bot / H])
        ax["D"].set_position([(left + w + LEFT_GAP) / W, bottom / H, w / W, h_bot / H])
        yb = bottom + h_bot + gap_v
        ax["A"].set_position([left / W, yb / H, w / W, (H - top - yb) / H])
        ax["B"].set_position([(left + w + LEFT_GAP) / W, yb / H, w / W, (H - top - yb) / H])

    place([0.47, 0.05, 0.62, 0.30])
    drawn = {"C": share_panel(ax["C"], d["y_main"], "Y", S.BLUE, "C", "Cells with 3 or more Y UMIs", subs["C"],
                              C_YLIM, [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0], None,
                              list(d["y_sep"].cell_name), sep=d["y_sep"]),
             "D": share_panel(ax["D"], d["xist"], "XIST", S.ORANGE, "D", "Cells with 3 or more XIST UMIs", subs["D"],
                              D_YLIM, [D_ZERO, 1e-9, 1e-7, 1e-5, 1e-3], [1e-8, 1e-6, 1e-4], d["top_xist"])}
    dot = lambda c: Line2D([], [], ls="", marker="o", ms=4, mfc=c, mec="white", mew=0.4)       # noqa: E731
    cross = Line2D([], [], ls="", marker="x", ms=3.6, mec=S.GRAY, mew=0.7)
    diamond = Line2D([], [], ls="", marker="D", ms=3.6, mfc=S.GRAY, mec="white", mew=0.4)
    leg = fig.legend([(dot(S.BLUE), dot(S.ORANGE)), cross, diamond],
                     ["observed", "expected from the plate background", "Partial-Y lines shown for comparison"],
                     handler_map={tuple: HandlerTuple(ndivide=None, pad=0.15)}, loc="lower center",
                     bbox_to_anchor=(0.5, 0.0), ncol=3, columnspacing=1.6, handlelength=1.8, labelcolor=S.TEXT)
    leg.set_gid("legend")
    holder["legend"] = leg
    V.fit_page(fig, place, [0.47, 0.05, 0.03, 0.30])
    labels = {}
    labels.update(share_labels(ax["C"], d["y_main"], "Y", "C", list(d["y_sep"].cell_name), drawn["C"], sep=d["y_sep"]))
    labels.update(share_labels(ax["D"], d["xist"], "XIST", "D", d["top_xist"], drawn["D"]))
    V.size_frame(fig)
    return fig, dict(ax=ax, d=d, labels=labels, subs=subs)
