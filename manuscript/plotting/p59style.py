"""Matplotlib settings and label placement for the fellowship figures."""
import hashlib
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
GRAY = "#8f8e88"          # neutral marks, 'other', 'no label'
LIGHT = "#c9c8c2"         # de-emphasised marks (null distributions, background points)
BAND = "#ebeae6"          # tolerance bands, shaded regions
GRID = "#ecebe7"
TEXT, TEXT2 = "#0b0b0b", "#52514e"
FULL, HALF = 6.5, 3.15    # inches: full text width, half width

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5, "legend.frameon": False, "legend.handletextpad": 0.4, "legend.borderaxespad": 0.2,
    "axes.edgecolor": TEXT2, "axes.labelcolor": TEXT, "axes.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "axes.titleweight": "bold", "axes.titlelocation": "left", "axes.titlepad": 6,
    "xtick.color": TEXT2, "ytick.color": TEXT2, "xtick.labelcolor": TEXT, "ytick.labelcolor": TEXT,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.direction": "out", "ytick.direction": "out",
    "lines.linewidth": 1.0, "lines.markersize": 4, "patch.linewidth": 0,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.03, "figure.dpi": 150,
    "axes.unicode_minus": True,
})


def panel(ax, letter, title=None, sub=None, x=-0.02):
    """Bold panel letter and title above the axes, left-aligned; optional grey subtitle under it."""
    ax.set_title("")
    t = f"{letter}  {title}" if title else letter
    ax.annotate(t, xy=(x, 1.0), xycoords="axes fraction", xytext=(0, 16 if sub else 6), textcoords="offset points",
                ha="left", va="bottom", fontsize=8.5, fontweight="bold", color=TEXT)
    if sub:
        ax.annotate(sub, xy=(x, 1.0), xycoords="axes fraction", xytext=(0, 5), textcoords="offset points",
                    ha="left", va="bottom", fontsize=7.5, color=TEXT2)


def hgrid(ax):
    ax.yaxis.grid(True, color=GRID, linewidth=0.5, zorder=0)
    ax.set_axisbelow(True)


def cutoff(ax, value, axis="x", label=None, where=0.98, color=TEXT2):
    """Dashed reference line with a small label at its end, in text colour."""
    if axis == "x":
        ax.axvline(value, color=color, linewidth=0.7, linestyle=(0, (3, 2)), zorder=1)
        if label:
            ax.annotate(label, xy=(value, where), xycoords=("data", "axes fraction"), xytext=(3, 0),
                        textcoords="offset points", ha="left", va="top", fontsize=7, color=TEXT2)
    else:
        ax.axhline(value, color=color, linewidth=0.7, linestyle=(0, (3, 2)), zorder=1)
        if label:
            ax.annotate(label, xy=(where, value), xycoords=("axes fraction", "data"), xytext=(0, 2),
                        textcoords="offset points", ha="right", va="bottom", fontsize=7, color=TEXT2)


def place_labels(ax, xs, ys, texts, fontsize=7, avoid_xy=None, color=TEXT, leader=True, radius=(8, 14, 22, 32, 44)):
    """Place each label near its point without overlapping other labels, the points in avoid_xy, or the axes edge.
    Tries 8 directions at growing distances; draws a thin leader line when the label sits away from its point.
    Returns the list of text artists (their final positions can be read back for checks)."""
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    trans = ax.transData
    pts = np.array([trans.transform((x, y)) for x, y in (avoid_xy if avoid_xy is not None else zip(xs, ys))])
    axbb = ax.get_window_extent(rend)
    taken = []
    out = []
    dirs = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]
    for x, y, s in zip(xs, ys, texts):
        px, py = trans.transform((x, y))
        best = None
        for r in radius:
            for dx, dy in dirs:
                ox, oy = dx * r, dy * r
                ha = "left" if dx > 0 else ("right" if dx < 0 else "center")
                va = "bottom" if dy > 0 else ("top" if dy < 0 else "center")
                t = ax.annotate(s, xy=(x, y), xytext=(ox, oy), textcoords="offset points", ha=ha, va=va,
                                fontsize=fontsize, color=color, zorder=5)
                bb = t.get_window_extent(rend).expanded(1.06, 1.15)
                inside = axbb.contains(bb.x0, bb.y0) and axbb.contains(bb.x1, bb.y1)
                hit_lab = any(bb.overlaps(o) for o in taken)
                near = (pts[:, 0] > bb.x0 - 2) & (pts[:, 0] < bb.x1 + 2) & (pts[:, 1] > bb.y0 - 2) & (pts[:, 1] < bb.y1 + 2)
                if inside and not hit_lab and not near.any():
                    best = (t, bb, r)
                    break
                t.remove()
            if best:
                break
        if best is None:  # fall back to the farthest candidate up-right, accept overlap with points only
            t = ax.annotate(s, xy=(x, y), xytext=(radius[-1], radius[-1]), textcoords="offset points", ha="left",
                            va="bottom", fontsize=fontsize, color=color, zorder=5)
            best = (t, t.get_window_extent(rend), radius[-1])
        t, bb, r = best
        if leader and r > radius[0]:
            t.arrow_patch = None
            ax.annotate("", xy=(x, y), xytext=t.get_position(), textcoords="offset points",
                        arrowprops=dict(arrowstyle="-", color=TEXT2, linewidth=0.5, shrinkA=0, shrinkB=2))
        taken.append(bb)
        out.append(t)
    return out


def save(fig, stem):
    """Write stem.pdf and stem.png; return their sha256 values."""
    os.makedirs(os.path.dirname(stem) or ".", exist_ok=True)
    fig.savefig(stem + ".pdf", metadata={"CreationDate": None, "ModDate": None})  # no timestamp: byte-stable PDFs
    fig.savefig(stem + ".png", dpi=300)
    plt.close(fig)
    return {ext: hashlib.sha256(open(stem + "." + ext, "rb").read()).hexdigest() for ext in ("pdf", "png")}
