"""Power for effects planted in measured responses."""

import re

import p59style as S
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
from matplotlib.lines import Line2D      # noqa: E402

import layout as V                     # noqa: E402

TERMS = [("factY", "Y term"), ("factX", "XIST term")]
# file series -> (label, colour, marker, filled, dashed)
SERIES = {
    "panel-level test: p < 0.025": ("panel-level test, p < 0.025", S.GRAY, "s", True, False),
    "E3 on level 3: 95% CI excludes 0": ("pooled difference excludes 0, level 3", S.BLUE, "o", True, False),
    "E3 on level 5: 95% CI excludes 0": ("pooled difference excludes 0, level 5", S.BLUE, "o", False, True),
    "planted compounds at q < 0.1 (share of planted)": ("share of planted compounds at q < 0.1", S.ORANGE, "o", True,
                                                        False),
    "at least one planted compound at q < 0.1": ("any planted compound at q < 0.1", S.AQUA, "D", True, False),
}
FIG_H = 4.60
LEG_ROOM = 0.95        # inches under the bottom plots: tick labels, axis label and the three-row key


def load(inputs):
    d = V.read_tsv(inputs / "power/full_run/figS3_power_values.tsv")
    d["ln"] = d["_line"]                                      # itertuples drops names that start with _
    assert set(d.figure.astype(str)) == {"S3"}
    panels = list(dict.fromkeys(d.panel))
    assert len(panels) == 3 and [p.split(":")[0] for p in panels] == ["A", "B", "C"], panels
    assert set(d.series) <= set(SERIES), set(d.series) - set(SERIES)
    assert set(d.term) == {t for t, _ in TERMS}
    assert ((d.y_lo95 <= d.y) & (d.y <= d.y_hi95)).all() and d.y.between(0, 1).all()
    assert d.n_rep.nunique() == 1 and d.B.nunique() == 1
    return dict(d=d, panels=panels, n_rep=int(d.n_rep.iloc[0]), B=int(d.B.iloc[0]))


def column_title(pname):
    """A short title and a subtitle that repeats the file's own description of the planted effect (its numbers are
    the file's)."""
    desc = pname.split(": ", 1)[1]
    if pname.startswith("A"):
        return "Every response larger", "all compounds and genes"
    if pname.startswith("B"):
        m = re.fullmatch(r"(\d+)% larger response in k compounds", desc)
        assert m, desc
        return "Larger response", f"{m.group(1)}% larger, in k compounds"
    m = re.fullmatch(r"\+(\d+) SD in a (\d+)-gene module in k compounds", desc)
    assert m, desc
    return "Gene-module shift", f"+{m.group(1)} SD in {m.group(2)} genes, in k compounds"


def build(inputs):
    src = load(inputs)
    d, panels = src["d"], src["panels"]
    fig = plt.figure(figsize=(S.FULL, FIG_H))
    axes = {}
    info = dict(src=src, rows=[], axes=axes)
    for r, (term, tname) in enumerate(TERMS):
        for c, pname in enumerate(panels):
            ax = fig.add_axes([0.1, 0.1, 0.2, 0.2])
            axes[(term, pname)] = ax
            sub = d[(d.term == term) & (d.panel == pname)]
            cat = c > 0                                       # B, C: k compounds at evenly spaced positions
            ks = sorted(sub.x.unique())
            pos = {k: i for i, k in enumerate(ks)}
            for series, g in sub.groupby("series", sort=False):
                lab, col, mk, filled, dashed = SERIES[series]
                g = g.sort_values("x")
                x = np.array([pos[v] for v in g.x]) if cat else 100 * g.x.to_numpy(float)
                z = 2.6 if col == S.GRAY else 3.0              # the grey panel test under the coloured readouts
                ax.vlines(x, g.y_lo95, g.y_hi95, color=col, linewidth=0.8, zorder=z - 0.5,
                          gid=f"{term}:{pname[0]}:{series}:interval")
                ax.plot(x, g.y, color=col, linewidth=1.2, linestyle=(0, (3, 2)) if dashed else "-", marker=mk,
                        markersize=4.2, markerfacecolor=col if filled else "white", markeredgecolor=col if not filled
                        else "white", markeredgewidth=0.8 if not filled else 0.5, zorder=z,
                        gid=f"{term}:{pname[0]}:{series}:line")
                for xv, row in zip(x, g.itertuples()):
                    info["rows"].append(dict(term=term, panel=pname, series=series, x_file=row.x, x_drawn=float(xv),
                                             y=row.y, y_lo95=row.y_lo95, y_hi95=row.y_hi95,
                                             source_line=int(row.ln)))
            ax.set_ylim(-0.03, 1.03)
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1], ["0", "0.25", "0.5", "0.75", "1"] if c == 0 else [""] * 5)
            S.hgrid(ax)
            if cat:
                ax.set_xlim(-0.4, len(ks) - 0.6)
                ax.set_xticks(range(len(ks)), [f"{int(k)}" for k in ks])
            else:
                xs = sorted(100 * sub.x.unique())
                ax.set_xlim(-4, max(xs) + 4)
                ax.set_xticks(xs, [f"{v:g}" for v in xs])
            if r == len(TERMS) - 1:
                ax.set_xlabel("larger by (%)" if not cat else "compounds with the effect (k)")
            if c == 0:
                ax.set_ylabel(f"{tname}\nshare detected", linespacing=1.15)
            if r == 0:
                title, sub = column_title(pname)
                S.panel(ax, pname[0], title, sub=sub)
                ax.texts[-2].set_gid(f"{pname[0]}:title")
                ax.texts[-1].set_gid(f"{pname[0]}:subtitle")
    # one key for the whole figure, under the plots
    handles, labels = [], []
    for series in SERIES:
        if series not in set(d.series):
            continue
        lab, col, mk, filled, dashed = SERIES[series]
        handles.append(Line2D([], [], color=col, linewidth=1.2, linestyle=(0, (3, 2)) if dashed else "-", marker=mk,
                              markersize=4.2, markerfacecolor=col if filled else "white",
                              markeredgecolor=col if not filled else "white", markeredgewidth=0.8 if not filled else 0.5))
        labels.append(lab)
    leg = fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2, frameon=False,
                     handlelength=2.2, columnspacing=1.4, fontsize=7.5, labelcolor=S.TEXT)
    leg.set_gid("legend")
    W, H = S.FULL, FIG_H
    holder = {"legend": leg}
    def place(m):
        left, right, bottom, top = m
        holder["legend"].set_bbox_to_anchor((0.5, bottom / H), transform=fig.transFigure)
        bottom = bottom + LEG_ROOM
        gap_x, gap_y = 0.22, 0.30
        w = (W - left - right - 2 * gap_x) / 3
        h = (H - top - bottom - gap_y) / 2
        for r, (term, _) in enumerate(TERMS):
            for c, pname in enumerate(panels):
                y0 = bottom + (len(TERMS) - 1 - r) * (h + gap_y)
                axes[(term, pname)].set_position([(left + c * (w + gap_x)) / W, y0 / H, w / W, h / H])

    V.fit_page(fig, place, [0.55, 0.05, 0.03, 0.32])
    V.size_frame(fig)
    return fig, info
