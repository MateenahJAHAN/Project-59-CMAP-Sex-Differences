"""Marker agreement across independent assays."""


import p59style as S
import matplotlib.pyplot as plt          # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

import layout as V                     # noqa: E402

STR_SOURCES = [("Tahoe-100M", "Tahoe-100M"), ("EmeraldBay", "EmeraldBay"), ("sci-Plex", "sci-Plex"),
               ("CCLE", "CCLE"), ("CMP", "Cell Model Passports")]
STR_COLS = [("n_compared", "lines"), ("n_agree", "agree"), ("n_type_A", "AMELY,\nY score low"),
            ("n_type_B", "no AMELY,\nY expressed")]
PATTERNS = ["Y only", "XIST only", "both", "neither"]
TICK = {"Y only": "Y\nonly", "XIST only": "XIST\nonly", "both": "both", "neither": "neither"}
Y_LEVELS, STR_LEVELS = ["Y score high", "Y score low"], ["X,Y", "X only"]
STR_TICK = {"X,Y": "X\nand Y", "X only": "X\nonly"}
CONFIDENT = ("high", "medium")
PITCH = 0.31           # inches per tile
GAP_FRAC = 0.06        # surface gap between tiles, as a fraction of the pitch
FIG_H = 4.56
NOTE_PT = 7
LEG_H = 0.16           # inches: the legend row above the panels
HEAD_H = 0.40          # inches: panel letter, title and subtitle
COL_LAB_H = 0.44       # inches: column title and column labels above the tiles
ROW_PT = 7.5
ROW_LAB_W = 0.76       # inches: row title and row labels left of the tiles
NOTE_GAP = 0.14        # inches between the tiles and the note at their right


def load(inputs):
    rows = V.read_tsv(inputs / "figures/rebuilt_plot_data/fig3_agreement_plotted.tsv", keep_default_na=False)
    lines = rows.loc[rows.element == "line"].copy()
    lines["same"] = lines["same"].astype(str).map({"True": True, "False": False})
    counts = V.read_tsv(inputs / "part_a/tables/figure3_agreement_by_source.tsv")
    return dict(lines=lines, str_counts=counts.rename(columns={"_line": "line"}).set_index("source").to_dict("index"))


PANELS = [
    dict(key="A", title="Tahoe-100M and external RNA-seq", rows="Tahoe-100M", cols="external RNA-seq",
         rlev=PATTERNS, clev=PATTERNS, unit="lines"),
    dict(key="B", title="EmeraldBay and Tahoe-100M", rows="EmeraldBay", cols="Tahoe-100M",
         rlev=PATTERNS, clev=PATTERNS, unit="shared lines"),
    dict(key="C", title="sci-Plex and external RNA-seq", rows="sci-Plex", cols="external RNA-seq",
         rlev=PATTERNS, clev=PATTERNS, unit="lines in single experiments"),
    dict(key="D", title="STR Amelogenin and the Y call, by source", rows="", cols="", rlev=[], clev=[],
         unit="lines", table=True),
]


def texts(cfg, d):
    k, n = int(d.same.sum()), len(d)
    if cfg["key"] == "A":
        sub = f"{k} of {n} lines show the same pattern"
    elif cfg["key"] == "B":
        conf = d[d.note == "confident in both"]
        sub = f"{k} of {n} shared lines; {int(conf.same.sum())} of {len(conf)} with confident calls in both"
    elif cfg["key"] == "C":
        n_lines = d.accession.nunique()
        n_exp = d.name.str.split(" in ").str[1].nunique()
        sub = f"{k} of {n} calls show the same pattern ({n_lines} lines, {n_exp} experiments)"
    else:
        sub = "AMELY (on Yp) in public STR profiles, against each source's Y call"
    off = d[~d.same].sort_values("name", key=lambda s: s.str.lower())
    def listing(names):
        names = [n.replace(" ", "\u00a0") for n in names]          # a line's name is never broken across lines
        return names[0] if len(names) == 1 else ", ".join(names[:-1]) + " and " + names[-1]
    if cfg["key"] == "B":
        note = ""
    elif cfg["key"] == "D":
        kinds = {("Y score high", "X only"): "no AMELY, Y expressed", ("Y score low", "X,Y"): "AMELY, Y score low"}
        parts = [f"{n.replace(' ', chr(160))} ({kinds[(a, b)]})" for n, a, b in zip(off.name, off.row_call, off.col_call)]
        note = "Tahoe-100M lines that differ: " + (" and ".join(parts) if parts else "none")
    else:
        note = "Different: " + listing(off.name.tolist()) if len(off) else "None different"
    return sub, note


def draw_table(fig, ax, counts, max_w, gap=0.12, row_h=0.19, head_h=0.34):
    """Panel D: one row per source, the numbers right-aligned under two-line headers. The axes spans the table in
    inches (x) and rows (y); returns its size in inches."""
    rend = fig.canvas.get_renderer()

    def wid(s_, **kw):
        t = fig.text(0, 0, s_, **kw)
        w_ = t.get_window_extent(rend).width / fig.dpi
        t.remove()
        return w_
    names = [nm for _, nm in STR_SOURCES]
    vals = {key: [f"{counts[key][c]:,}" for c, _ in STR_COLS] for key, _ in STR_SOURCES}
    name_w = max(wid(nm, fontsize=ROW_PT) for nm in names)
    col_w = [max([wid(h, fontsize=NOTE_PT)] + [wid(vals[key][j], fontsize=ROW_PT) for key, _ in STR_SOURCES])
             for j, (_, h) in enumerate(STR_COLS)]
    rights, x = [], name_w
    for w_ in col_w:
        x += gap + w_
        rights.append(x)
    width = rights[-1]
    assert width <= max_w, f"table {width:.2f} in wider than its cell {max_w:.2f} in"
    n = len(STR_SOURCES)
    height = head_h + n * row_h
    ax.set_xlim(0, width)
    ax.set_ylim(height, 0)
    ax.axis("off")
    out = dict(width=width, height=height, cells=[])
    for (c, h), xr in zip(STR_COLS, rights):
        ax.text(xr, head_h - 0.06, h, ha="right", va="bottom", fontsize=NOTE_PT, color=S.TEXT2, linespacing=1.1,
                multialignment="right", gid=f"D:head:{c}")
    ax.plot([0, width], [head_h - 0.02, head_h - 0.02], color=S.GRID, linewidth=0.6, gid="D:rule")
    for i, (key, nm) in enumerate(STR_SOURCES):
        y = head_h + (i + 0.5) * row_h
        ax.text(0, y, nm, ha="left", va="center", fontsize=ROW_PT, color=S.TEXT, gid=f"D:row:{key}")
        for (c, _), xr, v in zip(STR_COLS, rights, vals[key]):
            ax.text(xr, y, v, ha="right", va="center", fontsize=ROW_PT, color=S.TEXT, gid=f"D:cell:{key}:{c}")
            out["cells"].append(dict(source=key, column=c, value=counts[key][c], text=v, line=counts[key]["line"]))
    return out


def build(inputs):
    src = load(inputs)
    lines = src["lines"]
    fig = plt.figure(figsize=(S.FULL, FIG_H))
    info = {"src": src, "panels": {}}
    axes = {}
    for cfg in PANELS:
        p = cfg["key"]
        d = lines[lines.panel == p]
        ax = fig.add_axes([0.1, 0.1, 0.2, 0.2])
        axes[p] = ax
        if cfg.get("table"):
            sub, note = texts(cfg, d)
            info["panels"][p] = dict(ax=ax, tiles=[], sub=sub, note=note, cfg=cfg)
            continue
        nr, nc = len(cfg["rlev"]), len(cfg["clev"])
        ax.set_xlim(-0.5, nc - 0.5)
        ax.set_ylim(nr - 0.5, -0.5)
        tiles = []
        w = 1 - GAP_FRAC
        for i, rl in enumerate(cfg["rlev"]):
            for j, cl in enumerate(cfg["clev"]):
                n = int(((d.row_call == rl) & (d.col_call == cl)).sum())
                same = (i == j)
                face = S.BAND if n == 0 else (S.BLUE if same else S.ORANGE)
                ax.add_patch(Rectangle((j - w / 2, i - w / 2), w, w, facecolor=face, edgecolor="none", linewidth=0,
                                       zorder=1, gid=f"{p}:tile:{i}:{j}"))
                if n:
                    ax.text(j, i, f"{n}", ha="center", va="center", fontsize=7.5, color="white", zorder=2,
                            gid=f"{p}:count:{i}:{j}")
                tiles.append(dict(row=rl, col=cl, n=n, same=same))
        ax.set_xticks(range(nc), [TICK.get(c, STR_TICK.get(c, c)) for c in cfg["clev"]])
        ax.set_yticks(range(nr), cfg["rlev"])
        ax.xaxis.tick_top()
        ax.tick_params(axis="both", length=0, pad=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        ax.set_xlabel(cfg["cols"], labelpad=4)
        ax.xaxis.set_label_position("top")
        ax.set_ylabel(cfg["rows"], labelpad=4)
        sub, note = texts(cfg, d)
        info["panels"][p] = dict(ax=ax, tiles=tiles, sub=sub, note=note, cfg=cfg)
    W, H = S.FULL, FIG_H
    pad = plt.rcParams["savefig.pad_inches"]
    table = draw_table(fig, axes["D"], src["str_counts"], (W - 2 * pad) / 2)
    info["table"] = table
    # titles and notes are placed in inches from the page, beside each panel's tiles
    heads, notes = {}, {}
    for cfg in PANELS:
        p = cfg["key"]
        hax = fig.add_axes([0, 0, 0.01, 0.001])
        hax.axis("off")
        S.panel(hax, p, cfg["title"], sub=info["panels"][p]["sub"], x=0)
        hax.texts[-2].set_gid(f"{p}:title")
        hax.texts[-1].set_gid(f"{p}:subtitle")
        heads[p] = hax
        if cfg.get("table"):
            room = (W - 2 * pad) / 2 - 0.10
        else:
            room = (W - 2 * pad) / 2 - ROW_LAB_W - len(cfg["clev"]) * PITCH - NOTE_GAP - 0.10
        wrapped = V.wrap_to_width(fig, info["panels"][p]["note"], room, fontsize=NOTE_PT)
        notes[p] = fig.text(0, 0, wrapped, fontsize=NOTE_PT, color=S.TEXT2, ha="left", va="top", linespacing=1.25,
                            gid=f"{p}:note")
    leg_handles = [Rectangle((0, 0), 1, 1, facecolor=S.BLUE, edgecolor="none"),
                   Rectangle((0, 0), 1, 1, facecolor=S.ORANGE, edgecolor="none"),
                   Rectangle((0, 0), 1, 1, facecolor=S.BAND, edgecolor="none")]
    legend = fig.legend(leg_handles, ["same call", "different call", "no line"], loc="upper right",
                        bbox_to_anchor=(1, 1), ncol=3, frameon=False, handlelength=1.0,
                        handleheight=1.0, handletextpad=0.4, columnspacing=1.2, borderpad=0, borderaxespad=0,
                        fontsize=7.5, labelcolor=S.TEXT)
    legend.set_gid("legend")
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()

    def place(m):
        left, right, bottom, top = m
        legend.set_bbox_to_anchor(((W - right) / W, (H - top) / H), transform=fig.transFigure)
        col_w = (W - left - right) / 2
        row_h = (H - top - LEG_H - bottom) / 2
        for k, cfg in enumerate(PANELS):
            p = cfg["key"]
            r, cidx = divmod(k, 2)
            x0 = left + cidx * col_w
            y_top = H - top - LEG_H - r * row_h
            nr, nc = len(cfg["rlev"]), len(cfg["clev"])
            # header (letter, title, subtitle) at the cell's top, then the column labels, then the tiles
            heads[p].set_position([x0 / W, (y_top - HEAD_H) / H, 0.01 / W, 0.001 / H])
            if cfg.get("table"):
                t_top = y_top - HEAD_H - 0.12
                axes[p].set_position([x0 / W, (t_top - table["height"]) / H, table["width"] / W,
                                      table["height"] / H])
                notes[p].set_position((x0 / W, (t_top - table["height"] - 0.10) / H))
                continue
            ax_left = x0 + ROW_LAB_W
            ax_top = y_top - HEAD_H - COL_LAB_H
            axes[p].set_position([ax_left / W, (ax_top - nr * PITCH) / H, nc * PITCH / W, nr * PITCH / H])
            notes[p].set_position(((ax_left + nc * PITCH + NOTE_GAP) / W, ax_top / H))

    V.fit_page(fig, place, [pad, pad, pad, pad])
    V.size_frame(fig)
    info.update(axes=axes, heads=heads, notes=notes, legend=legend)
    return fig, info
