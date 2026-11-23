"""A549 marker expression after depth matching."""

import json

import p59style as S
import matplotlib.pyplot as plt          # noqa: E402
from matplotlib import transforms        # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

import layout as V                     # noqa: E402

LINES = [("A549", S.BLUE, "A549"), ("NCI-H460", S.GRAY, "NCI-H460, control")]
WORDS = {"ten": 10}
EST = [("a3000", "Cells with no Y UMI", "observed share among cells\nwith at least {n} UMIs"),
       ("b2000", "Negative-binomial mixture", "fitted Y-negative weight\nat {n} UMIs"),
       ("c", "Two-Poisson mixture", "fitted Y-negative\nweight")]
B_XLIM, B_XTICKS = (0.1, 0.6), [0.2, 0.4, 0.6]
FIG_H = 4.30
ROW_PT = 7.5


def load(inputs):
    cur = V.read_tsv(inputs / "part_a/tables/partA_local_depth_curve.tsv")
    curves = {}
    for name, *_ in LINES:
        c = cur[(cur.cell_name == name) & (cur.marker == "Y")].sort_values("depth_lo")
        assert len(c) > 0 and c.median_depth.is_monotonic_increasing, name
        curves[name] = c
    P = V.read_tsv(inputs / "part_a/tables/A3_pairs.tsv", dtype=str, keep_default_na=False)
    summary = json.loads((inputs / "part_a/receipts/A3_summary.json").read_text())
    tol, nrep = summary["tolerance"], summary["n_replicates"]
    A = P[(P.line == "A549") & (P.estimator.isin([e for e, *_ in EST]))]
    groups = []
    for p in A.pair.unique():
        g = A[A.pair == p]
        ds, thin = g.dataset.iloc[0], g.thinned_side.iloc[0]
        rows = {}
        for _, r in g.iterrows():
            if r.used == "True":
                t, o, d = float(r.share_tahoe), float(r.share_other), float(r.abs_diff)
                assert abs(abs(t - o) - d) < 1.5e-6 and r.agree == str(d <= tol)
                rows[r.estimator] = dict(used=True, t=t, o=o, d=d, lo=float(r.thinned_rep_min),
                                         hi=float(r.thinned_rep_max), line=int(r._line), agree=d <= tol)
            else:
                rows[r.estimator] = dict(used=False, reason=r.not_used_reason, line=int(r._line))
        key = f"{p}:A549:{ds}"
        maxd = summary["pairs"][key]["max_abs_diff"]
        assert abs(maxd - max(v["d"] for v in rows.values() if v["used"])) < 1.5e-6
        groups.append(dict(pair=int(p), ds=ds, thin=thin, n_o=int(g.n_cells_other.iloc[0]),
                           n_t=int(g.n_cells_tahoe.iloc[0]), rows=rows, maxd=maxd,
                           agree=summary["pairs"][key]["agree"]))
    assert [g["ds"] for g in groups] == ["EmeraldBay", "sciPlex2", "sciPlex3", "sciPlex4"]
    est_n = {"a3000": int(A[A.estimator == "a3000"].estimator.iloc[0][1:]),
             "b2000": int(A[A.estimator == "b2000"].estimator.iloc[0][1:])}
    return dict(curves=curves, groups=groups, tol=tol, nrep=nrep, est_n=est_n)


def build(inputs):
    src = load(inputs)
    fig = plt.figure(figsize=(S.FULL, FIG_H))
    W, H = S.FULL, FIG_H
    pad = V.content_pad()                      # text anchored at the page edge starts here
    info = {"src": src, "rows": []}
    # ---------------------------------------------------------------- A
    axA = fig.add_axes([0.1, 0.6, 0.6, 0.3])
    for name, colour, label in LINES:
        c = src["curves"][name]
        axA.plot(c.median_depth, c.homogeneous_nb_zero_frac, color=colour, linewidth=1.0, linestyle=(0, (3, 2)),
                 zorder=2, gid=f"A:expected:{name}")
        axA.plot(c.median_depth, c.zero_frac, color=colour, linewidth=1.4, marker="o", markersize=4,
                 markerfacecolor=colour, markeredgecolor="white", markeredgewidth=0.6, zorder=3,
                 gid=f"A:observed:{name}")
        last = c.iloc[-1]
        # at the line's last point; a line that ends on the axis gets its label raised above the axis
        low = last.zero_frac < 0.04
        axA.annotate(label, xy=(last.median_depth, last.zero_frac), xytext=(6, 5 if low else 0),
                     textcoords="offset points", ha="left", va="bottom" if low else "center", fontsize=7,
                     color=S.TEXT, gid=f"A:label:{name}")
    axA.set_xscale("log")
    axA.set_xlim(400, 30000)
    axA.set_ylim(0, 0.75)
    axA.xaxis.set_major_locator(FixedLocator([1000, 3000, 10000]))
    axA.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:,.0f}"))
    axA.xaxis.set_minor_locator(NullLocator())
    axA.set_yticks([0, 0.2, 0.4, 0.6], ["0", "0.2", "0.4", "0.6"])
    S.hgrid(axA)
    axA.set_xlabel("UMIs per cell (median of each depth bin)")
    axA.set_ylabel("share with no Y UMI")
    S.panel(axA, "A", "Cells without a Y UMI, by depth",
            sub="Tahoe-100M vehicle cells; solid: observed; dashed: expected if every cell expressed Y")
    axA.texts[-2].set_gid("A:title")
    axA.texts[-1].set_gid("A:subtitle")
    # ---------------------------------------------------------------- B
    groups, tol = src["groups"], src["tol"]
    ys = {g["ds"]: -i for i, g in enumerate(groups)}
    axB = {}
    for k, (est, title, sub) in enumerate(EST):
        ax = fig.add_axes([0.2 + 0.25 * k, 0.1, 0.2, 0.25])
        axB[est] = ax
        ax.set_xlim(*B_XLIM)
        ax.set_ylim(-len(groups) + 0.5, 0.5)
        for g in groups:
            y = ys[g["ds"]]
            r = g["rows"][est]
            if not r["used"]:
                ax.text(sum(B_XLIM) / 2, y, "not usable", ha="center", va="center", fontsize=7, color=S.TEXT2,
                        gid=f"B:{est}:{g['ds']}:not_usable")
                info["rows"].append(dict(panel="B", estimator=est, dataset=g["ds"], used=False, reason=r["reason"],
                                         source_line=r["line"]))
                continue
            ax.add_patch(Rectangle((r["t"] - tol, y - 0.32), 2 * tol, 0.64, facecolor=S.BAND, edgecolor="none",
                                   linewidth=0, zorder=0, gid=f"B:{est}:{g['ds']}:band"))
            thin_t = g["thin"] == "Tahoe"
            ax.plot([r["lo"], r["hi"]], [y, y], color=S.BLUE if thin_t else S.ORANGE, linewidth=0.9,
                    solid_capstyle="butt", zorder=2, gid=f"B:{est}:{g['ds']}:range")
            ax.plot([r["t"], r["t"]], [y - 0.24, y + 0.24], color=S.BLUE, linewidth=2.0, solid_capstyle="butt",
                    zorder=3, gid=f"B:{est}:{g['ds']}:tahoe")
            ax.plot([r["o"]], [y], marker="o", markersize=4.5, markerfacecolor=S.ORANGE, markeredgecolor="white",
                    markeredgewidth=0.6, linestyle="none", zorder=4, gid=f"B:{est}:{g['ds']}:other")
            info["rows"].append(dict(panel="B", estimator=est, dataset=g["ds"], used=True, tahoe=r["t"], other=r["o"],
                                     diff=r["d"], range_lo=r["lo"], range_hi=r["hi"], thinned=g["thin"],
                                     within=r["agree"], source_line=r["line"]))
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xticks(B_XTICKS, [f"{v:.1f}" for v in B_XTICKS])
        ax.set_xlabel("estimator output")
        n = src["est_n"].get(est)
        ax.annotate(title, xy=(0, 1), xycoords="axes fraction", xytext=(0, 23), textcoords="offset points",
                    ha="left", va="bottom", fontsize=7.5, fontweight="bold", color=S.TEXT, gid=f"B:{est}:head")
        ax.annotate(sub.format(n=f"{n:,}" if n else ""), xy=(0, 1), xycoords="axes fraction", xytext=(0, 4),
                    textcoords="offset points", ha="left", va="bottom", fontsize=7, color=S.TEXT2,
                    linespacing=1.1, gid=f"B:{est}:sub")
    # row labels (left of the first facet) and the largest difference (right of the last facet)
    first, last = axB[EST[0][0]], axB[EST[-1][0]]
    tr_l = transforms.blended_transform_factory(fig.transFigure, first.transData)
    tr_r = transforms.blended_transform_factory(last.transAxes, last.transData)
    for g in groups:
        fig.text(pad / W, ys[g["ds"]], f"{g['ds']}, {g['n_o']:,} cells", transform=tr_l, ha="left", va="center",
                 fontsize=ROW_PT, color=S.TEXT, gid=f"B:row:{g['ds']}")
        last.annotate(f"{g['maxd']:.3f}", xy=(1, ys[g["ds"]]), xycoords=tr_r, xytext=(16, 0),
                      textcoords="offset points", ha="left", va="center", fontsize=ROW_PT, color=S.TEXT,
                      gid=f"B:maxd:{g['ds']}")
        info["rows"].append(dict(panel="B", estimator="max", dataset=g["ds"], diff=g["maxd"], within=g["agree"]))
    last.annotate("largest\ndifference", xy=(1, 1), xycoords="axes fraction", xytext=(16, 4),
                  textcoords="offset points", ha="left", va="bottom", fontsize=7, color=S.TEXT2, linespacing=1.1,
                  gid="B:maxd_head")
    tax = fig.add_axes([0, 0, 0.01, 0.001])
    tax.axis("off")
    S.panel(tax, "B", "A549 estimator outputs after depth matching",
            sub=f"tick: Tahoe-100M ({groups[0]['n_t']:,} cells); dot: other culture; band: within {tol:g} of Tahoe; "
                f"line: range over {src['nrep']} thinning replicates", x=0)
    tax.texts[-2].set_gid("B:title")
    tax.texts[-1].set_gid("B:subtitle")
    rend = fig.canvas.get_renderer()

    def width(s, **kw):
        t = fig.text(0, 0, s, **kw)
        w_ = t.get_window_extent(rend).width / fig.dpi
        t.remove()
        return w_

    lab_w = max(width(f"{g['ds']}, {g['n_o']:,} cells", fontsize=ROW_PT) for g in groups)
    right_w = max([width(f"{g['maxd']:.3f}", fontsize=ROW_PT) for g in groups] + [width("difference", fontsize=7)])

    def place(m):
        left, right, bottom, top = m
        hA = 1.25
        axA.set_position([left / W, (H - top - hA) / H, (W - left - right - 1.05) / W, hA / H])
        b_left = pad + lab_w + 0.12
        avail = W - right - (16 / 72 + right_w) - b_left
        gap = 0.16
        wf = (avail - 2 * gap) / 3
        hB = len(groups) * 0.23
        for k, (est, *_r) in enumerate(EST):
            axB[est].set_position([(b_left + k * (wf + gap)) / W, bottom / H, wf / W, hB / H])
        tax.set_position([pad / W, (bottom + hB + 0.52) / H, 0.01 / W, 0.001 / H])

    V.fit_page(fig, place, [0.45, 0.05, 0.40, 0.35])
    V.size_frame(fig)
    info.update(axA=axA, axB=axB, tax=tax)
    return fig, info
