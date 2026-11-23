"""Shared layout and table readers for the fellowship figures."""
import io
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.transforms import Bbox

import p59style as S

# Subscripts (log$_2$) are written with mathtext; this makes mathtext use the style's one font, Arial.
plt.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
                     "mathtext.bf": "Arial:bold", "mathtext.sf": "Arial", "mathtext.cal": "Arial",
                     "mathtext.tt": "Arial", "mathtext.default": "regular"})

MINUS = "−"
MU = "µ"
LOG2 = "log$_2$"


# ------------------------------------------------------------------------------------------------ reading
def read_tsv(path, **kw):
    """A TSV as a DataFrame with _line, each row's line number in the file (header = line 1)."""
    raw = Path(path).read_bytes()
    return read_tsv_bytes(raw, str(path), **kw)


def read_tsv_bytes(raw, name, **kw):
    df = pd.read_csv(io.BytesIO(raw), sep="\t", low_memory=False, **kw)
    n_lines = raw.count(b"\n") + (0 if raw.endswith(b"\n") else 1)
    assert n_lines == len(df) + 1, f"{name}: {n_lines} lines for {len(df)} rows, so line numbers would be wrong"
    df.insert(0, "_line", np.arange(2, len(df) + 2))
    return df


def one_row(df, **match):
    sel = df
    for k, v in match.items():
        sel = sel[sel[k] == v]
    assert len(sel) == 1, (match, len(sel))
    return sel.iloc[0]


def fmt_minus(s):
    """True minus signs in a formatted number."""
    return s.replace("-", MINUS)


# ------------------------------------------------------------------------------------------------ page
SLACK = 0.012          # inches kept free inside the pad (see fit_page)


def content_pad():
    """Where drawn content starts, in inches from the page edge: the save pad plus SLACK."""
    return plt.rcParams["savefig.pad_inches"] + SLACK


def fit_page(fig, place, margins, tries=14, slack=SLACK):
    """Adjust the four outer margins (inches: left, right, bottom, top, each measured from the page edge) until the
    drawn content ends slack inches inside the pad on every side. place(margins) positions every axes. The slack
    absorbs the fraction of a point by which text runs wider in the PDF than on the Agg canvas; size_frame() then
    fixes the saved page at exactly the figure size."""
    W, H = fig.get_size_inches()
    pad = plt.rcParams["savefig.pad_inches"] + slack
    m = list(margins)
    for _ in range(tries):
        place(m)
        fig.canvas.draw()
        bb = fig.get_tightbbox(fig.canvas.get_renderer())
        d = [bb.x0 - pad, (W - pad) - bb.x1, bb.y0 - pad, (H - pad) - bb.y1]
        if max(abs(v) for v in d) < 1e-4:
            return m
        m = [m[0] - d[0], m[1] - d[1], m[2] - d[2], m[3] - d[3]]
    raise AssertionError(f"layout did not converge: {d}")


def size_frame(fig):
    """An invisible frame inset by the pad: every artist lies inside it, so the saved tight bbox is exactly the figure
    size whatever the renderer's sub-point text extents."""
    pad = plt.rcParams["savefig.pad_inches"]
    W, H = fig.get_size_inches()
    frame = plt.Rectangle((pad / W, pad / H), 1 - 2 * pad / W, 1 - 2 * pad / H, transform=fig.transFigure,
                          facecolor="none", edgecolor="none", linewidth=0)
    frame.set_gid("size_frame")
    fig.add_artist(frame)


def printed_texts(fig):
    """Every visible string the figure draws: (gid or '', text, box in inches from the lower left corner)."""
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    inch = fig.dpi_scale_trans.inverted()
    out = []
    seen = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            for loc, tk in zip(axis.get_majorticklocs(), axis.get_major_ticks()):
                seen.add(id(tk.label1))
                seen.add(id(tk.label2))
                if lo - 1e-9 <= loc <= hi + 1e-9 and tk.label1.get_visible() and tk.label1.get_text().strip():
                    out.append(("tick", tk.label1.get_text(), tk.label1.get_window_extent(rend).transformed(inch)))
            for tk in axis.get_minor_ticks():
                seen.add(id(tk.label1))
                seen.add(id(tk.label2))
    from matplotlib.text import Text
    for t in fig.findobj(Text):
        if id(t) in seen or not t.get_visible() or not t.get_text().strip():
            continue
        out.append((t.get_gid() or "", t.get_text(), t.get_window_extent(rend).transformed(inch)))
    return out


# ------------------------------------------------------------------------------------------------ direct labels
# The label search of figures_v2 (fig1_annotated_vs_measured.py, fig2_arm_scores.py), kept unchanged: labels are
# placed by p59style.place_labels, one at a time, clear of every marker, cut-off line, fixed text and earlier label.
DIRS = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]   # p59style.place_labels' directions
RADII = (8, 14, 22, 32, 44)                                                       # and its distances, in points


def marker_radius_px(fig, s, marker, edge):
    """Outer radius of a scatter marker of area s (pt^2), its white ring included, in pixels."""
    return (np.sqrt(s) / 2 * (1.0 if marker == "o" else np.sqrt(2)) + edge / 2) * fig.dpi / 72


def grid_box(bb, step=2.0):
    return [(x, y) for x in np.arange(bb.x0, bb.x1 + step, step) for y in np.arange(bb.y0, bb.y1 + step, step)]


def grid_segment(p, q, step=2.0):
    n = max(2, int(np.hypot(q[0] - p[0], q[1] - p[1]) / step) + 1)
    return list(zip(np.linspace(p[0], q[0], n), np.linspace(p[1], q[1], n)))


def box_dist(bb, x, y):
    return float(np.hypot(max(bb.x0 - x, 0, x - bb.x1), max(bb.y0 - y, 0, y - bb.y1)))


def seg_hits_box(p, q, bb, step=1.0):
    return any(bb.x0 <= x <= bb.x1 and bb.y0 <= y <= bb.y1 for x, y in grid_segment(p, q, step))


def ambiguous(bb, own, marks, leader, factor=2.0, close_pt=6.0, dpi=150):
    """True if the label box bb could be read as another point's label (see figures_v2 for the rule)."""
    ox, oy = own
    d_own = box_dist(bb, ox, oy)
    for qx, qy, rad in marks:
        if np.hypot(qx - ox, qy - oy) < 0.5:
            continue
        d_q = box_dist(bb, qx, qy)
        if leader:
            if d_q < min(d_own, close_pt * dpi / 72):
                return True
            continue
        if d_q < d_own:
            return True
        if bb.y0 - rad <= qy <= bb.y1 + rad and not bb.x0 <= qx <= bb.x1:
            gap = bb.x0 - qx if qx < bb.x0 else qx - bb.x1
            between = bb.y0 <= oy <= bb.y1 and (qx < ox < bb.x0 if qx < bb.x0 else bb.x1 < ox < qx)
            if not between and gap < factor * d_own:
                return True
        if bb.x0 - rad <= qx <= bb.x1 + rad and not bb.y0 <= qy <= bb.y1:
            gap = bb.y0 - qy if qy < bb.y0 else qy - bb.y1
            between = bb.x0 <= ox <= bb.x1 and (qy < oy < bb.y0 if qy < bb.y0 else bb.y1 < oy < qy)
            if not between and gap < factor * d_own:
                return True
    return False


def candidate_box(own, w, h, dx, dy, r, k):
    ax_, ay_ = own[0] + dx * r * k, own[1] + dy * r * k
    x0 = ax_ if dx > 0 else (ax_ - w if dx < 0 else ax_ - w / 2)
    y0 = ay_ if dy > 0 else (ay_ - h if dy < 0 else ay_ - h / 2)
    return Bbox.from_bounds(x0, y0, w, h)


def place_direct_labels(ax, items, groups, fixed_texts, edge, label_pt=7, cut_x=(), cut_y=(), extra_xy=None,
                        shared=None):
    """items: (x, y, text); groups: (xs, ys, s, marker) for every plotted series. extra_xy: more points (data
    coordinates) that no label may cover, such as points along drawn segments. shared: for each item, the points (data
    coordinates) of other lines that the same label names (points that coincide on the page); they count as the
    label's own points. Returns [(text, leader or None)]."""
    fig = ax.figure
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    T, inv, k = ax.transData, ax.transData.inverted(), fig.dpi / 72
    marks = [(px, py, marker_radius_px(fig, s, mk, edge)) for xs, ys, s, mk in groups
             for px, py in T.transform(np.c_[xs, ys])]
    ang = np.linspace(0, 2 * np.pi, 16, endpoint=False)
    base = []
    for px, py, rad in marks:
        base.append((px, py))
        base.extend(zip(px + (rad + 1.5) * np.cos(ang), py + (rad + 1.5) * np.sin(ang)))
    axbb = ax.get_window_extent(rend)
    for c in cut_x:
        X = T.transform((c, 0))[0]
        base.extend((X, y) for y in np.arange(axbb.y0, axbb.y1, 2.0))
    for c in cut_y:
        Y = T.transform((0, c))[1]
        base.extend((x, Y) for x in np.arange(axbb.x0, axbb.x1, 2.0))
    if extra_xy is not None and len(extra_xy):
        base.extend(map(tuple, T.transform(np.asarray(extra_xy))))
    boxes = [t.get_window_extent(rend) for t in fixed_texts]
    lines_px = [("x", T.transform((c, 0))[0]) for c in cut_x] + [("y", T.transform((0, c))[1]) for c in cut_y]
    leaders, placed = [], []
    for idx, (x, y, s) in enumerate(items):
        own = tuple(T.transform((x, y)))
        same = [tuple(T.transform(p)) for p in (shared[idx] if shared else [])]
        others = [(qx, qy, rad) for qx, qy, rad in marks if all(np.hypot(qx - a, qy - b) >= 0.5 for a, b in same)]
        probe = ax.annotate(s, xy=(x, y), xytext=(0, 0), textcoords="offset points", fontsize=label_pt)
        w, h = probe.get_window_extent(rend).width, probe.get_window_extent(rend).height
        probe.remove()
        obstacles = base + [p for b in boxes for p in grid_box(b)] + [p for seg in leaders for p in grid_segment(*seg)]
        chosen = None
        for r in RADII:
            for level in (True, False):
                block = [tuple(candidate_box(own, w, h, dx, dy, r, k).corners().mean(axis=0))
                         for dx, dy in DIRS if (dy == 0) != level]
                refused = []
                while chosen is None:
                    avoid = obstacles + block + refused
                    n0 = len(ax.texts)
                    (t,) = S.place_labels(ax, [x], [y], [s], fontsize=label_pt,
                                          avoid_xy=[tuple(p) for p in inv.transform(np.array(avoid))],
                                          radius=(r,) if r == RADII[0] else (1, r))
                    new = ax.texts[n0:]
                    lead = [a for a in new if a is not t]
                    bb = t.get_window_extent(rend)
                    off = t.get_position()
                    exhausted = (not (axbb.x0 <= bb.x0 and bb.x1 <= axbb.x1 and axbb.y0 <= bb.y0 and bb.y1 <= axbb.y1)
                                 or any(bb.x0 - 1 <= px <= bb.x1 + 1 and bb.y0 - 1 <= py <= bb.y1 + 1
                                        for px, py in avoid)
                                 or (abs(off[0]), abs(off[1])) not in {(r, 0), (r, r), (0, r)})
                    if exhausted:
                        for a in new:
                            a.remove()
                        break
                    seg = (own, (own[0] + off[0] * k, own[1] + off[1] * k)) if lead else None
                    bad = ambiguous(bb, own, others, leader=bool(lead), dpi=fig.dpi)
                    if seg and not bad:
                        bad = (any(np.hypot(qx - own[0], qy - own[1]) >= 0.5 and
                                   min(np.hypot(qx - a, qy - b) for a, b in grid_segment(*seg, step=1.0)) < rad + 1
                                   for qx, qy, rad in others)
                               or any(seg_hits_box(*seg, b) for b in boxes)
                               or any(min(p[i], q[i]) <= v <= max(p[i], q[i])
                                      for (p, q) in [seg] for ax_, v in lines_px for i in [0 if ax_ == "x" else 1]))
                    if not bad:
                        chosen = (t, lead[0] if lead else None, bb, seg)
                        break
                    for a in new:
                        a.remove()
                    refused.append(tuple(bb.corners().mean(axis=0)))
                if chosen:
                    break
            if chosen:
                break
        if chosen is None:
            raise AssertionError(f"no clear place for the label {s!r}")
        t, lead, bb, seg = chosen
        boxes.append(bb)
        if seg:
            leaders.append(seg)
        placed.append((t, lead))
    return placed


def panel_lines(ax, letter, title, subs, x=-0.02, pitch=9.5):
    """p59style.panel() with more than one subtitle line: the same fonts, sizes, colours and offsets, the title
    raised by one line pitch (points) for each extra line. Returns [title, sub line 1, sub line 2, ...]."""
    ax.set_title("")
    n = len(subs)
    out = [ax.annotate(f"{letter}  {title}" if title else letter, xy=(x, 1.0), xycoords="axes fraction",
                       xytext=(0, 16 + (n - 1) * pitch), textcoords="offset points", ha="left", va="bottom",
                       fontsize=8.5, fontweight="bold", color=S.TEXT)]
    for i, s in enumerate(subs):
        out.append(ax.annotate(s, xy=(x, 1.0), xycoords="axes fraction", xytext=(0, 5 + (n - 1 - i) * pitch),
                               textcoords="offset points", ha="left", va="bottom", fontsize=7.5, color=S.TEXT2))
    return out


def wrap_to_width(fig, text, width_in, **text_kw):
    """Greedy word wrap of text into lines no wider than width_in inches when drawn with text_kw (measured)."""
    rend = fig.canvas.get_renderer()

    def w(s):
        t = fig.text(0, 0, s, **text_kw)
        out = t.get_window_extent(rend).width / fig.dpi
        t.remove()
        return out
    lines, cur = [], ""
    for word in text.split(" "):
        cand = word if not cur else cur + " " + word
        if cur and w(cand) > width_in:
            lines.append(cur)
            cur = word
        else:
            cur = cand
    lines.append(cur)
    return "\n".join(lines)
