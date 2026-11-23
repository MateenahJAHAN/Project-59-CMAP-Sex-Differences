"""Apply the release state rule to LINCS cells using Cellosaurus annotations and the original nine-point grids."""
import re

TASK_A_SOURCES = ("ccle", "sanger", "cmpbroad")


CELLO_RULE_LABEL = {"Male": "M", "Female": "F", "Sex unspecified": "Unknown", "": "Unknown"}


RULE_LABEL_NAME = {"M": "Male", "F": "Female", "Unknown": "Unknown"}


def index(rows, *cols):
    out = {}
    for r in rows:
        k = r[cols[0]] if len(cols) == 1 else tuple(r[c] for c in cols)
        assert k not in out, f"duplicate key {k}"
        out[k] = r
    return out


def a_call(label, y_high, xist_high):
    """Task A's rule table for one RNA-seq source (src/A/02_states.py, call_state). label: M, F or Unknown."""
    if label == "M":
        return "XY" if y_high else ("discordant" if xist_high else "Y-lost")
    if label == "F":
        return "discordant" if y_high else ("XX-like" if xist_high else "XIST-lost")
    return "XY" if y_high else ("XX-like" if xist_high else "ambiguous")


A_INVERSE = {lab: {a_call(lab, y, x): (y, None if y else x) for y in (True, False) for x in (True, False)}
             for lab in RULE_LABEL_NAME}


def a_combine(per, grids):
    """Task A's combined state, confidence and note (src/A/02_states.py, main)."""
    avail = [s for s in TASK_A_SOURCES if per[s] != "no_data"]
    if not avail:
        return "no_data", "none", ""
    calls = {s: per[s] for s in avail}
    if len(set(calls.values())) > 1:
        return "ambiguous", "low", "sources disagree: " + ", ".join(f"{k}={v}" for k, v in calls.items())
    st = calls[avail[0]]
    if st == "ambiguous":
        cf = "low"
    else:
        cf = "high" if any(len(set(grids[s])) == 1 for s in avail) else "medium"
    return st, cf, (f"single source: {avail[0]}" if len(avail) == 1 else "")


def task_a_recall(D):
    """Task A's calls from its own high and low calls, with the CMap label and with the Cellosaurus label.

    Task A wrote, for each source, the base calls (Yhi_*, Xhi_*) and the states at the 9 points of its grid of
    shifted cut-offs (A_state_sensitivity_per_line.tsv). The grid states give back the calls at each point: under any
    label the three outcomes map to three different states. Checks: the base calls give task A's per-source states;
    the grid calls agree with task A's values and cut-offs at every point; and with the CMap label the combined
    state, confidence and note equal task A's frozen values on every row. Returns cell_iname -> result."""
    grid = index(D["a_grid"], "cell_iname")
    cut = {(r["source"], r["axis"]): float(r["cutoff"]) for r in D["a_cutoffs"]}
    pat = re.compile(r"^state_(ccle|sanger|cmpbroad)_dy([+-]\d+)_dx([+-]\d+)$")
    points = {}
    for col in next(iter(grid.values())):
        if col != "cell_iname":
            m = pat.match(col)
            assert m, col
            points.setdefault(m.group(1), []).append((int(m.group(2)), int(m.group(3)), col))
    assert {s: len(v) for s, v in points.items()} == {s: 9 for s in TASK_A_SOURCES}
    assert all((0, 0) in {(dy, dx) for dy, dx, _ in v} for v in points.values())
    decimals = {len(a[c].split(".")[1]) for a in D["a_states"] for c in a
                if c.split("_")[0] in {"XIST", "y"} and c.endswith(TASK_A_SOURCES) and a[c]}
    assert len(decimals) == 1, decimals
    tol = 0.5 * 10 ** -decimals.pop()  # task A rounded its values; a value this close to a cut-off is not compared
    flag = {"True": True, "False": False, "": None}
    out, n_points, n_near = {}, 0, 0
    assert len(grid) == len(D["a_states"])
    for a in D["a_states"]:
        assert a["cmap_donor_sex"] in RULE_LABEL_NAME and a["cellosaurus_sex"] in CELLO_RULE_LABEL
        g = grid[a["cell_iname"]]
        lab = {"cmap": a["cmap_donor_sex"], "cellosaurus": CELLO_RULE_LABEL[a["cellosaurus_sex"]]}
        per = {k: {} for k in lab}
        grids = {k: {} for k in lab}
        for s in TASK_A_SOURCES:
            yhi, xhi = flag[a[f"Yhi_{s}"]], flag[a[f"Xhi_{s}"]]
            assert (yhi is None) == (xhi is None)
            base = "no_data" if yhi is None else a_call(lab["cmap"], yhi, xhi)
            assert base == a[f"state_{s}"], (a["cell_iname"], s)
            for k in lab:
                per[k][s] = "no_data" if yhi is None else a_call(lab[k], yhi, xhi)
                grids[k][s] = []
            for dy, dx, col in points[s]:
                st = g[col]
                if st == "no_data":
                    assert yhi is None
                    continue
                y_g, x_g = A_INVERSE[lab["cmap"]][st]
                if (dy, dx) == (0, 0):
                    assert st == base and y_g == yhi and (y_g or x_g == xhi)
                y, x = float(a[f"y_score_{s}"]), float(a[f"XIST_{s}"])
                for v, c, call, used in ((y, cut[(s, "y_score")] + dy, y_g, True), (x, cut[(s, "XIST")] + dx, x_g, not y_g)):
                    if not used:
                        continue
                    if abs(v - c) <= tol:
                        n_near += 1
                    else:
                        assert (v >= c) == call, (a["cell_iname"], s, dy, dx)
                        n_points += 1
                for k in lab:
                    grids[k][s].append(a_call(lab[k], y_g, x_g))
        res = {k: (per[k],) + a_combine(per[k], grids[k]) for k in lab}
        assert res["cmap"][1:] == (a["state"], a["confidence"], a["notes"]), a["cell_iname"]
        out[a["cell_iname"]] = res
    D["a_recall_points"] = (n_points, n_near)
    return out
