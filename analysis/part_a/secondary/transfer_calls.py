"""Construct per-accession Y/XIST calls for the cross-stock comparison from fresh state tables and cutoffs."""
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

WORK = Path(sys.argv[1]).resolve()
RES = WORK / "results"
ADIR = RES / "transfer"
ADIR.mkdir(parents=True, exist_ok=True)
A_FILE = RES / "A/A_cellline_sex_states.tsv"
C_FILE = RES / "C/C_cellline_sex_labels.tsv"

YE, YB = "Y expressed", "Y below cut-off"
XH, XB = "XIST high", "XIST below cut-off"


def read_tsv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v


checks = []


def check(name, ok, detail="", kind="pipeline"):
    checks.append(dict(check=name, kind=kind, passed=str(bool(ok)), detail=detail))


def norm_label(v):
    v = (v or "").strip()
    if v in ("Male", "Female"):
        return v
    if v in ("", "Unknown", "Sex unspecified", "no_SX_line"):
        return "Unknown"
    raise ValueError(f"unexpected label {v!r}")


def implied(state, label):
    """(Y call, XIST call) that a rule-table state (task A rule table, used by the atlases) implies; None = not fixed."""
    t = {("Male", "XY"): (YE, None), ("Male", "discordant"): (YB, XH), ("Male", "Y-lost"): (YB, XB),
         ("Female", "discordant"): (YE, None), ("Female", "XX-like"): (YB, XH), ("Female", "XIST-lost"): (YB, XB),
         ("Unknown", "XY"): (YE, None), ("Unknown", "XX-like"): (YB, XH), ("Unknown", "ambiguous"): (YB, XB)}
    return t[(label, state)]


def implied_taskC(state, label):
    """Task C rule (C_plan 3b): XY = Y hi X lo; XX-like = X hi Y lo; ambiguous = both hi, or both lo for a label
    other than Male/Female; Y-lost / XIST-lost = both lo."""
    if state == "XY":
        return YE, XB
    if state == "XX-like":
        return YB, XH
    if state in ("Y-lost", "XIST-lost"):
        return YB, XB
    if state == "ambiguous":
        return (YE, XH) if label in ("Male", "Female") else (None, None)
    return None, None  # not_measured


def call(v, cut, hi, lo):
    return None if v is None else (hi if v > cut else lo)


lines = defaultdict(dict)
accs = [r["accession"] for r in read_tsv(Path(__file__).resolve().parents[1] / "references/accessions.tsv")]
labels = defaultdict(dict)   # accession -> {source: label}
names = defaultdict(dict)

# ------------------------------------------------------------------ atlases
def atlas(tag, rows, name_col="cell_name"):
    n_mismatch_y = n_mismatch_x = 0
    for r in rows:
        a = r["cellosaurus_id"]
        lab = norm_label(r["sex_label"])
        labels[a][tag] = lab
        names[a][tag] = r[name_col]
        y = call(num(r["y_score"]), num(r["y_score_cutoff"]), YE, YB)
        x = call(num(r["xist"]), num(r["xist_cutoff"]), XH, XB)
        st = r.get("state_final", r.get("state"))
        iy, ix = implied(st, lab)
        note = []
        if iy is not None and iy != y:
            n_mismatch_y += 1
            note.append(f"score gives {y}, state {st} gives {iy}; state used")
            y = iy
        if ix is not None and ix != x:
            n_mismatch_x += 1
            note.append(f"score gives {x}, state {st} gives {ix}; state used")
            x = ix
        yield a, r, y, x, "; ".join(note)
    check(f"{tag}: Y call from the score equals the call implied by the released state", n_mismatch_y == 0,
          f"{n_mismatch_y} mismatches")
    check(f"{tag}: XIST call from the score equals the call implied by the released state (where it fixes it)",
          n_mismatch_x == 0, f"{n_mismatch_x} mismatches")


s1 = read_tsv(RES / "supplement/supplementary_table_S1.tsv")
tc = {r["score"]: r["cutoff"] for r in read_tsv(RES / "tahoe/latest/partA_local_cutoffs.tsv") if r["level"] == "lineplate"}
other = read_tsv(RES / "F1/F1_cutoffs.tsv")
ec = {r["score"]: r["cutoff"] for r in other if r["dataset"] == "EmeraldBay"}
sc = {r["score"]: r["cutoff"] for r in other if r["dataset"].startswith("sci-Plex")}


def source_rows(source, cut):
    selected = []
    for r in s1:
        if r["source"] != source:
            continue
        selected.append(dict(r, experiment=source, state_final=r["state"], confidence_final=r["confidence"],
                             arm_state_final=r["arm_state"], y_score_cutoff=cut["y_score"],
                             yp_score_cutoff=cut["yp_score"], xist_cutoff=cut["XIST"]))
    return selected


tahoe = source_rows("Tahoe-100M", tc)
assert len(tahoe) == 50
for a, r, y, x, note in list(atlas("tahoe", tahoe)):
    L = lines[a]
    L.update(tahoe_state=r["state_final"], tahoe_confidence=r["confidence_final"], tahoe_arm_state=r["arm_state_final"],
             tahoe_y_score=r["y_score"], tahoe_y_cutoff=r["y_score_cutoff"], tahoe_yp_score=r["yp_score"],
             tahoe_yp_cutoff=r["yp_score_cutoff"], tahoe_xist=r["xist"], tahoe_xist_cutoff=r["xist_cutoff"],
             tahoe_n_cells_called=r["n_vehicle_cells_called"], tahoe_Y=y, tahoe_X=x, tahoe_note=note)

em = source_rows("EmeraldBay", ec)
assert len(em) == 52 and len({r["cellosaurus_id"] for r in em}) == 52
for a, r, y, x, note in list(atlas("emeraldbay", em)):
    L = lines[a]
    L.update(em_state=r["state"], em_confidence=r["confidence"], em_y_score=r["y_score"],
             em_y_cutoff=r["y_score_cutoff"], em_xist=r["xist"], em_xist_cutoff=r["xist_cutoff"],
             em_n_cells=r["n_vehicle_cells"], em_Y=y, em_X=x, em_note=note)

sp = [r for source in ("sciPlex2", "sciPlex3", "sciPlex4") for r in source_rows(source, sc)]
assert len(sp) == 6
sp_by = defaultdict(list)
for a, r, y, x, note in list(atlas("sciplex", sp)):
    sp_by[a].append((r["experiment"], r["state"], r["confidence"], y, x, note))
for a, recs in sp_by.items():
    ys, xs = {t[3] for t in recs}, {t[4] for t in recs}
    lines[a].update(
        sciplex_Y=ys.pop() if len(ys) == 1 else "mixed", sciplex_X=xs.pop() if len(xs) == 1 else "mixed",
        sciplex_experiments="; ".join(f"{e}: {st}, {cf}, {y}, {x}" for e, st, cf, y, x, _ in sorted(recs)),
        sciplex_note="; ".join(t[5] for t in recs if t[5]))

# ------------------------------------------------------------------ task A (external, LINCS2020 lines)
A = read_tsv(A_FILE)
acut = read_tsv(RES / "A/A_cutoffs.tsv")
AC = {(r["source"], r["axis"]): float(r["cutoff"]) for r in acut}
bad_flag = 0
taskA = defaultdict(list)
for r in A:
    a = r["cellosaurus_ac"]
    if not a:
        continue
    for src, ycol, xcol in (("ccle", "y_score_ccle", "XIST_ccle"), ("sanger", "y_score_sanger", "XIST_sanger"),
                            ("cmpbroad", "y_score_cmpbroad", "XIST_cmpbroad")):
        yf, xf = r[f"Yhi_{src}"], r[f"Xhi_{src}"]
        assert yf in ("True", "False", "") and xf in ("True", "False", "")
        ys, xs = num(r[ycol]), num(r[xcol])
        if yf:
            if (ys > AC[(src, "y_score")]) != (yf == "True"):
                bad_flag += 1
        if xf:
            if (xs > AC[(src, "XIST")]) != (xf == "True"):
                bad_flag += 1
    taskA[a].append(r)
check("task A: every Yhi/Xhi flag equals its score against task A's cut-off", bad_flag == 0, f"{bad_flag} mismatches")


def flag_call(v, hi, lo):
    return None if v == "" else (hi if v == "True" else lo)


for a, rs in taskA.items():
    out = {}
    for src in ("ccle", "sanger", "cmpbroad"):
        yv = {flag_call(r[f"Yhi_{src}"], YE, YB) for r in rs}
        xv = {flag_call(r[f"Xhi_{src}"], XH, XB) for r in rs}
        out[f"taskA_{src}_Y"] = yv.pop() if len(yv) == 1 else "mixed"
        out[f"taskA_{src}_X"] = xv.pop() if len(xv) == 1 else "mixed"
    lab = {norm_label(r["cellosaurus_sex"]) for r in rs}
    assert len(lab) == 1, a
    labels[a]["taskA"] = lab.pop()
    names[a]["taskA"] = rs[0]["cellosaurus_name"]
    out["taskA_state"] = "; ".join(sorted({r["state"] for r in rs}))
    out["taskA_confidence"] = "; ".join(sorted({r["confidence"] for r in rs}))
    lines[a].update({k: (v or "") for k, v in out.items()})
check("task A: rows with the same accession give the same calls",
      all("mixed" not in [lines[a].get(f"taskA_{s}_{z}") for s in ("ccle", "sanger", "cmpbroad") for z in "YX"]
          for a in taskA), "")

# ------------------------------------------------------------------ task C (external, GDSC2 and PRISM lines)
C = read_tsv(C_FILE)
ccut = read_tsv(RES / "C/C_expression_cutoffs.tsv")
CC = {(r["expr_source"], r["variable"]): float(r["cutoff_log2"]) for r in ccut}
assert len({r["RRID"] for r in C}) == len(C)
bad_state = []
for r in C:
    a = r["RRID"]
    labels[a]["taskC"] = norm_label(r["cellosaurus_sex"])
    names[a]["taskC"] = r["model_name"]
    out = {}
    for src, pre in (("CMP", "cmp"), ("CCLE", "ccle")):
        ys, xs = num(r[f"{pre}_yscore"]), num(r[f"{pre}_xist"])
        y = call(ys, CC[(src, "yscore")], YE, YB)
        x = call(xs, CC[(src, "xist")], XH, XB)
        iy, ix = implied_taskC(r[f"{pre}_state"], r["annotated_sex"])
        if r[f"{pre}_state"] == "not_measured":
            if y is not None:
                bad_state.append(f"{a} {src}: not_measured but has a score")
        else:
            if (iy is not None and iy != y) or (ix is not None and ix != x):
                bad_state.append(f"{a} {src}: state {r[pre + '_state']} vs calls {y}, {x}")
        out[f"taskC_{pre}_Y"], out[f"taskC_{pre}_X"] = y or "", x or ""
        out[f"taskC_{pre}_state"] = r[f"{pre}_state"]
    out["taskC_state"] = r["expressed_state"]
    out["taskC_state_source"] = r["expressed_state_source"]
    lines[a].update(out)
check("task C: calls from the scores agree with task C's per-source states", not bad_state,
      f"{len(bad_state)} mismatches: {bad_state[:10]}")

# ------------------------------------------------------------------ CCLE and CMP calls, and labels
for a in accs:
    L = lines[a]
    if L.get("taskC_ccle_Y"):
        L["ccle_Y"], L["ccle_X"], L["ccle_caller"] = L["taskC_ccle_Y"], L["taskC_ccle_X"], "task C"
    elif L.get("taskA_ccle_Y"):
        L["ccle_Y"], L["ccle_X"], L["ccle_caller"] = L["taskA_ccle_Y"], L["taskA_ccle_X"], "task A"
    if L.get("taskC_cmp_Y"):
        L["cmp_Y"], L["cmp_X"], L["cmp_caller"] = L["taskC_cmp_Y"], L["taskC_cmp_X"], "task C"
    elif L.get("taskA_sanger_Y"):
        L["cmp_Y"], L["cmp_X"], L["cmp_caller"] = L["taskA_sanger_Y"], L["taskA_sanger_X"], "task A (Sanger rows)"
    elif L.get("taskA_cmpbroad_Y"):
        L["cmp_Y"], L["cmp_X"], L["cmp_caller"] = L["taskA_cmpbroad_Y"], L["taskA_cmpbroad_X"], "task A (Broad rows)"
    labs = set(labels[a].values())
    L["label"] = labs.pop() if len(labs) == 1 else ("conflict" if labs else "")
    L["label_by_file"] = "; ".join(f"{k}: {v}" for k, v in sorted(labels[a].items()))
    nm = names[a]
    L["atlas_name"] = nm.get("tahoe") or nm.get("emeraldbay") or nm.get("sciplex") or ""
    L["project_name"] = L["atlas_name"] or nm.get("taskA") or nm.get("taskC") or ""

COLS = ["accession", "project_name", "atlas_name", "label", "label_by_file",
        "tahoe_state", "tahoe_confidence", "tahoe_arm_state", "tahoe_y_score", "tahoe_y_cutoff", "tahoe_yp_score",
        "tahoe_yp_cutoff", "tahoe_xist", "tahoe_xist_cutoff", "tahoe_n_cells_called", "tahoe_Y", "tahoe_X", "tahoe_note",
        "em_state", "em_confidence", "em_y_score", "em_y_cutoff", "em_xist", "em_xist_cutoff", "em_n_cells", "em_Y",
        "em_X", "em_note", "sciplex_Y", "sciplex_X", "sciplex_experiments", "sciplex_note",
        "ccle_Y", "ccle_X", "ccle_caller", "cmp_Y", "cmp_X", "cmp_caller",
        "taskA_state", "taskA_confidence", "taskA_ccle_Y", "taskA_ccle_X", "taskA_sanger_Y", "taskA_sanger_X",
        "taskA_cmpbroad_Y", "taskA_cmpbroad_X", "taskC_state", "taskC_state_source", "taskC_cmp_state", "taskC_ccle_state",
        "taskC_cmp_Y", "taskC_cmp_X", "taskC_ccle_Y", "taskC_ccle_X"]
outside = [a for a in lines if a not in set(accs)]
outside_with_call = [a for a in outside if any(lines[a].get(k) in (YE, YB) for k in
                                               ("taskA_ccle_Y", "taskA_sanger_Y", "taskA_cmpbroad_Y",
                                                "taskC_cmp_Y", "taskC_ccle_Y"))]
check("every line with a Y call from any source is in the fetched accession list", not outside_with_call,
      f"{len(outside)} task A or task C lines are outside the list, {len(outside_with_call)} of them with a Y call")
rows = []
for a in accs:
    L = lines[a]
    rows.append({c: (a if c == "accession" else (L.get(c) if L.get(c) is not None else "")) for c in COLS})
no_call = [r["accession"] for r in rows if not any(r[k] in (YE, YB, "mixed") for k in
                                                  ("tahoe_Y", "em_Y", "sciplex_Y", "ccle_Y", "cmp_Y"))]
check("every line in the accession list has a Y call from at least one source", not no_call,
      f"{len(no_call)} without: {no_call[:10]}")
with open(ADIR / "y_calls_by_line.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=COLS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
with open(ADIR / "checks_y_calls.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["check", "kind", "passed", "detail"], delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(checks)
print(f"wrote y_calls_by_line.tsv ({len(rows)} rows) and checks_y_calls.tsv")
for c in checks:
    print(("PASS " if c["passed"] == "True" else "FAIL ") + c["check"] + (f" ({c['detail']})" if c["detail"] else ""))
failed = [c["check"] for c in checks if c["passed"] != "True"]
if failed:
    raise ValueError("RNA call checks failed: " + ", ".join(failed))
