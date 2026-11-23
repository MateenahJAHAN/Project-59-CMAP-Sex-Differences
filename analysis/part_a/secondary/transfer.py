"""Compare the 41 P1 factorial cell lines with independent RNA measurements and STR records."""
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

WORK = Path(sys.argv[1]).resolve()
RES = WORK / "results"
ADIR = RES / "transfer"
GSETS = RES / "G/G_sets_lines.tsv"
GCOUNTS = RES / "G/G_sets_counts.tsv"
LTAB = WORK / "model_inputs/taskD/D_line_table.tsv"
STATUS = "post hoc (added at revision, after the Part B unblinding)"
YE, YB, XH, XB = "Y expressed", "Y below cut-off", "XIST high", "XIST below cut-off"
GROUPS = ["Y+ only", "X+ only", "neither"]


def read_tsv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(name, rows, cols=None):
    cols = cols or list(rows[0].keys())
    with open(ADIR / name, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {name}: {len(rows)} rows")


checks = []


def check(name, ok, detail="", kind="pipeline"):
    checks.append(dict(check=name, kind=kind, passed=str(bool(ok)), detail=detail))


DEV2 = " [added after the first results, PLAN.md Deviation 2]"


def group(y, x):
    if y in ("mixed",) or x in ("mixed",):
        return "mixed"
    if y not in (YE, YB) or x not in (XH, XB):
        return ""
    return {(YE, XB): "Y+ only", (YB, XH): "X+ only", (YB, XB): "neither", (YE, XH): "both"}[(y, x)]


def group_flags(fy, fx):
    return {(1, 0): "Y+ only", (0, 1): "X+ only", (0, 0): "neither", (1, 1): "both"}[(fy, fx)]


def taskc_group(s):
    return {"XY": "Y+ only", "XX-like": "X+ only", "Y-lost": "neither", "XIST-lost": "neither",
            "ambiguous": "undetermined"}.get(s, "")


# ------------------------------------------------------------------ inputs and their hashes
sets = [r for r in read_tsv(GSETS) if r["panel"] == "P1" and r["family"] == "factorial"]
assert len(sets) == 41 and len({r["cell_iname"] for r in sets}) == 41
lt = {r["cell_iname"]: r for r in read_tsv(LTAB)}
yc = {r["accession"]: r for r in read_tsv(ADIR / "y_calls_by_line.tsv")}
amel = {r["accession"]: r for r in read_tsv(RES / "str/amelogenin_by_record.tsv")}

rows = []
for s in sets:
    L = lt[s["cell_iname"]]
    assert L["lineage"] == s["lineage"], s
    a = L["cellosaurus_ac"]
    fy, fx = int(float(L["fact_Y"])), int(float(L["fact_X"]))
    assert L["taskA_Yhi"] in ("True", "False") and L["taskA_Xhi"] in ("True", "False")
    assert (fy == 1) == (L["taskA_Yhi"] == "True") and (fx == 1) == (L["taskA_Xhi"] == "True"), s
    assert L["elig_factorial"] == "True" and L["source_conflict"] in ("False", "True")
    Y = yc[a]
    pb = group_flags(fy, fx)
    tahoe = group(Y["tahoe_Y"], Y["tahoe_X"])
    em = group(Y["em_Y"], Y["em_X"])
    sp = group(Y["sciplex_Y"], Y["sciplex_X"])
    ccle = group(Y["taskA_ccle_Y"], Y["taskA_ccle_X"])
    if Y["taskA_sanger_Y"]:
        cmp_, cmp_src = group(Y["taskA_sanger_Y"], Y["taskA_sanger_X"]), "Sanger rows"
    elif Y["taskA_cmpbroad_Y"]:
        cmp_, cmp_src = group(Y["taskA_cmpbroad_Y"], Y["taskA_cmpbroad_X"]), "Broad rows"
    else:
        cmp_, cmp_src = "", ""
    tc = taskc_group(L["taskC_state"])
    assert L["taskC_state"] == Y["taskC_state"], (s, L["taskC_state"], Y["taskC_state"])
    tcmp = group(Y["taskC_cmp_Y"], Y["taskC_cmp_X"])
    src = L["taskA_y_score_source"]
    # the Part B call must equal task A's call from the source it names
    own = {"ccle": ccle, "sanger": cmp_ if cmp_src == "Sanger rows" else None,
           "cmpbroad": cmp_ if cmp_src == "Broad rows" else None}.get(src)
    acls = amel[a]["amelogenin_class"]
    rows.append(dict(line=s["cell_iname"], accession=a, cellosaurus_name=amel[a]["record_id"], lineage=s["lineage"],
                     group=pb, fact_Y=fy, fact_X=fx, partB_call_source=src, partB_call_equals_own_source=str(own == pb),
                     taskA_state=L["taskA_state"], taskA_confidence=L["taskA_confidence"],
                     tahoe_group=tahoe, tahoe_state=Y["tahoe_state"], tahoe_confidence=Y["tahoe_confidence"],
                     emeraldbay_group=em, emeraldbay_state=Y["em_state"], emeraldbay_confidence=Y["em_confidence"],
                     sciplex_group=sp, sciplex_experiments=Y["sciplex_experiments"],
                     ccle_group=ccle, cmp_group=cmp_, cmp_rows=cmp_src,
                     taskC_state=L["taskC_state"], taskC_group=tc, source_conflict=L["source_conflict"],
                     taskC_state_source=Y["taskC_state_source"], taskC_cmp_state=Y["taskC_cmp_state"],
                     taskC_cmp_group=tcmp,
                     amelogenin_class=acls, amelogenin_results=amel[a]["amelogenin_groups"],
                     cmap_donor_sex=L["cmap_donor_sex"], label_cellosaurus=Y["label"]))
check("every Part B call equals task A's call from the source it names (taskA_y_score_source)",
      all(r["partB_call_equals_own_source"] == "True" for r in rows),
      str([r["line"] for r in rows if r["partB_call_equals_own_source"] != "True"]))
check("no Part B line is in both groups (Y+ and X+)", all(r["group"] != "both" for r in rows))
check("the line table's taskC_state equals task C's expressed_state for every line", True,
      "asserted row by row above")

# ------------------------------------------------------------------ flags (PLAN.md, section 2.5)
for r in rows:
    reasons, notes, alt = [], [], set()
    for nm, col in (("tahoe", "tahoe_group"), ("emeraldbay", "emeraldbay_group"), ("sciplex", "sciplex_group")):
        g = r[col]
        if g == "mixed":
            reasons.append(f"{nm}_mixed")
        elif g and g != r["group"]:
            reasons.append(f"atlas_{nm}_differs")
            alt.add(g)
    for nm, conf in (("tahoe", r["tahoe_confidence"]), ("emeraldbay", r["emeraldbay_confidence"])):
        if conf == "low":
            notes.append(f"atlas_{nm}_low_confidence")
    if r["ccle_group"] and r["cmp_group"] and r["ccle_group"] != r["cmp_group"]:
        reasons.append("ccle_cmp_differ")
        alt |= {g for g in (r["ccle_group"], r["cmp_group"]) if g != r["group"]}
    if r["taskC_group"] == "undetermined":
        reasons.append("taskC_ambiguous")
    elif r["taskC_group"] and r["taskC_group"] != r["group"]:
        reasons.append("taskC_differs")
        alt.add(r["taskC_group"])
    if r["amelogenin_class"] == "X only" and r["group"] == "Y+ only":
        reasons.append("str_X_only_but_Y+")
        alt.add("X+ only" if r["fact_X"] == 1 else "neither")
    if r["taskA_confidence"] != "high":
        reasons.append(f"taskA_confidence_{r['taskA_confidence']}")
    if r["taskC_cmp_group"] == "mixed" or (r["taskC_cmp_group"] and r["taskC_cmp_group"] != r["group"]):
        notes.append("taskC_cmp_differs")  # Deviation 2: reported, does not set the flag
    if r["amelogenin_class"] == "X,Y" and r["group"] != "Y+ only":
        notes.append("str_XY_but_not_Y+")
    if r["amelogenin_class"] == "mixed across sources":
        notes.append("str_mixed")
    if r["amelogenin_class"] in ("no STR Amelogenin", "no record"):
        notes.append("str_none")
    r["flag"] = "yes" if reasons else "no"
    r["flag_reasons"] = ";".join(reasons) if reasons else "none"
    r["note_reasons"] = ";".join(notes) if notes else "none"
    r["alt_groups"] = ";".join(sorted(alt)) if alt else ""
    avail = ["task A " + {"ccle": "CCLE", "sanger": "CMP (Sanger rows)", "cmpbroad": "CMP (Broad rows)"}[r["partB_call_source"]]
             + " (Part B call)"]
    if r["ccle_group"] and r["partB_call_source"] != "ccle":
        avail.append("task A CCLE")
    if r["cmp_group"] and not (r["partB_call_source"] in ("sanger", "cmpbroad")):
        avail.append(f"task A CMP ({r['cmp_rows']})")
    if r["taskC_state"] and r["taskC_state"] != "not_measured":
        avail.append("task C")
    for nm, col in (("Tahoe-100M", "tahoe_group"), ("EmeraldBay", "emeraldbay_group"), ("sci-Plex", "sciplex_group")):
        if r[col]:
            avail.append(nm)
    if r["amelogenin_class"] not in ("no STR Amelogenin", "no record"):
        avail.append("STR Amelogenin")
    r["sources_available"] = "; ".join(avail)

    def show(g):
        return g if g else "no call"
    r["calls_by_source"] = " | ".join([
        f"Part B (task A, {r['partB_call_source']}; {r['taskA_state']}, {r['taskA_confidence']}): {r['group']}",
        f"CCLE: {show(r['ccle_group'])}", f"CMP{' (' + r['cmp_rows'] + ')' if r['cmp_rows'] else ''}: {show(r['cmp_group'])}",
        f"task C: {r['taskC_state'] or 'no row'}" + (f" (from {r['taskC_state_source']})" if r["taskC_state_source"] else ""),
        f"task C CMP: {show(r['taskC_cmp_group'])}",
        f"Tahoe-100M: {show(r['tahoe_group'])}" + (f" ({r['tahoe_state']}, {r['tahoe_confidence']})" if r['tahoe_group'] else ""),
        f"EmeraldBay: {show(r['emeraldbay_group'])}" + (f" ({r['emeraldbay_state']}, {r['emeraldbay_confidence']})" if r['emeraldbay_group'] else ""),
        f"sci-Plex: {show(r['sciplex_group'])}",
        f"STR Amelogenin: {r['amelogenin_class']}"])
    r["status"] = STATUS

COLS = ["line", "accession", "lineage", "group", "sources_available", "calls_by_source", "flag", "flag_reasons",
        "note_reasons", "alt_groups", "cellosaurus_name", "fact_Y", "fact_X", "partB_call_source", "taskA_state",
        "taskA_confidence", "ccle_group", "cmp_group", "cmp_rows", "taskC_state", "taskC_state_source", "taskC_group",
        "taskC_cmp_state", "taskC_cmp_group", "source_conflict",
        "tahoe_group", "tahoe_state", "tahoe_confidence", "emeraldbay_group", "emeraldbay_state",
        "emeraldbay_confidence", "sciplex_group", "sciplex_experiments", "amelogenin_class", "amelogenin_results",
        "label_cellosaurus", "cmap_donor_sex", "status"]
rows.sort(key=lambda r: (r["lineage"], GROUPS.index(r["group"]), r["line"]))
write_tsv("partB_lines_transfer.tsv", rows, COLS)

# ------------------------------------------------------------------ group sizes, by lineage, checked against G
gs = []
by_lin = defaultdict(Counter)
for r in rows:
    by_lin[r["lineage"]][r["group"]] += 1
gc = {r["lineage"]: r for r in read_tsv(GCOUNTS) if r["panel"] == "P1" and r["family"] == "factorial"}
ok = set(gc) == set(by_lin)
for lin in sorted(by_lin):
    c = by_lin[lin]
    g = gc.get(lin, {})
    same = (int(g.get("n_lines", -1)) == sum(c.values()) and int(g.get("n_g1", -1)) == c["Y+ only"]
            and int(g.get("n_g2", -1)) == c["X+ only"] and int(g.get("n_g3", -1)) == c["neither"]
            and int(g.get("n_g4", -1)) == 0)
    ok &= same
    gs.append(dict(lineage=lin, n_lines=sum(c.values()), n_Y_only=c["Y+ only"], n_X_only=c["X+ only"],
                   n_neither=c["neither"], equals_G_sets_counts=str(same), status=STATUS))
tot = Counter(r["group"] for r in rows)
gs.append(dict(lineage="all", n_lines=len(rows), n_Y_only=tot["Y+ only"], n_X_only=tot["X+ only"],
               n_neither=tot["neither"], equals_G_sets_counts=str(ok), status=STATUS))
check("group sizes by lineage equal G_sets_counts.tsv (P1 factorial: g1 = Y+ only, g2 = X+ only, g3 = neither, g4 = 0)",
      ok)
write_tsv("partB_group_sizes.tsv", gs)

# ------------------------------------------------------------------ agreement
ag = []


def yx(g):
    return {"Y+ only": (1, 0), "X+ only": (0, 1), "neither": (0, 0), "both": (1, 1)}.get(g)


for nm, col in (("Tahoe-100M", "tahoe_group"), ("EmeraldBay", "emeraldbay_group"), ("sci-Plex", "sciplex_group")):
    have = [r for r in rows if r[col] and r[col] != "mixed"]
    ag.append(dict(comparison=f"{nm} vs Part B call", n_lines=len(rows), n_with_both=len(have),
                   n_same_group=sum(1 for r in have if r[col] == r["group"]),
                   n_same_Y=sum(1 for r in have if yx(r[col])[0] == r["fact_Y"]),
                   n_same_X=sum(1 for r in have if yx(r[col])[1] == r["fact_X"]),
                   differing=", ".join(f"{r['line']} ({r[col]} vs {r['group']})" for r in have if r[col] != r["group"]),
                   status=STATUS))
anyat = [r for r in rows if any(r[c] for c in ("tahoe_group", "emeraldbay_group", "sciplex_group"))]
ag.append(dict(comparison="any atlas vs Part B call (every atlas call the same group)", n_lines=len(rows),
               n_with_both=len(anyat),
               n_same_group=sum(1 for r in anyat if all(r[c] == r["group"] for c in ("tahoe_group", "emeraldbay_group",
                                                                                     "sciplex_group") if r[c])),
               n_same_Y="", n_same_X="",
               differing=", ".join(r["line"] for r in anyat if any(r[c] and r[c] != r["group"] for c in (
                   "tahoe_group", "emeraldbay_group", "sciplex_group"))), status=STATUS))
both = [r for r in rows if r["ccle_group"] and r["cmp_group"]]
ag.append(dict(comparison="CCLE vs CMP (task A sources)", n_lines=len(rows), n_with_both=len(both),
               n_same_group=sum(1 for r in both if r["ccle_group"] == r["cmp_group"]),
               n_same_Y=sum(1 for r in both if yx(r["ccle_group"])[0] == yx(r["cmp_group"])[0]),
               n_same_X=sum(1 for r in both if yx(r["ccle_group"])[1] == yx(r["cmp_group"])[1]),
               differing=", ".join(f"{r['line']} (CCLE {r['ccle_group']}, CMP {r['cmp_group']})" for r in both
                                   if r["ccle_group"] != r["cmp_group"]), status=STATUS))
tcd = [r for r in rows if r["taskC_group"] in GROUPS + ["both"]]
ag.append(dict(comparison="task C state vs Part B call", n_lines=len(rows), n_with_both=len(tcd),
               n_same_group=sum(1 for r in tcd if r["taskC_group"] == r["group"]), n_same_Y="", n_same_X="",
               differing=", ".join(f"{r['line']} (task C {r['taskC_state']})" for r in rows
                                   if r["taskC_group"] and r["taskC_group"] != r["group"]), status=STATUS))
tcm = [r for r in rows if r["taskC_cmp_group"] and r["taskC_cmp_group"] != "mixed"]
ag.append(dict(comparison="CMP (task C cmp_state) vs Part B call" + DEV2, n_lines=len(rows), n_with_both=len(tcm),
               n_same_group=sum(1 for r in tcm if r["taskC_cmp_group"] == r["group"]),
               n_same_Y=sum(1 for r in tcm if yx(r["taskC_cmp_group"])[0] == r["fact_Y"]),
               n_same_X=sum(1 for r in tcm if yx(r["taskC_cmp_group"])[1] == r["fact_X"]),
               differing=", ".join(f"{r['line']} (task C CMP {r['taskC_cmp_state']}, Part B {r['group']})" for r in rows
                                   if r["taskC_cmp_group"] and r["taskC_cmp_group"] != r["group"]), status=STATUS))
for srcname in ("CMP", "CCLE"):
    sub = [r for r in tcd if r["taskC_state_source"] == srcname]
    ag.append(dict(comparison=f"task C state vs Part B call, task C state from {srcname}" + DEV2, n_lines=len(rows),
                   n_with_both=len(sub), n_same_group=sum(1 for r in sub if r["taskC_group"] == r["group"]),
                   n_same_Y="", n_same_X="",
                   differing=", ".join(r["line"] for r in sub if r["taskC_group"] != r["group"]), status=STATUS))
check("every determinate task C state comes from CMP or CCLE", all(r["taskC_state_source"] in ("CMP", "CCLE")
                                                                    for r in tcd))
st = [r for r in rows if r["amelogenin_class"] in ("X only", "X,Y")]
ag.append(dict(comparison="STR Amelogenin vs Part B Y call (agree = X,Y with Y+, or X only without Y+)",
               n_lines=len(rows), n_with_both=len(st),
               n_same_group="", n_same_Y=sum(1 for r in st if (r["amelogenin_class"] == "X,Y") == (r["fact_Y"] == 1)),
               n_same_X="",
               differing=", ".join(f"{r['line']} ({r['amelogenin_class']}, {r['group']})" for r in st
                                   if (r["amelogenin_class"] == "X,Y") != (r["fact_Y"] == 1)), status=STATUS))
write_tsv("partB_agreement.tsv", ag)

# ------------------------------------------------------------------ counts used by the text
cnt = []


def add(k, v, meaning):
    cnt.append(dict(count=k, value=v, meaning=meaning, status=STATUS))


add("n_lines", len(rows), "Part B primary-panel lines")
for g, k in (("Y+ only", "n_Y_only"), ("X+ only", "n_X_only"), ("neither", "n_neither")):
    add(k, tot[g], f"lines in group {g}")
add("n_lineages", len(by_lin), "lineages")
for s in ("ccle", "sanger", "cmpbroad"):
    add(f"n_partB_call_from_{s}", sum(1 for r in rows if r["partB_call_source"] == s), f"Part B calls taken from {s}")
add("n_with_any_atlas", len(anyat), "lines with at least one atlas call")
for nm, col in (("tahoe", "tahoe_group"), ("emeraldbay", "emeraldbay_group"), ("sciplex", "sciplex_group")):
    add(f"n_with_{nm}", sum(1 for r in rows if r[col]), f"lines with a {nm} call")
add("n_with_two_atlases", sum(1 for r in rows if sum(bool(r[c]) for c in ("tahoe_group", "emeraldbay_group",
                                                                          "sciplex_group")) >= 2),
    "lines with calls from two or more atlases")
add("n_with_taskC_cmp", len(tcm), "lines with a task C CMP call (Deviation 2)")
add("n_taskC_cmp_same_group", sum(1 for r in tcm if r["taskC_cmp_group"] == r["group"]),
    "lines whose task C CMP group equals the Part B group (Deviation 2)")
add("n_taskC_state_from_cmp", sum(1 for r in tcd if r["taskC_state_source"] == "CMP"),
    "lines whose task C state comes from CMP (Deviation 2)")
add("n_taskC_state_from_ccle", sum(1 for r in tcd if r["taskC_state_source"] == "CCLE"),
    "lines whose task C state comes from CCLE (Deviation 2)")
add("n_ccle_and_taskA_cmp", len(both), "lines with both a task A CCLE and a task A CMP call")
add("n_ccle_taskA_cmp_same_group", sum(1 for r in both if r["ccle_group"] == r["cmp_group"]),
    "of those, lines with the same group")
add("n_taskC_determinate", len(tcd), "lines with a determinate task C state")
add("n_taskC_same_group", sum(1 for r in tcd if r["taskC_group"] == r["group"]),
    "of those, lines whose task C group equals the Part B group")
add("n_str_compared", len(st), "Part B lines with Amelogenin X only or X,Y")
add("n_str_agree_Y", sum(1 for r in st if (r["amelogenin_class"] == "X,Y") == (r["fact_Y"] == 1)),
    "of those, lines whose Amelogenin result agrees with the Part B Y call")
for nm, col in (("tahoe", "tahoe_group"), ("emeraldbay", "emeraldbay_group"), ("sciplex", "sciplex_group")):
    have = [r for r in rows if r[col] and r[col] != "mixed"]
    add(f"n_{nm}_compared", len(have), f"lines with a {nm} group (not mixed)")
    add(f"n_{nm}_same_group", sum(1 for r in have if r[col] == r["group"]), f"of those, the same group as Part B")
    add(f"n_{nm}_same_Y", sum(1 for r in have if yx(r[col])[0] == r["fact_Y"]), f"of those, the same Y call")
add("n_any_atlas_all_same", sum(1 for r in anyat if all(r[c] == r["group"] for c in (
    "tahoe_group", "emeraldbay_group", "sciplex_group") if r[c])),
    "lines whose every atlas call (a mixed sci-Plex call counts as different) equals the Part B group")
add("n_any_atlas_all_same_Y", sum(1 for r in anyat if all(
    (r[c] == "mixed" and r["sciplex_experiments"] and all(
        ("Y expressed" in e) == (r["fact_Y"] == 1) for e in r["sciplex_experiments"].split("; ")))
    or (r[c] != "mixed" and yx(r[c])[0] == r["fact_Y"])
    for c in ("tahoe_group", "emeraldbay_group", "sciplex_group") if r[c])),
    "lines whose every atlas Y call (each sci-Plex experiment counted) equals the Part B Y call")
add("n_flagged", sum(1 for r in rows if r["flag"] == "yes"), "flagged lines (discordance or low confidence)")
for reason in ("atlas_tahoe_differs", "atlas_emeraldbay_differs", "atlas_sciplex_differs", "ccle_cmp_differ",
               "taskC_differs", "taskC_ambiguous", "str_X_only_but_Y+", "taskA_confidence_medium",
               "taskA_confidence_low", "tahoe_mixed", "emeraldbay_mixed", "sciplex_mixed"):
    add(f"n_reason_{reason}", sum(1 for r in rows if reason in r["flag_reasons"].split(";")),
        f"lines flagged for {reason}")
for reason in ("str_XY_but_not_Y+", "str_mixed", "str_none", "atlas_tahoe_low_confidence",
               "atlas_emeraldbay_low_confidence", "taskC_cmp_differs"):
    add(f"n_note_{reason}", sum(1 for r in rows if reason in r["note_reasons"].split(";")),
        f"lines with the note {reason} (does not set the flag)")
add("n_flag_discordance", sum(1 for r in rows if any(x for x in r["flag_reasons"].split(";")
                                                      if x != "none" and not x.startswith("taskA_confidence"))),
    "lines flagged for a discordance reason")
add("n_flag_lowconf_only", sum(1 for r in rows if r["flag"] == "yes" and all(
    x.startswith("taskA_confidence") for x in r["flag_reasons"].split(";"))), "lines flagged for low confidence only")
for c in ("X only", "X,Y", "mixed across sources", "no STR Amelogenin", "no record", "other"):
    add(f"n_amelogenin_{c.replace(' ', '_').replace(',', '')}", sum(1 for r in rows if r["amelogenin_class"] == c),
        f"Part B lines with Amelogenin class {c}")
for g in GROUPS:
    for c in ("X only", "X,Y", "mixed across sources"):
        add(f"n_{g.replace('+', '').replace(' ', '_')}_amelogenin_{c.replace(' ', '_').replace(',', '')}",
            sum(1 for r in rows if r["group"] == g and r["amelogenin_class"] == c),
            f"Part B lines in group {g} with Amelogenin class {c}")
write_tsv("partB_counts.tsv", cnt)

# ------------------------------------------------------------------ supplementary-table-ready version
supp = [dict(line=r["line"], cellosaurus_name=r["cellosaurus_name"], accession=r["accession"], lineage=r["lineage"],
             group=r["group"], partB_call_source=r["partB_call_source"], taskA_confidence=r["taskA_confidence"],
             ccle=r["ccle_group"] or "no call", cmp=(r["cmp_group"] + f" ({r['cmp_rows']})") if r["cmp_group"] else "no call",
             task_C_state=(r["taskC_state"] + (f" (from {r['taskC_state_source']})" if r["taskC_state_source"] else ""))
             if r["taskC_state"] else "no row",
             task_C_cmp=r["taskC_cmp_group"] or "no call",
             tahoe_100m=(f"{r['tahoe_group']} ({r['tahoe_confidence']})" if r["tahoe_group"] else "not in atlas"),
             emeraldbay=(f"{r['emeraldbay_group']} ({r['emeraldbay_confidence']})" if r["emeraldbay_group"] else "not in atlas"),
             sci_plex=r["sciplex_group"] or "not in atlas", str_amelogenin=r["amelogenin_class"],
             flag=r["flag"], flag_reasons=r["flag_reasons"], notes=r["note_reasons"], status=STATUS) for r in rows]
write_tsv("supp_table_partB_transfer.tsv", supp)
write_tsv("checks_analysis2.tsv", checks, ["check", "kind", "passed", "detail"])
for c in checks:
    print(("PASS " if c["passed"] == "True" else "FAIL ") + c["check"] + (f" ({c['detail']})" if c["detail"] else ""))
failed = [c["check"] for c in checks if c["passed"] != "True"]
if failed:
    raise ValueError("Cross-stock checks failed: " + ", ".join(failed))

# Use the group membership and reasons derived above for the sensitivity analysis.
flags = [dict(cell_iname=r["line"], group=r["group"], flag_reasons=r["flag_reasons"], alt_groups=r["alt_groups"])
         for r in rows if r["flag"] == "yes"]
expected = read_tsv(Path(__file__).resolve().parents[3] / "verification/source_replay/power/sensitivity_lines.tsv")
if sorted(flags, key=lambda r: r["cell_iname"]) != sorted(expected, key=lambda r: r["cell_iname"]):
    raise ValueError(f"Regenerated flagged lines differ from the published sensitivity set: {flags}")
write_tsv("sensitivity_lines.tsv", flags, ["cell_iname", "group", "flag_reasons", "alt_groups"])
