#!/usr/bin/env python3
"""Step 2 (post hoc, revision): read the STR Amelogenin result of each cached Cellosaurus record and class it
(PLAN.md, section 1.3, with Deviation 1 for the line format).

Only the cell-line block is read: from the top of the file to the first line that is exactly '//' (the API's txt output
appends the cited references after it, each ending in '//').

Classes: X only | X,Y | mixed across sources | no STR Amelogenin | no record | other.

Usage: parse_str.py <output_dir> <cellosaurus.txt> <accessions.tsv>
Writes <analysis_dir>/amelogenin_by_record.tsv
"""
import csv
import hashlib
import re
import sys
from pathlib import Path

ADIR = Path(sys.argv[1]).resolve()
RAW = ADIR / "raw"
OUT = ADIR / "amelogenin_by_record.tsv"


def read_tsv(p):
    with open(p, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def split_outside_parens(s, sep=";"):
    parts, depth, cur = [], 0, []
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == sep and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return [p.strip() for p in parts if p.strip()]


GROUP = re.compile(r"^(?P<alleles>[^()]*?)\s*(?:\((?P<sources>.*)\))?$")


def group_result(g):
    m = GROUP.match(g)
    if not m:
        return "unparsed", "", g
    alleles = [a.strip() for a in m.group("alleles").split(",") if a.strip()]
    sources = (m.group("sources") or "").strip()
    s = set(alleles)
    if not alleles or not s <= {"X", "Y"}:
        return "unparsed", sources, g
    if s == {"X"}:
        return "X", sources, g
    if s == {"X", "Y"}:
        return "X,Y", sources, g
    return "Y only", sources, g


def release_records(accessions, release):
    """Use the archived database release for each requested primary or secondary accession."""
    records = {}
    block = []
    with release.open(encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("ID   "):
                block = [line]
            elif block:
                block.append(line)
                if line.strip() == "//":
                    ids = {a.strip() for field in block if field[:5] in ("AC   ", "AS   ")
                           for a in field[5:].split(";")}
                    matches = ids & accessions
                    if matches:
                        data = "".join(block).encode()
                        for key in matches:
                            if key in records:
                                raise ValueError(f"Duplicate accession {key}")
                            records[key] = data
                    block = []
    return records


with open(sys.argv[3]) as f:
    accs = [r["accession"] for r in csv.DictReader(f, delimiter="\t")]
records = release_records(set(accs), Path(sys.argv[2]))
rows = []
for a in accs:
    row = dict(accession=a, http_status="200" if a in records else "missing", record_ac="", record_id="",
               record_secondary_acs="", accession_is_secondary="", record_sx="", str_sources="",
               n_amelogenin_lines=0, amelogenin_raw="", amelogenin_groups="", amelogenin_class="",
               problematic_comment="", record_sha256="", parse_note="")
    if a not in records:
        row["amelogenin_class"] = "no record"
        rows.append(row)
        continue
    b = records[a]
    row["record_sha256"] = hashlib.sha256(b).hexdigest()
    lines = b.decode("utf-8").splitlines()
    block = []
    for ln in lines:
        if ln.strip() == "//":
            break
        block.append(ln)
    ids = [ln[5:].strip() for ln in block if ln.startswith("ID   ")]
    acs = [ln[5:].strip() for ln in block if ln.startswith("AC   ")]
    ass = [ln[5:].strip() for ln in block if ln.startswith("AS   ")]
    sxs = [ln[5:].strip() for ln in block if ln.startswith("SX   ")]
    assert len(ids) == 1 and len(acs) == 1, (a, ids, acs)
    row["record_id"], row["record_ac"] = ids[0], acs[0]
    sec = [x.strip() for s in ass for x in s.split(";") if x.strip()]
    row["record_secondary_acs"] = "; ".join(sec)
    if acs[0] != a:
        row["accession_is_secondary"] = str(a in sec)
        row["parse_note"] += f"record AC {acs[0]} differs from the requested accession; "
    row["record_sx"] = "; ".join(sxs)
    row["str_sources"] = "; ".join(ln.split(":", 1)[1].strip() for ln in block if ln.startswith("ST   Source(s):"))
    row["problematic_comment"] = str(any(ln.startswith("CC   Problematic cell line:") for ln in block))
    amel = [ln.split(":", 1)[1].strip() for ln in block if ln.startswith("ST   Amelogenin:")]
    row["n_amelogenin_lines"] = len(amel)
    row["amelogenin_raw"] = " | ".join(amel)
    if not amel:
        row["amelogenin_class"] = "no STR Amelogenin"
        rows.append(row)
        continue
    results, gtxt = [], []
    for v in amel:
        for g in split_outside_parens(v):
            res, src, raw = group_result(g)
            results.append(res)
            gtxt.append(f"{res} [{src}]" if src else res)
    row["amelogenin_groups"] = " | ".join(gtxt)
    s = set(results)
    if s == {"X"}:
        cls = "X only"
    elif s == {"X,Y"}:
        cls = "X,Y"
    elif s == {"X", "X,Y"}:
        cls = "mixed across sources"
    else:
        cls = "other"
        row["parse_note"] += f"group results {sorted(s)}; "
    row["amelogenin_class"] = cls
    rows.append(row)

cols = list(rows[0].keys())
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=cols, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(f"wrote {OUT.name}: {len(rows)} rows")
