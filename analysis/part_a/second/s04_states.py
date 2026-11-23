"""Step 3: cell table checks, pseudobulk marker values, KDE cut-offs, foreign-cell removal and state calls.
Writes results/F2/F2_cell_markers_partNN.tsv.gz, F2_cutoffs.tsv, F2_states.tsv (no comparison with any
other table is made here)."""
import glob, json, os, re, sys, time
import numpy as np, pandas as pd
from scipy.stats import gaussian_kde
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import runlog, MARKERS, Y8, YP, YQ

t0 = time.time()
R, D = "results/F2", "data/F2"
ck = json.load(open(f"{D}/reduce_ckpt.json"))
LINES, PLATES = ck["lines"], ck["plates"]
NP = len(PLATES)

# ---------- cell table ----------
cells = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{D}/cells_parts/shard_*.parquet"))], ignore_index=True)
idx = pd.read_csv(f"{R}/F2_dmso_index.tsv", sep="\t")
assert len(cells) == int(idx.vehicle_rows.sum()), (len(cells), int(idx.vehicle_rows.sum()))
assert cells.barcode.is_unique
n_vehicle_rows = len(cells)
obs = pd.read_parquet(f"{D}/obs_dmso.parquet")
full = obs[obs.pass_filter == "full"][["BARCODE_SUB_LIB_ID", "cell_name", "cell_line", "tscp_count", "gene_count", "pcnt_mito", "phase"]]
a = cells[cells.in_analysis_set].merge(full, left_on="barcode", right_on="BARCODE_SUB_LIB_ID", how="left", validate="1:1", indicator=True)
assert (a._merge == "both").all() and len(a) == int(cells.in_analysis_set.sum())
assert (a.cell_line == a.cell_line_id).all()
a = a.drop(columns=["_merge", "BARCODE_SUB_LIB_ID"]).rename(columns={"cell_name": "line"})
umi_mismatch = int((a.umis != a.tscp_count).sum())
gene_mismatch = int((a.n_genes != a.gene_count).sum())
print("analysis-set cells", len(a), "UMI != tscp_count:", umi_mismatch, "n_genes != gene_count:", gene_mismatch)
mm = a[a.umis != a.tscp_count].groupby("shard").size().rename("cells_umis_ne_tscp").reset_index()
mm.to_csv(f"{R}/F2_umi_mismatch_by_shard.tsv", sep="\t", index=False)
zf = cells.groupby("shard").agg(vehicle_cells=("barcode", "size"), cells_with_zero_entry=("n_zero_entries", lambda v: int((v > 0).sum())),
                                zero_entries=("n_zero_entries", "sum"), nongene_entries=("n_nongene_entries", "sum")).reset_index()
zf.to_csv(f"{R}/F2_zero_entries_by_shard.tsv", sep="\t", index=False)
assert (cells.n_nongene_entries == 1).all(), cells.n_nongene_entries.value_counts()
print("files with any zero entry:", int((zf.zero_entries > 0).sum()), "cells with zero entry:", int(zf.cells_with_zero_entry.sum()))
a["Y8"] = a[Y8].sum(axis=1); a["Yp"] = a[YP].sum(axis=1); a["Yq"] = a[YQ].sum(axis=1)
a["plate_i"] = a.plate.map({p: i for i, p in enumerate(PLATES)})
assert a.plate_i.notna().all()
a = a.sort_values(["line", "plate_i", "barcode"]).reset_index(drop=True)
# per-cell marker table in parts under 20 MB
out_cols = ["barcode", "line", "cell_line_id", "plate", "sample", "umis", "tscp_count", "n_genes", "mito_umis", "phase"] + MARKERS
for f in glob.glob(f"{R}/F2_cell_markers_part*.tsv.gz"):
    os.remove(f)
per = 450000
for k, s in enumerate(range(0, len(a), per)):
    a.iloc[s:s + per][out_cols].to_csv(f"{R}/F2_cell_markers_part{k + 1:02d}.tsv.gz", sep="\t", index=False, compression="gzip")
for f in glob.glob(f"{R}/F2_cell_markers_part*.tsv.gz"):
    assert os.path.getsize(f) < 20e6, f

# ---------- pseudobulk checks ----------
PB = np.load(f"{D}/pb_line_plate.npy")
gm = pd.read_parquet(f"{D}/meta/gene_metadata.parquet")
sym = gm.gene_symbol.to_numpy()
mcol = {m: int(np.flatnonzero(sym == m)[0]) for m in MARKERS}
unit = a.line.map({l: i for i, l in enumerate(LINES)}).to_numpy() * NP + a.plate_i.to_numpy()
cu = pd.DataFrame({"unit": unit, "umis": a.umis, **{m: a[m] for m in MARKERS}}).groupby("unit").sum()
assert (PB.sum(axis=1)[cu.index] == cu.umis.to_numpy()).all(), "pseudobulk totals != cell sums"
for m in MARKERS:
    assert (PB[cu.index, mcol[m]] == cu[m].to_numpy()).all(), m
assert PB.sum() == a.umis.sum()


# ---------- helpers ----------
def l2cpm(c, tot):
    return np.log2(1e6 * c / tot + 1)


def unit_values(df, by):
    g = df.groupby(by)
    s = g[["umis"] + Y8 + ["XIST"]].sum()
    v = pd.DataFrame(index=s.index)
    v["n_cells"] = g.size()
    v["umis"] = s.umis
    v["xist"] = l2cpm(s.XIST, s.umis)
    for m in Y8:
        v[m + "_l2cpm"] = l2cpm(s[m], s.umis)
    v["y_score"] = v[[m + "_l2cpm" for m in Y8]].mean(axis=1)
    v["yp_score"] = v[[m + "_l2cpm" for m in YP]].mean(axis=1)
    v["yq_score"] = v[[m + "_l2cpm" for m in YQ]].mean(axis=1)
    return v


FALLBACK = {"xist": 6.08, "y_score": 3.31, "yp_score": 3.31, "yq_score": 3.31}


def kde_cut(x, name):
    x = np.asarray(x, float)
    grid = np.linspace(x.min(), x.max(), 1024)
    d = gaussian_kde(x, bw_method="scott")(grid)
    modes = [i for i in range(1024) if (i == 0 and d[0] > d[1]) or (i == 1023 and d[1023] > d[1022])
             or (0 < i < 1023 and d[i] > d[i - 1] and d[i] > d[i + 1])]
    if len(modes) < 2:
        return dict(cutoff=FALLBACK[name], mode_low=np.nan, mode_high=np.nan, n_modes=len(modes), fallback=True, n_values=len(x))
    top = sorted(sorted(modes, key=lambda i: -d[i])[:2])
    c = top[0] + int(np.argmin(d[top[0]:top[1] + 1]))
    return dict(cutoff=grid[c], mode_low=grid[top[0]], mode_high=grid[top[1]], n_modes=len(modes), fallback=False, n_values=len(x),
                density_at_cut=d[c], density_mode_low=d[top[0]], density_mode_high=d[top[1]])


def call(sex, xist, y, cx, cy):
    yh, xh = y > cy, xist > cx
    if sex == "Female" and yh:
        return "discordant"
    if yh:
        return "XY"
    if sex == "Male" and xh:
        return "discordant"
    if xh:
        return "XX-like"
    if sex == "Male":
        return "Y-lost"
    if sex == "Female":
        return "XIST-lost"
    return "ambiguous"


def arm(yp, yq, cp, cq):
    return {(True, True): "Y intact", (True, False): "Yq lost", (False, True): "Yp lost", (False, False): "Y lost"}[(yp > cp, yq > cq)]


# ---------- Cellosaurus sex ----------
tc = pd.read_csv("inputs/C_cellline_sex_labels.tsv", sep="\t")
tcr = tc[tc.RRID.notna()]
dup = tcr.RRID[tcr.RRID.duplicated()].unique()
cl = a[["line", "cell_line_id"]].drop_duplicates()
assert len(cl) == 50 and cl.line.is_unique and cl.cell_line_id.is_unique
pie = json.load(open(f"{D}/pie/descriptions.json"))
norm = lambda s: re.sub(r"[^A-Z0-9]", "", str(s).upper())
pie_rows = []
for k, txt in pie.items():
    f = {}
    for part in txt.split(";\n"):
        if ":" in part:
            kk, vv = part.split(":", 1)
            f[kk.strip()] = vv.strip().rstrip(";").strip()
    names = {norm(k), norm(f.get("Context Name", "")), norm(f.get("Recommended Name", ""))}
    names |= {norm(s) for s in f.get("Synonyms", "").split(",") if s.strip()}
    cv = re.findall(r"CVCL_[A-Z0-9]{4}", k + " " + f.get("Context Name", ""))  # accession as the context key or name only
    pie_rows.append(dict(context=k, recommended=f.get("Recommended Name"), sex=f.get("Sex"), names=names - {""}, cvcl=set(cv)))
sexrows = []
for line, cv in zip(cl.line, cl.cell_line_id):
    r = tcr[tcr.RRID == cv]
    assert len(r) <= 1 or cv in dup
    by_id = [p for p in pie_rows if cv in p["cvcl"]]
    by_name = [p for p in pie_rows if norm(line) in p["names"]]
    pm = by_id if by_id else by_name
    pm_how = "cellosaurus_id" if by_id else ("normalised_name" if by_name else "")
    assert len(pm) <= 1, (line, [p["context"] for p in pm])
    pie_sex = pm[0]["sex"] if pm else None
    if len(r):
        sx, src = r.cellosaurus_sex.iloc[0], "TaskC_cellosaurus_sex"
    elif pie_sex:
        sx, src = pie_sex, "PIE_context_Sex"
    else:
        sx, src = "Unknown", "none"
    sexrows.append(dict(line=line, cell_line_id=cv, cellosaurus_sex=sx, sex_source=src, in_task_c=len(r) > 0,
                        pie_context=pm[0]["context"] if pm else "", pie_match=pm_how, pie_sex=pie_sex or ""))
sexdf = pd.DataFrame(sexrows)
assert sexdf.cellosaurus_sex.isin(["Male", "Female", "Unknown", "Sex unspecified"]).all(), sexdf.cellosaurus_sex.unique()
SEX = dict(zip(sexdf.line, sexdf.cellosaurus_sex))

# ---------- pass 1: all cells ----------
ncell = a.groupby("line").size()
eligible = sorted(ncell[ncell >= 200].index)
lp1 = unit_values(a[a.line.isin(eligible)], ["line", "plate"])
cut1 = {s: kde_cut(lp1[s], s) for s in ["xist", "y_score", "yp_score", "yq_score"]}
ln1 = unit_values(a[a.line.isin(eligible)], "line")
yneg = set(ln1.index[ln1.y_score <= cut1["y_score"]["cutoff"]])
xneg = set(ln1.index[ln1.xist <= cut1["xist"]["cutoff"]])
a["foreign_Y"] = a.line.isin(yneg) & (a.Y8 >= 3)
a["foreign_XIST"] = a.line.isin(xneg) & (a.XIST >= 3)
a["foreign"] = a.foreign_Y | a.foreign_XIST
a[["barcode", "line", "plate", "umis", "Y8", "XIST", "foreign_Y", "foreign_XIST"]].to_parquet(f"{D}/cells_flags.parquet", index=False)
a.to_parquet(f"{D}/cells_analysis.parquet", index=False)

# ---------- pass 2: after foreign-cell removal ----------
b = a[~a.foreign & a.line.isin(eligible)]
lp = unit_values(b, ["line", "plate"])
ln = unit_values(b, "line")
cut = {s: kde_cut(lp[s], s) for s in ["xist", "y_score", "yp_score", "yq_score"]}
cutL = {s: kde_cut(ln[s], s) for s in ["xist", "y_score", "yp_score", "yq_score"]}
crow = []
for lab, cc in [("pass1_all_cells_line_x_plate", cut1), ("final_line_x_plate", cut), ("final_line_level", cutL)]:
    for s, d in cc.items():
        crow.append(dict(cutoff_set=lab, score=s, **d))
pd.DataFrame(crow).to_csv(f"{R}/F2_cutoffs.tsv", sep="\t", index=False)


def calls(cc, tag):
    cx, cy, cp, cq = cc["xist"]["cutoff"], cc["y_score"]["cutoff"], cc["yp_score"]["cutoff"], cc["yq_score"]["cutoff"]
    rows = []
    for line in eligible:
        sx = SEX[line]
        v = ln.loc[line]
        st = call(sx, v.xist, v.y_score, cx, cy)
        pv = lp.loc[line]
        pst = [call(sx, r.xist, r.y_score, cx, cy) for _, r in pv.iterrows()]
        agree = sum(p == st for p in pst)
        grid = [call(sx, v.xist, v.y_score, cx + dx, cy + dy) for dx in (-1, 0, 1) for dy in (-1, 0, 1)]
        gs = all(g == st for g in grid)
        conf = "high" if (agree == len(pst) and gs) else ("medium" if agree >= 0.8 * len(pst) else "low")
        parm = [arm(r.yp_score, r.yq_score, cp, cq) for _, r in pv.iterrows()]
        am = arm(v.yp_score, v.yq_score, cp, cq)
        rows.append({"line": line, f"state{tag}": st, f"confidence{tag}": conf, f"plates_agree{tag}": agree, f"grid_stable{tag}": gs,
                     f"grid_states{tag}": ",".join(sorted(set(grid))), f"plate_states{tag}": ",".join(f"{p}:{s}" for p, s in zip(pv.index, pst)),
                     f"arm_state{tag}": am, f"arm_plates_agree{tag}": sum(x == am for x in parm)})
    return pd.DataFrame(rows).set_index("line")


S = calls(cut, "")
SL = calls(cutL, "_linelevel")
# pass-1 calls (all cells, pass-1 cut-offs) for the record
ln_keep = ln
ln = ln1; lp_keep = lp; lp = lp1
S1 = calls(cut1, "_pass1")[["state_pass1", "confidence_pass1"]]
ln = ln_keep; lp = lp_keep
cmeta = pd.read_parquet(f"{D}/meta/cell_line_metadata.parquet")[["cell_name", "Cell_ID_DepMap", "Cell_ID_Cellosaur", "Organ"]].drop_duplicates()
assert cmeta.cell_name.is_unique
st = sexdf.set_index("line").join(cmeta.set_index("cell_name"), how="left")
assert (st.Cell_ID_Cellosaur == st.cell_line_id).all()
fc = a.groupby("line").agg(n_vehicle_cells=("barcode", "size"), n_foreign_removed=("foreign", "sum"),
                           n_foreign_Y=("foreign_Y", "sum"), n_foreign_XIST=("foreign_XIST", "sum"))
st = st.join(fc).join(ln[["n_cells", "umis", "xist", "y_score", "yp_score", "yq_score"] + [m + "_l2cpm" for m in Y8]].rename(
    columns={"n_cells": "n_cells_called", "umis": "umis_called"}))
st["n_plates"] = lp.reset_index().groupby("line").size()
st["eligible"] = st.n_vehicle_cells >= 200
st = st.join(S).join(SL).join(S1)
st["y_negative_pass1"] = st.index.isin(yneg)
st["xist_negative_pass1"] = st.index.isin(xneg)
st.index.name = "cell_name"
st = st.reset_index()
st.to_csv(f"{R}/F2_states.tsv", sep="\t", index=False)
lp.reset_index().to_csv(f"{R}/F2_line_plate_values.tsv", sep="\t", index=False)
sexdf.to_csv(f"{D}/sex_sources.tsv", sep="\t", index=False)
json.dump(dict(n_vehicle_rows=n_vehicle_rows, n_analysis=len(a), umi_mismatch=umi_mismatch, gene_mismatch=gene_mismatch,
               n_foreign=int(a.foreign.sum()), n_called=len(b), n_lines=len(ncell), n_eligible=len(eligible)),
          open(f"{D}/consort_states.json", "w"))
print(st[["cell_name", "cellosaurus_sex", "sex_source", "n_vehicle_cells", "n_foreign_removed", "xist", "y_score", "state", "confidence",
          "plates_agree", "n_plates", "grid_stable", "arm_state", "state_linelevel"]].to_string())
print(pd.DataFrame(crow)[["cutoff_set", "score", "cutoff", "mode_low", "mode_high", "n_modes", "fallback", "n_values"]].to_string())
print(st.state.value_counts().to_string())
runlog("s04_states", t0, f"cells {len(a)}; called {len(b)}")
