"""Step 1 of task E: panel rebuild (task D rule) and well identification.

Reads only metadata (siginfo, instinfo, geneinfo, gctx id lists). Reads no expression value
and no sex-state, Y, X or donor-sex column.

Outputs
  results/E/E_panel_definitions_rebuilt.tsv  one row per panel x compound
  results/E/E_panel_signatures.tsv.gz        signature list per panel (sig_id, line, compound, dose label)
  results/E/E_wells.tsv.gz                   one row per well used (treated or DMSO), with file row index
  results/E/E_well_counts.tsv                per panel x plate counts, plus panel totals
"""
import sys, os, time
import numpy as np, pandas as pd

ROOT = sys.argv[1]
META = f"{ROOT}/data/E/meta"
IDS = f"{ROOT}/data/E/ids"
OUT = f"{ROOT}/results/E"
os.makedirs(OUT, exist_ok=True)
FORBIDDEN = {"donor_sex", "gender", "sex"}
t0 = time.time()

SIG_COLS = ["sig_id", "pert_id", "cmap_name", "cell_iname", "pert_dose", "pert_dose_unit", "pert_idose",
            "pert_time", "pert_type", "project_code", "qc_pass", "nsample", "det_plates", "distil_ids"]
assert not FORBIDDEN & set(SIG_COLS)
sig = pd.read_csv(f"{META}/siginfo_beta.txt", sep="\t", usecols=SIG_COLS, low_memory=False)
assert sig.sig_id.is_unique


def panel(d, n_lines, thr=0.8):
    cnt = d.groupby("cmap_name").cell_iname.nunique()
    keep = sorted(cnt[cnt >= thr * n_lines].index)
    return d[d.cmap_name.isin(keep)].copy(), keep


# P1: ASG, 24 h, dose labels "1.11 uM" and "0.12 uM" pooled; compounds in >= 80% of the panel's lines
a = sig[(sig.project_code == "ASG") & (sig.pert_type == "trt_cp") & (sig.pert_time == 24)
        & sig.pert_idose.isin(["1.11 uM", "0.12 uM"])]
# P2: CPC, 10 uM, 6 h
c = sig[(sig.project_code == "CPC") & (sig.pert_type == "trt_cp") & (sig.pert_time == 6)
        & (sig.pert_idose == "10 uM")]
panels = {}
for name, d in [("P1", a), ("P2", c)]:
    nl = d.cell_iname.nunique()
    d2, keep = panel(d, nl)
    d2["panel"] = name
    panels[name] = d2
    fill = d2.groupby(["cmap_name", "cell_iname"]).ngroups / (len(keep) * nl)
    print(name, "lines", nl, "compounds", len(keep), "signatures", len(d2), "fill", round(fill, 4), flush=True)
psig = pd.concat(panels.values(), ignore_index=True)
assert psig.sig_id.is_unique

# panel definitions
rows = []
for (p, cm), g in psig.groupby(["panel", "cmap_name"]):
    nl_panel = panels[p].cell_iname.nunique()
    rows.append(dict(panel=p, cmap_name=cm, n_lines=g.cell_iname.nunique(), n_lines_panel=nl_panel,
                     n_signatures=len(g), pert_ids="|".join(sorted(g.pert_id.unique())),
                     dose_labels="|".join(sorted(g.pert_idose.unique())),
                     numeric_pert_doses="|".join(str(x) for x in sorted(g.pert_dose.unique()))))
pdef = pd.DataFrame(rows)
pdef.to_csv(f"{OUT}/E_panel_definitions_rebuilt.tsv", sep="\t", index=False)
psig[["panel", "sig_id", "cell_iname", "cmap_name", "pert_id", "pert_idose", "pert_dose", "pert_time",
      "qc_pass", "nsample", "det_plates", "distil_ids"]].to_csv(f"{OUT}/E_panel_signatures.tsv.gz", sep="\t", index=False)

# treated wells: siginfo.distil_ids (pipe-separated) -> instinfo.sample_id
tw = psig[["panel", "sig_id", "cell_iname", "cmap_name", "pert_id", "pert_idose", "pert_time", "distil_ids", "nsample"]].copy()
tw["sample_id"] = tw.distil_ids.str.split("|")
n_expected = int(tw.nsample.sum())
tw = tw.explode("sample_id").drop(columns="distil_ids")
assert len(tw) == n_expected, (len(tw), n_expected)
assert tw.sample_id.is_unique, "a well is linked to more than one panel signature"

INST_COLS = ["sample_id", "det_plate", "det_well", "pert_type", "cmap_name", "pert_id", "cell_iname",
             "project_code", "pert_time", "pert_idose", "pert_dose", "qc_pass"]
assert not FORBIDDEN & set(INST_COLS)
inst = pd.read_csv(f"{META}/instinfo_beta.txt", sep="\t", usecols=INST_COLS, low_memory=False)
assert inst.sample_id.is_unique
n0 = len(tw)
tw = tw.merge(inst, on="sample_id", how="left", suffixes=("", "_inst"), validate="one_to_one")
assert len(tw) == n0
assert tw.det_plate.notna().all(), "treated well missing from instinfo"
for col in ["cell_iname", "cmap_name", "pert_id", "pert_idose"]:
    bad = (tw[col] != tw[col + "_inst"]).sum()
    assert bad == 0, (col, bad)
assert (tw.pert_time == tw.pert_time_inst).all()
assert (tw.pert_type == "trt_cp").all()
assert (tw.sample_id.str.split(":").str[0] == tw.det_plate).all()
tw["role"] = "treated"

# one line per plate?
lpp = tw.groupby("det_plate").cell_iname.nunique()
print("plates with >1 line among treated wells:", int((lpp > 1).sum()), flush=True)

# DMSO wells on the same det_plate
plates = set(tw.det_plate)
dm = inst[(inst.pert_type == "ctl_vehicle") & (inst.cmap_name == "DMSO") & inst.det_plate.isin(plates)].copy()
plate_panel = tw.groupby("det_plate").panel.agg(lambda x: "|".join(sorted(set(x))))
assert not plate_panel.str.contains(r"\|").any(), "a plate holds wells of both panels"
dm["panel"] = dm.det_plate.map(plate_panel)
plate_line = tw.groupby("det_plate").cell_iname.first()
assert (dm.cell_iname == dm.det_plate.map(plate_line)).all(), "DMSO well line differs from plate line"
dm["role"] = "dmso"
dm["sig_id"] = ""

# map to gctx row indices
for key, df in [("l3_trt", tw), ("l3_ctl", dm)]:
    ids = pd.Series(np.arange(sum(1 for _ in open(f"{IDS}/{key}_col_ids.txt"))),
                    index=[x.rstrip("\n") for x in open(f"{IDS}/{key}_col_ids.txt")])
    df["file"] = key
    df["row_index"] = df.sample_id.map(ids)
    miss = df.row_index.isna().sum()
    print(key, "wells", len(df), "missing from level-3 file", int(miss), flush=True)
    df["in_file"] = df.row_index.notna()

keep_cols = ["panel", "role", "sample_id", "det_plate", "det_well", "cell_iname", "cmap_name", "pert_id",
             "pert_idose", "pert_dose", "pert_time", "qc_pass", "sig_id", "file", "row_index", "in_file"]
wells = pd.concat([tw[keep_cols], dm[keep_cols]], ignore_index=True)
assert wells.sample_id.is_unique
wells["row_index"] = wells.row_index.astype("Int64")
wells.to_csv(f"{OUT}/E_wells.tsv.gz", sep="\t", index=False)

# counts per panel x plate
w = wells[wells.in_file]
cnt = (w.groupby(["panel", "det_plate", "cell_iname", "role"]).size().unstack("role", fill_value=0)
       .reset_index().rename(columns={"treated": "n_treated", "dmso": "n_dmso"}))
cnt["n_treated_all"] = cnt.det_plate.map(wells[wells.role == "treated"].groupby("det_plate").size()).fillna(0).astype(int)
cnt["plate_kept"] = cnt.n_dmso >= 2
tot = []
for p, g in cnt.groupby("panel"):
    tw_p = wells[(wells.panel == p) & (wells.role == "treated")]
    tot.append(dict(panel=p, det_plate="TOTAL", cell_iname=f"{g.cell_iname.nunique()} lines",
                    n_treated=int(g.n_treated.sum()), n_dmso=int(g.n_dmso.sum()),
                    n_treated_all=len(tw_p), plate_kept=f"{int(g.plate_kept.sum())} of {len(g)} plates"))
    tot.append(dict(panel=p, det_plate="TOTAL_KEPT", cell_iname=f"{g[g.plate_kept].cell_iname.nunique()} lines",
                    n_treated=int(g[g.plate_kept].n_treated.sum()), n_dmso=int(g[g.plate_kept].n_dmso.sum()),
                    n_treated_all=np.nan, plate_kept=""))
    tot.append(dict(panel=p, det_plate="TOTAL_DROPPED_lt2_DMSO", cell_iname=f"{g[~g.plate_kept].cell_iname.nunique()} lines",
                    n_treated=int(g[~g.plate_kept].n_treated.sum()), n_dmso=int(g[~g.plate_kept].n_dmso.sum()),
                    n_treated_all=np.nan, plate_kept=f"{int((~g.plate_kept).sum())} plates"))
cnt = pd.concat([pd.DataFrame(tot), cnt[["panel", "det_plate", "cell_iname", "n_treated", "n_dmso", "n_treated_all", "plate_kept"]]],
                ignore_index=True)
cnt.to_csv(f"{OUT}/E_well_counts.tsv", sep="\t", index=False)
print(pd.DataFrame(tot).to_string(), flush=True)
print(f"step1 done {time.time()-t0:.1f}s")
