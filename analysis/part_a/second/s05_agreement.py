"""Step 4: agreement with Task C, Task A and the pilot table; label audit; PIE mapping.
Run only after results/F2/F2_states.tsv and F2_hashes.tsv were pushed."""
import sys, time
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import runlog

t0 = time.time()
R = "results/F2"
st = pd.read_csv(f"{R}/F2_states.tsv", sep="\t")
st = st[st.eligible]
tc = pd.read_csv("inputs/C_cellline_sex_labels.tsv", sep="\t")
ta = pd.read_csv("inputs/A_cellline_sex_states.tsv", sep="\t")
pi = pd.read_csv("inputs/tahoe_pilot_states.tsv", sep="\t")
assert tc.RRID.dropna().is_unique and pi.cell_line_id.is_unique

refs = {}
c = tc[tc.RRID.notna()][["RRID", "model_name", "expressed_state", "expressed_state_source"]].rename(columns={"RRID": "cell_line_id", "expressed_state": "ref_state"})
c = c[c.ref_state.notna() & (c.ref_state != "not_measured")]
refs["TaskC_expressed_state"] = c
a = ta[ta.cellosaurus_ac.notna() & (ta.state != "no_data")][["cellosaurus_ac", "cell_iname", "state"]].rename(columns={"cellosaurus_ac": "cell_line_id", "state": "ref_state"})
dupA = a.cell_line_id[a.cell_line_id.duplicated()].unique().tolist()
a = a.groupby("cell_line_id").agg(cell_iname=("cell_iname", ",".join), ref_state=("ref_state", lambda v: v.iloc[0] if v.nunique() == 1 else "conflict")).reset_index()
refs["TaskA_state"] = a
p = pi[["cell_line_id", "cell_name", "tahoe_state", "confidence"]].rename(columns={"tahoe_state": "ref_state", "cell_name": "pilot_name", "confidence": "pilot_confidence"})
refs["pilot_tahoe_state"] = p

THR = {"TaskC_expressed_state": (0.90, 0.80)}
summ, detail = [], []
for name, ref in refs.items():
    for lev, col in [("line_x_plate_cutoffs", "state"), ("line_level_cutoffs", "state_linelevel")]:
        m = st[["cell_name", "cell_line_id", col, "confidence" if col == "state" else "confidence_linelevel"]].merge(ref, on="cell_line_id", how="inner", validate="1:1")
        m = m[m.ref_state != "conflict"]
        m["agree"] = m[col] == m.ref_state
        n, k = len(m), int(m.agree.sum())
        row = dict(comparison=name, cutoffs=lev, n_lines_both=n, n_agree=k, fraction=k / n if n else np.nan,
                   disagreeing_lines=";".join(f"{r.cell_name} ({r[col]} vs {r.ref_state})" for _, r in m[~m.agree].iterrows()))
        if name in THR:
            row["threshold_pass"] = THR[name][0]
            row["pass"] = bool(k / n >= THR[name][0])
            row["below_culture_specific_threshold_0.80"] = bool(k / n < THR[name][1])
        summ.append(row)
        for _, r in m.iterrows():
            detail.append(dict(comparison=name, cutoffs=lev, cell_name=r.cell_name, cell_line_id=r.cell_line_id, f2_state=r[col], ref_state=r.ref_state, agree=r.agree))
S = pd.DataFrame(summ)
S.to_csv(f"{R}/F2_agreement.tsv", sep="\t", index=False)
pd.DataFrame(detail).to_csv(f"{R}/F2_agreement_by_line.tsv", sep="\t", index=False)
print(S.drop(columns=["disagreeing_lines"]).to_string())
print("Task A duplicate accessions:", dupA)

# ---------- label audit ----------
usual = {"Male": "XY", "Female": "XX-like"}
la = st[["cell_name", "cell_line_id", "cellosaurus_sex", "sex_source", "pie_context", "pie_match", "pie_sex", "state", "confidence",
         "state_linelevel", "confidence_linelevel", "arm_state", "arm_state_linelevel"]].copy()
la["labelled"] = la.cellosaurus_sex.isin(["Male", "Female"])
la["usual_pattern"] = la.cellosaurus_sex.map(usual)
la["matches_label_line_x_plate"] = np.where(la.labelled, la.state == la.usual_pattern, np.nan)
la["matches_label_line_level"] = np.where(la.labelled, la.state_linelevel == la.usual_pattern, np.nan)
la["pie_sex_equals_cellosaurus_sex"] = np.where(la.pie_sex.notna() & (la.pie_sex != ""), la.pie_sex == la.cellosaurus_sex, np.nan)
la.to_csv(f"{R}/F2_label_audit.tsv", sep="\t", index=False)
lab = la[la.labelled]
ls = pd.DataFrame([
    dict(item="lines eligible", value=len(la)),
    dict(item="lines with Male or Female Cellosaurus label", value=len(lab)),
    dict(item="labelled lines without the label's usual pattern, line x plate cut-offs", value=int((lab.matches_label_line_x_plate == 0).sum())),
    dict(item="labelled lines without the label's usual pattern, line-level cut-offs", value=int((lab.matches_label_line_level == 0).sum())),
    dict(item="lines mapped to a PIE context", value=int((la.pie_context.notna() & (la.pie_context != "")).sum())),
    dict(item="lines mapped to PIE by Cellosaurus accession", value=int((la.pie_match == "cellosaurus_id").sum())),
    dict(item="lines mapped to PIE by normalised name", value=int((la.pie_match == "normalised_name").sum())),
    dict(item="mapped lines whose PIE Sex equals the Cellosaurus sex used", value=int((la.pie_sex_equals_cellosaurus_sex == 1).sum())),
    dict(item="lines whose sex came from PIE (no Task C row)", value=int((la.sex_source == "PIE_context_Sex").sum())),
] + [dict(item=f"state count, line x plate cut-offs: {k}", value=int(v)) for k, v in la.state.value_counts().items()]
  + [dict(item=f"state count, line-level cut-offs: {k}", value=int(v)) for k, v in la.state_linelevel.value_counts().items()]
  + [dict(item=f"labelled Male: arm state {k}", value=int(v)) for k, v in la[la.cellosaurus_sex == "Male"].arm_state.value_counts().items()])
ls.to_csv(f"{R}/F2_label_audit_summary.tsv", sep="\t", index=False)
print(ls.to_string())
runlog("s05_agreement", t0)

# ---------- per-plate check (positive control: markers separate states in every plate) ----------
st_all = pd.read_csv(f"{R}/F2_states.tsv", sep="\t")
rows = []
for _, r in st_all.iterrows():
    for item in r.plate_states.split(","):
        p, s = item.split(":")
        rows.append(dict(cell_name=r.cell_name, plate=p, plate_state=s, line_state=r.state))
pp = pd.DataFrame(rows)
pp["agree"] = pp.plate_state == pp.line_state
ps = pp.groupby("plate").agg(lines=("cell_name", "size"), lines_agreeing=("agree", "sum")).reset_index()
ps["plate_n"] = ps.plate.str.replace("plate", "").astype(int)
ps = ps.sort_values("plate_n")
ps.to_csv(f"{R}/F2_plate_agreement.tsv", sep="\t", index=False)
lines_all = pp.groupby("cell_name").agree.all()
pd.DataFrame([dict(item="lines whose call is the same on every plate", value=int(lines_all.sum()), of=len(lines_all)),
              dict(item="line x plate units whose plate call equals the line call", value=int(pp.agree.sum()), of=len(pp))]).to_csv(
    f"{R}/F2_plate_agreement_summary.tsv", sep="\t", index=False)
print(ps.to_string()); print(lines_all[~lines_all])
ta_dup = ta[ta.cellosaurus_ac == "CVCL_1511"][["cell_iname", "cellosaurus_ac", "state"]]
ta_dup.to_csv(f"{R}/F2_taskA_duplicate_accessions.tsv", sep="\t", index=False)
print(ta_dup)
