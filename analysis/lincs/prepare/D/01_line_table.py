"""Task D script 1: the line table. Reads no L1000 expression matrix.

Inputs: Locally rebuilt A and C state tables, CMap
cellinfo and siginfo (metadata), Cell Model Passports model list and mutation
summary, CCLE 2018 MAF, Cellosaurus, CCLE RPKM (NR3C1 row only).

Outputs (results/D/):
  D_line_table.tsv        one row per panel line, with states (frozen; for task G)
  D_line_covariates.tsv   same lines, no state / Y / X / sex column (for analysis scripts)
  D_eligible_lines.tsv    per panel and contrast, eligible lines (no group label)
  D_lineage_counts.tsv    per panel, contrast and lineage, group counts
  D_acceptance_checks.tsv the pilot acceptance checks
  D_line_table_sources.tsv SHA-256 of the two source tables
"""
import os, re, gzip, hashlib, time
import numpy as np, pandas as pd

ROOT = os.environ["P59_LINCS_WORK"]
RES = os.path.join(ROOT, "results", "D")
DAT = os.path.join(ROOT, "data", "D")
LOG = os.path.join(ROOT, "logs", "D_runlog.txt")
T0 = time.time()


def log(msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} [01_line_table +{time.time()-T0:.0f}s] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def read_states(path):
    fullpath = os.path.join(ROOT, path)
    with open(fullpath, "rb") as stream:
        fingerprint = hashlib.file_digest(stream, "sha256").hexdigest()
    return pd.read_csv(fullpath, sep="\t", low_memory=False), fingerprint


def merge_checked(left, right, on, how="left", name=""):
    n0 = len(left)
    assert right[on].dropna().is_unique, f"{name}: right key {on} not unique"
    out = left.merge(right, on=on, how=how, validate="m:1")
    assert len(out) == n0, f"{name}: row count changed {n0} -> {len(out)}"
    return out


# ---------------- panels (same rule as script 2) ----------------
si = pd.read_csv(os.path.join(DAT, "cmap/siginfo_beta.txt"), sep="\t", low_memory=False,
                 usecols=["sig_id", "pert_type", "cell_iname", "project_code", "cmap_name", "pert_itime", "pert_idose", "qc_pass"])
PANELS = {"P1": ("ASG", "24 h", ["1.11 uM", "0.12 uM"], "1.11 uM"),
          "P1b": ("ASG", "24 h", ["1.11 uM", "0.12 uM"], "0.12 uM"),
          "P2": ("CPC", "6 h", ["10 uM"], "10 uM")}
panel_sigs, panel_lines = {}, {}
for pn, (proj, tm, gl, lab) in PANELS.items():
    b = si[(si.project_code == proj) & (si.pert_type == "trt_cp") & (si.pert_itime == tm) & si.pert_idose.isin(gl)]
    lines = sorted(b.cell_iname.unique())
    cnt = b.groupby("cmap_name").cell_iname.nunique()
    keep = cnt[cnt >= 0.8 * len(lines)].index
    s = b[b.cmap_name.isin(keep) & (b.pert_idose == lab)]
    panel_sigs[pn] = s
    panel_lines[pn] = lines
pl = pd.read_csv(os.path.join(RES, "D_panel_lines.tsv"), sep="\t")
for pn in PANELS:
    assert sorted(pl[pl.panel == pn].cell_iname) == panel_lines[pn], pn
all_lines = sorted(set().union(*[set(v) for v in panel_lines.values()]))
log(f"panel lines: P1 {len(panel_lines['P1'])}, P1b {len(panel_lines['P1b'])}, P2 {len(panel_lines['P2'])}, union {len(all_lines)}")

# ---------------- CMap cellinfo ----------------
ci = pd.read_csv(os.path.join(DAT, "cmap/cellinfo_beta.txt"), sep="\t", dtype=str, keep_default_na=False)
assert ci.cell_iname.is_unique
for c in ["cellosaurus_id", "donor_sex", "cell_lineage", "ccle_name"]:
    assert c in ci.columns, c
T = pd.DataFrame({"cell_iname": all_lines})
T = merge_checked(T, ci[["cell_iname", "cell_lineage", "donor_sex", "cellosaurus_id", "ccle_name", "primary_disease"]]
                  .rename(columns={"cell_lineage": "lineage", "donor_sex": "cmap_donor_sex", "ccle_name": "cmap_ccle_name",
                                   "cellosaurus_id": "cmap_cellosaurus_id"}), "cell_iname", name="cellinfo")
assert T.lineage.notna().all()
for pn in PANELS:
    T[f"in_{pn}"] = T.cell_iname.isin(panel_lines[pn])

# ---------------- Task A ----------------
A, commitA = read_states("results/A/A_cellline_sex_states.tsv")
assert A.cell_iname.is_unique
acol = ["cell_iname", "state", "confidence", "y_score", "y_score_source", "Yhi", "Xhi", "XIST_ccle", "XIST_sanger", "XIST_cmpbroad",
        "cmap_donor_sex", "ccle_name", "sanger_model_id", "cellosaurus_ac"]
Aa = A[acol].rename(columns={c: f"taskA_{c}" for c in acol if c != "cell_iname"})
T = merge_checked(T, Aa, "cell_iname", name="taskA")
allowedA = {"XY", "XX-like", "Y-lost", "XIST-lost", "discordant", "ambiguous", "no_data"}
assert set(T.taskA_state.dropna()) <= allowedA, set(T.taskA_state.dropna())
T["taskA_state"] = T.taskA_state.fillna("no_data")
# CMap donor sex check
cmapA = T.taskA_cmap_donor_sex.fillna("NA")
mism = T[(cmapA != "NA") & (cmapA != T.cmap_donor_sex)]
log(f"Task A rows matched {T.taskA_confidence.notna().sum()}; donor-sex mismatches with cellinfo: {len(mism)}")

# ---------------- CMP model list ----------------
ml = pd.read_csv(os.path.join(DAT, "cmp/model_list_20240110.csv"), low_memory=False)
for c in ["model_id", "ploidy_wes", "ploidy_snp6", "ploidy_wgs", "msi_status", "mutational_burden", "cancer_type", "BROAD_ID", "CCLE_ID", "RRID"]:
    assert c in ml.columns, c
assert ml.model_id.is_unique
by_rrid = ml.dropna(subset=["RRID"]).groupby("RRID").model_id.apply(list).to_dict()
by_ccle = ml.dropna(subset=["CCLE_ID"]).groupby("CCLE_ID").model_id.apply(list).to_dict()
mids, keys = [], []
for _, r in T.iterrows():
    m, k = None, ""
    rr = r.cmap_cellosaurus_id
    if rr and rr in by_rrid and len(by_rrid[rr]) == 1:
        m, k = by_rrid[rr][0], "RRID"
    elif r.cmap_ccle_name and r.cmap_ccle_name in by_ccle and len(by_ccle[r.cmap_ccle_name]) == 1:
        m, k = by_ccle[r.cmap_ccle_name][0], "CCLE_ID"
    elif isinstance(r.taskA_sanger_model_id, str) and r.taskA_sanger_model_id in set(ml.model_id):
        m, k = r.taskA_sanger_model_id, "taskA_sanger_model_id"
    mids.append(m); keys.append(k)
T["model_id"] = mids; T["cmp_match_key"] = keys
# consistency with Task A's mapping where both exist
both = T.model_id.notna() & T.taskA_sanger_model_id.notna()
n_dis = int((T.loc[both, "model_id"] != T.loc[both, "taskA_sanger_model_id"]).sum())
log(f"CMP models matched {T.model_id.notna().sum()} of {len(T)}; keys {T.cmp_match_key.value_counts().to_dict()}; disagreements with Task A mapping {n_dis}")
mlc = ml[["model_id", "model_name", "RRID", "BROAD_ID", "CCLE_ID", "ploidy_wes", "ploidy_snp6", "ploidy_wgs", "msi_status",
          "mutational_burden", "cancer_type"]].rename(columns={"model_name": "cmp_model_name", "cancer_type": "cmp_cancer_type"})
dup = T.model_id.dropna()
assert dup.is_unique, dup[dup.duplicated()].tolist()
T = merge_checked(T, mlc, "model_id", name="model_list")
T["ploidy"] = T.ploidy_wes.where(T.ploidy_wes.notna(), T.ploidy_snp6.where(T.ploidy_snp6.notna(), T.ploidy_wgs))
T["ploidy_source"] = np.where(T.ploidy_wes.notna(), "wes", np.where(T.ploidy_snp6.notna(), "snp6", np.where(T.ploidy_wgs.notna(), "wgs", "")))

# ---------------- CMP mutations ----------------
mu = pd.read_csv(os.path.join(DAT, "cmp/mutations_summary_latest.csv.gz"), low_memory=False)
assert {"model_id", "gene_symbol", "coding", "protein_mutation"} <= set(mu.columns)
seq = set(mu.model_id)
tp = mu[(mu.gene_symbol == "TP53") & (mu.coding == "t")]
tp_m = set(tp.model_id)
br = mu[(mu.gene_symbol == "BRAF") & mu.protein_mutation.astype(str).str.startswith("p.V600")]
br_m = set(br.model_id)


def call(mid, mut, seqset):
    if not isinstance(mid, str) or mid not in seqset:
        return "unknown"
    return "mutant" if mid in mut else "wild_type"


T["TP53_cmp"] = [call(m, tp_m, seq) for m in T.model_id]
T["BRAF_V600_cmp"] = [call(m, br_m, seq) for m in T.model_id]
T["TP53_cmp_protein"] = T.model_id.map(tp.groupby("model_id").protein_mutation.apply(lambda s: ";".join(sorted(set(map(str, s)))))).fillna("")
log(f"CMP TP53 calls: {T.TP53_cmp.value_counts().to_dict()}; BRAF V600: {T.BRAF_V600_cmp.value_counts().to_dict()}")

# ---------------- CCLE 2018 MAF ----------------
maf = pd.read_csv(os.path.join(DAT, "ccle/CCLE_DepMap_18q3_maf_20180718.txt"), sep="\t", low_memory=False,
                  usecols=["Hugo_Symbol", "Variant_Classification", "Protein_Change", "Tumor_Sample_Barcode", "Broad_ID"])
NONSILENT = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
             "Splice_Site", "Nonstop_Mutation", "Start_Codon_SNP", "Start_Codon_Del", "Start_Codon_Ins",
             "De_novo_Start_OutOfFrame", "De_novo_Start_InFrame"}
classes = maf.Variant_Classification.value_counts()
classes.rename("n_rows").to_csv(os.path.join(RES, "D_ccle_maf_variant_classes.tsv"), sep="\t")
log(f"MAF classes found: {sorted(classes.index)}; non-silent set used: {sorted(NONSILENT & set(classes.index))}")
maf_ach = set(maf.Broad_ID.dropna()); maf_bar = set(maf.Tumor_Sample_Barcode.dropna())
tpm = maf[(maf.Hugo_Symbol == "TP53") & maf.Variant_Classification.isin(NONSILENT)]
brm = maf[(maf.Hugo_Symbol == "BRAF") & maf.Protein_Change.astype(str).str.startswith("p.V600")]


def call_maf(ach, ccle, sub):
    if isinstance(ach, str) and ach in maf_ach:
        return "mutant" if ach in set(sub.Broad_ID) else "wild_type"
    if isinstance(ccle, str) and ccle in maf_bar:
        return "mutant" if ccle in set(sub.Tumor_Sample_Barcode) else "wild_type"
    return "unknown"


ccle_name = T.CCLE_ID.where(T.CCLE_ID.notna(), T.cmap_ccle_name.replace("", np.nan))
T["ccle_name"] = ccle_name
T["TP53_ccle"] = [call_maf(a, c, tpm) for a, c in zip(T.BROAD_ID, ccle_name)]
T["BRAF_V600_ccle"] = [call_maf(a, c, brm) for a, c in zip(T.BROAD_ID, ccle_name)]
T["TP53_cov"] = np.where(T.TP53_cmp != "unknown", T.TP53_cmp, T.TP53_ccle)
T["TP53_cov_source"] = np.where(T.TP53_cmp != "unknown", "CMP", np.where(T.TP53_ccle != "unknown", "CCLE_MAF", ""))
log(f"CCLE TP53 calls: {T.TP53_ccle.value_counts().to_dict()}; TP53_cov: {T.TP53_cov.value_counts().to_dict()}")
for c in ["TP53_cmp", "TP53_ccle", "TP53_cov", "BRAF_V600_cmp", "BRAF_V600_ccle"]:
    assert set(T[c]) <= {"mutant", "wild_type", "unknown"}

# ---------------- Task C ----------------
C, commitC = read_states("results/C/C_cellline_sex_labels.tsv")
Cc = C[["model_id", "BROAD_ID", "CCLE_ID", "RRID", "expressed_state"]].copy()
tc = []
for _, r in T.iterrows():
    hit = None; key = ""
    for k, v in [("BROAD_ID", r.BROAD_ID), ("CCLE_ID", r.ccle_name), ("RRID", r.RRID if isinstance(r.RRID, str) else r.cmap_cellosaurus_id)]:
        if isinstance(v, str) and v:
            m = Cc[Cc[k] == v]
            if len(m) == 1:
                hit = m.iloc[0].expressed_state; key = k; break
            assert len(m) <= 1, (k, v)
    tc.append((hit, key))
T["taskC_state"] = [h if h is not None else "no_row" for h, _ in tc]
T["taskC_match_key"] = [k for _, k in tc]
CALLED = {"XY", "XX-like", "Y-lost", "XIST-lost", "discordant"}
T["source_conflict"] = [bool(a in CALLED and c in CALLED and a != c) for a, c in zip(T.taskA_state, T.taskC_state)]
log(f"Task C matched {(T.taskC_state!='no_row').sum()}; source conflicts {T.source_conflict.sum()}")

# ---------------- Cellosaurus doubling time ----------------
dt = {}
acc = None
pat = re.compile(r"(\d+(?:\.\d+)?)(?:\s*(?:-|to)\s*(\d+(?:\.\d+)?))?\s*(hours?|days?|h)\b", re.I)
with open(os.path.join(DAT, "cello/cellosaurus.txt"), encoding="utf-8", errors="replace") as f:
    for line in f:
        if line.startswith("AC   "):
            acc = line[5:].strip()
        elif line.startswith("CC   Doubling time:") and acc:
            txt = line[len("CC   Doubling time:"):]
            for part in txt.split(";"):
                part = re.sub(r"\(.*?\)", "", part)
                m = pat.search(part)
                if m:
                    v = float(m.group(1)) if not m.group(2) else (float(m.group(1)) + float(m.group(2))) / 2
                    if m.group(3).lower().startswith("d"):
                        v *= 24
                    dt.setdefault(acc, []).append(v)
rr = T.RRID.where(T.RRID.notna(), T.cmap_cellosaurus_id.replace("", np.nan))
T["cellosaurus_ac"] = rr
T["doubling_time_h"] = [float(np.median(dt[a])) if isinstance(a, str) and a in dt else np.nan for a in rr]
T["n_doubling_values"] = [len(dt[a]) if isinstance(a, str) and a in dt else 0 for a in rr]
log(f"doubling time present {T.doubling_time_h.notna().sum()} of {len(T)}")

# ---------------- NR3C1 (CCLE RPKM) ----------------
with gzip.open(os.path.join(DAT, "ccle/CCLE_RNAseq_genes_rpkm_20180929.gct.gz"), "rt") as f:
    v = f.readline().strip(); dims = f.readline().strip(); hdr = f.readline().rstrip("\n").split("\t")
    nr = None
    for line in f:
        p = line.split("\t", 2)
        if p[1] == "NR3C1":
            assert nr is None, "two NR3C1 rows"
            nr = line.rstrip("\n").split("\t")
assert hdr[0] == "Name" and hdr[1] == "Description" and nr is not None
nrv = pd.Series(np.log2(np.array(nr[2:], dtype=float) + 1), index=hdr[2:])
T["NR3C1_log2rpkm1"] = T.ccle_name.map(nrv)
log(f"CCLE RPKM header '{v}' dims '{dims}'; NR3C1 gene id {nr[0]}; NR3C1 present for {T.NR3C1_log2rpkm1.notna().sum()} lines")

# ---------------- qc_pass fraction ----------------
for pn in PANELS:
    s = panel_sigs[pn]
    q = s.groupby("cell_iname").qc_pass.apply(lambda x: float((x == 1).mean()))
    T[f"qc_pass_frac_{pn}"] = T.cell_iname.map(q)
    T[f"n_sig_{pn}"] = T.cell_iname.map(s.groupby("cell_iname").size()).fillna(0).astype(int)

# ---------------- contrasts ----------------
def tf(x):
    return x if isinstance(x, bool) else (str(x).strip().lower() == "true")


T["elig_factorial"] = T.taskA_state.isin(CALLED)
T["fact_Y"] = np.where(T.elig_factorial, T.taskA_Yhi.map(tf), np.nan)
T["fact_X"] = np.where(T.elig_factorial, T.taskA_Xhi.map(tf), np.nan)
T["elig_strict"] = T.taskA_state.isin(["XY", "XX-like"])
T["strict_XY"] = np.where(T.elig_strict, (T.taskA_state == "XY").astype(float), np.nan)
T["elig_annotated"] = T.cmap_donor_sex.isin(["M", "F"])
T["annot_M"] = np.where(T.elig_annotated, (T.cmap_donor_sex == "M").astype(float), np.nan)

# ---------------- outputs ----------------
STATE_COLS = [c for c in T.columns if any(k in c for k in ["state", "sex", "Yhi", "Xhi", "XIST", "y_score", "conflict",
                                                             "confidence", "fact_", "strict_XY", "annot_M"])]
order = ["cell_iname", "lineage", "in_P1", "in_P1b", "in_P2",
         "taskA_state", "taskA_confidence", "taskA_y_score", "taskA_y_score_source", "taskA_Yhi", "taskA_Xhi",
         "taskA_XIST_ccle", "taskA_XIST_sanger", "taskA_XIST_cmpbroad", "taskC_state", "taskC_match_key", "source_conflict",
         "cmap_donor_sex", "taskA_cmap_donor_sex", "model_id", "cmp_model_name", "cmp_match_key", "RRID", "cellosaurus_ac", "BROAD_ID",
         "ccle_name", "TP53_cmp", "TP53_cmp_protein", "TP53_ccle", "TP53_cov", "TP53_cov_source", "BRAF_V600_cmp", "BRAF_V600_ccle",
         "ploidy", "ploidy_source", "ploidy_wes", "ploidy_snp6", "ploidy_wgs", "msi_status", "mutational_burden", "cmp_cancer_type",
         "primary_disease", "doubling_time_h", "n_doubling_values", "NR3C1_log2rpkm1",
         "qc_pass_frac_P1", "qc_pass_frac_P1b", "qc_pass_frac_P2", "n_sig_P1", "n_sig_P1b", "n_sig_P2",
         "elig_factorial", "fact_Y", "fact_X", "elig_strict", "strict_XY", "elig_annotated", "annot_M"]
assert set(order) <= set(T.columns), set(order) - set(T.columns)
assert T.cell_iname.is_unique
T[order].to_csv(os.path.join(RES, "D_line_table.tsv"), sep="\t", index=False)

cov_cols = ["cell_iname", "lineage", "in_P1", "in_P1b", "in_P2", "model_id", "RRID", "BROAD_ID", "ccle_name",
            "TP53_cmp", "TP53_ccle", "TP53_cov", "TP53_cov_source", "BRAF_V600_cmp", "BRAF_V600_ccle", "ploidy", "ploidy_source",
            "msi_status", "mutational_burden", "cmp_cancer_type", "doubling_time_h", "NR3C1_log2rpkm1",
            "qc_pass_frac_P1", "qc_pass_frac_P1b", "qc_pass_frac_P2"]
assert not (set(cov_cols) & set(STATE_COLS)), set(cov_cols) & set(STATE_COLS)
T[cov_cols].to_csv(os.path.join(RES, "D_line_covariates.tsv"), sep="\t", index=False)

elig, counts = [], []
CON = {"factorial": ("elig_factorial", ["fact_Y", "fact_X"]), "strict": ("elig_strict", ["strict_XY"]),
       "annotated": ("elig_annotated", ["annot_M"])}
for pn in PANELS:
    P = T[T[f"in_{pn}"]]
    for cn, (ec, gcols) in CON.items():
        E = P[P[ec]]
        for l in E.cell_iname:
            elig.append([pn, cn, l])
        for lin, g in E.groupby("lineage"):
            if cn == "factorial":
                y, x = g.fact_Y.astype(bool), g.fact_X.astype(bool)
                counts.append([pn, cn, lin, len(g), int((y & ~x).sum()), int((~y & x).sum()), int((~y & ~x).sum()), int((y & x).sum())])
            else:
                k = int(g[gcols[0]].sum())
                counts.append([pn, cn, lin, len(g), k, len(g) - k, np.nan, np.nan])
pd.DataFrame(elig, columns=["panel", "contrast", "cell_iname"]).to_csv(os.path.join(RES, "D_eligible_lines.tsv"), sep="\t", index=False)
cdf = pd.DataFrame(counts, columns=["panel", "contrast", "lineage", "n_eligible", "n_g1", "n_g2", "n_g3", "n_g4"])
cdf["group_definitions"] = cdf.contrast.map({"factorial": "g1=Y+X-;g2=Y-X+;g3=Y-X-;g4=Y+X+", "strict": "g1=XY;g2=XX-like",
                                             "annotated": "g1=annotated M;g2=annotated F"})
cdf.to_csv(os.path.join(RES, "D_lineage_counts.tsv"), sep="\t", index=False)

# ---------------- acceptance checks ----------------
chk = []
P1 = T[T.in_P1]; P2 = T[T.in_P2]
for st in ["XY", "XX-like"]:
    g = P1[P1.taskA_state == st]
    vc = g.TP53_cmp.value_counts()
    chk.append(["P1 (ASG) TP53 CMP among " + st, f"mutant {vc.get('mutant',0)}, wild type {vc.get('wild_type',0)}, unknown {vc.get('unknown',0)}",
                {"XY": "11, 2, 1", "XX-like": "5, 5, 3"}[st]])
for pn, P, pil in [("P1 (ASG)", P1, "45 of 52"), ("P2 (CPC)", P2, "42 of 47")]:
    called = P[P.elig_factorial]
    called_amb = P[P.taskA_state.isin(CALLED | {"ambiguous"})]
    chk.append([f"{pn} ploidy present among called lines (state in XY, XX-like, Y-lost, XIST-lost, discordant)",
                f"{called.ploidy.notna().sum()} of {len(called)}", pil])
    chk.append([f"{pn} ploidy present among lines with any Task A call incl. ambiguous",
                f"{called_amb.ploidy.notna().sum()} of {len(called_amb)}", pil])
pd.DataFrame(chk, columns=["check", "this_session", "pilot"]).to_csv(os.path.join(RES, "D_acceptance_checks.tsv"), sep="\t", index=False)
pd.DataFrame([["task-A", "results/A/A_cellline_sex_states.tsv", commitA], ["task-C", "results/C/C_cellline_sex_labels.tsv", commitC]],
             columns=["analysis", "path", "sha256"]).to_csv(os.path.join(RES, "D_line_table_sources.tsv"), sep="\t", index=False)
for r in chk:
    log(" | ".join(map(str, r)))
log("done")
