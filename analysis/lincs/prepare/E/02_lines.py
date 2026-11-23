"""Step 1b of task E: line annotations (lineage, TP53 status). Reads no expression value.

Blinding: reads only the columns listed in the usecols lists below. Asserts that no sex-state,
Y, X, donor-sex or gender column is among them.

Lineage: cellinfo_beta.txt column cell_lineage.
TP53: Cell Model Passports mutations_summary_latest.csv.gz.
  mutant  = at least one TP53 row with coding == 't' (all such effects in this file are non-silent:
            missense, nonsense, frameshift, ess_splice, inframe, start_lost, stop_lost)
  wild type = model present in the file (any gene) and no such TP53 row
  unknown = line not mapped to a CMP model, or model absent from the file
Mapping LINCS line -> CMP model_id (model_list_20240110.csv), in order:
  1) cellinfo cellosaurus_id == model_list RRID
  2) cellinfo ccle_name == model_list CCLE_ID
  3) normalised name (upper case, alphanumerics only) of cell_iname == model_name or a synonym
Output: results/E/E_line_annotations.tsv
"""
import sys, re
import numpy as np, pandas as pd

ROOT = sys.argv[1]
FORBIDDEN = {"donor_sex", "gender", "sex", "ethnicity", "donor_ethnicity"}
CELL_COLS = ["cell_iname", "cellosaurus_id", "cell_lineage", "primary_disease", "ccle_name"]
ML_COLS = ["model_id", "model_name", "synonyms", "RRID", "CCLE_ID", "BROAD_ID", "cancer_type"]
MUT_COLS = ["model_id", "gene_symbol", "coding", "effect", "protein_mutation"]
for cols in (CELL_COLS, ML_COLS, MUT_COLS):
    assert not FORBIDDEN & set(cols)

cell = pd.read_csv(f"{ROOT}/data/E/meta/cellinfo_beta.txt", sep="\t", usecols=CELL_COLS, dtype=str)
assert cell.cell_iname.is_unique
ml = pd.read_csv(f"{ROOT}/data/E/cmp/model_list_20240110.csv", usecols=ML_COLS, dtype=str)
assert ml.model_id.is_unique
mut = pd.read_csv(f"{ROOT}/data/E/cmp/mutations_summary_latest.csv.gz", usecols=MUT_COLS, dtype=str)
assert set(mut.coding.dropna()) <= {"t", "f"}

wells = pd.read_csv(f"{ROOT}/results/E/E_wells.tsv.gz", sep="\t", usecols=["panel", "cell_iname"], dtype=str)
lines = wells.drop_duplicates().groupby("cell_iname").panel.agg(lambda x: "|".join(sorted(set(x)))).reset_index()
lines = lines.rename(columns={"panel": "panels"})
n0 = len(lines)
lines = lines.merge(cell, on="cell_iname", how="left", validate="one_to_one")
assert len(lines) == n0 and lines.cell_lineage.notna().all()


def norm(s):
    return re.sub(r"[^A-Z0-9]", "", str(s).upper())


rrid = {r: m for r, m in zip(ml.RRID, ml.model_id) if isinstance(r, str)}
ccle = {r: m for r, m in zip(ml.CCLE_ID, ml.model_id) if isinstance(r, str)}
names = {}
for m, n, syn in zip(ml.model_id, ml.model_name, ml.synonyms):
    for x in [n] + (str(syn).split(";") if isinstance(syn, str) else []):
        names.setdefault(norm(x), set()).add(m)
mid, how = [], []
for _, r in lines.iterrows():
    if isinstance(r.cellosaurus_id, str) and r.cellosaurus_id in rrid:
        mid.append(rrid[r.cellosaurus_id]); how.append("RRID")
    elif isinstance(r.ccle_name, str) and r.ccle_name in ccle:
        mid.append(ccle[r.ccle_name]); how.append("CCLE_ID")
    elif len(names.get(norm(r.cell_iname), set())) == 1:
        mid.append(next(iter(names[norm(r.cell_iname)]))); how.append("name")
    else:
        mid.append(np.nan); how.append("unmapped")
lines["cmp_model_id"] = mid
lines["cmp_map_method"] = how
lines = lines.merge(ml[["model_id", "cancer_type"]].rename(columns={"model_id": "cmp_model_id", "cancer_type": "cmp_cancer_type"}),
                    on="cmp_model_id", how="left", validate="many_to_one")
assert len(lines) == n0
dup = lines.cmp_model_id.dropna().duplicated().sum()
print("lines mapped to the same CMP model as another line:", int(dup))

in_file = set(mut.model_id)
tp = mut[(mut.gene_symbol == "TP53") & (mut.coding == "t")]
tp_mut = tp.groupby("model_id").protein_mutation.agg(lambda x: ";".join(sorted(set(map(str, x)))))
status = []
for m in lines.cmp_model_id:
    if not isinstance(m, str) or m not in in_file:
        status.append("unknown")
    elif m in tp_mut.index:
        status.append("mutant")
    else:
        status.append("wild_type")
lines["tp53_status"] = status
lines["tp53_protein_changes"] = lines.cmp_model_id.map(tp_mut).fillna("")
assert set(lines.tp53_status) <= {"mutant", "wild_type", "unknown"}
lines = lines[["cell_iname", "panels", "cell_lineage", "primary_disease", "cellosaurus_id", "ccle_name",
               "cmp_model_id", "cmp_map_method", "cmp_cancer_type", "tp53_status", "tp53_protein_changes"]]
lines.to_csv(f"{ROOT}/results/E/E_line_annotations.tsv", sep="\t", index=False)
for p in ["P1", "P2"]:
    d = lines[lines.panels.str.contains(p)]
    print(p, len(d), d.cmp_map_method.value_counts().to_dict(), d.tp53_status.value_counts().to_dict(),
          "lineages", d.cell_lineage.nunique())
