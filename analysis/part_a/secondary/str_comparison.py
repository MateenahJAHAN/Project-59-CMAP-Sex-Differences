"""Compare expressed Y scores in five sources to Cellosaurus amelogenin STR records."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

work = Path(sys.argv[1]).resolve()
res = work / "results"
out = res / "str"


def read(path):
    return pd.read_csv(path, sep="\t")


def classify(score, cutoff):
    return pd.Series(np.where(score > cutoff, "Y expressed", "Y below cut-off"), index=score.index).where(score.notna())


records = read(out / "amelogenin_by_record.tsv").set_index("accession")
assert records.index.is_unique and len(records) == 1102
table = records[["amelogenin_class"]].copy()
f1 = res / "F1"
local = res / "tahoe/latest"
datasets = {"Tahoe-100M": (local / "partA_local_states.tsv", "local"),
            "EmeraldBay": (f1 / "F1_emerald_states.tsv", "EmeraldBay"),
            "sci-Plex": (f1 / "F1_sciplex_states.tsv", "CCLE")}
for dataset, (path, source) in datasets.items():
    states = read(path)
    if source == "local":
        cut = read(local / "partA_local_cutoffs.tsv").set_index("score").loc["y_score", "cutoff"]
        table["in_tahoe"] = table.index.isin(states.cellosaurus_id)
    else:
        cuts = read(f1 / "F1_cutoffs.tsv")
        rows = cuts[(cuts.dataset.eq("EmeraldBay") if dataset == "EmeraldBay" else cuts.dataset.ne("EmeraldBay"))
                    & cuts.score.eq("y_score")]
        cut, = rows.cutoff.unique()
    states["call"] = classify(states.y_score, cut)
    calls = states.groupby("cellosaurus_id").call.agg(lambda values: values.iloc[0] if values.nunique() == 1 else "mixed")
    table[dataset + "_Y"] = calls

a = read(res / "A/A_cellline_sex_states.tsv").dropna(subset=["cellosaurus_ac"]).set_index("cellosaurus_ac")
c = read(res / "C/C_cellline_sex_labels.tsv").dropna(subset=["RRID"]).set_index("RRID")
acuts = read(res / "A/A_cutoffs.tsv").set_index(["source", "axis"]).cutoff
ccuts = read(res / "C/C_expression_cutoffs.tsv").set_index(["expr_source", "variable"]).cutoff_log2
for dataset, fallback in (("CCLE", ["ccle"]), ("CMP", ["sanger", "cmpbroad"])):
    name = dataset + "_Y"
    table[name] = classify(c[dataset.lower() + "_yscore"], ccuts[dataset, "yscore"])
    for source in fallback:
        calls = classify(a["y_score_" + source], acuts[source, "y_score"])
        if calls.groupby(level=0).nunique().gt(1).any():
            raise ValueError("Multiple calls for one accession")
        table[name] = table[name].fillna(calls.groupby(level=0).first())

summary = []
columns = ["amelogenin_class", "in_tahoe"]
for dataset in (*datasets, "CCLE", "CMP"):
    call = table[dataset + "_Y"]
    allele = table.amelogenin_class
    eligible = allele.isin(["X only", "X,Y"]) & call.isin(["Y expressed", "Y below cut-off"])
    decision = pd.Series("not compared", index=table.index)
    decision[eligible] = "agree"
    decision[eligible & allele.eq("X,Y") & call.eq("Y below cut-off")] = "type A"
    decision[eligible & allele.eq("X only") & call.eq("Y expressed")] = "type B"
    table[dataset + "_vs_STR"] = decision
    columns.extend([dataset + "_Y", dataset + "_vs_STR"])
    summary.append(dict(source=dataset, n_compared=int(eligible.sum()), n_agree=int(decision.eq("agree").sum()),
                        n_type_A=int(decision.eq("type A").sum()), n_type_B=int(decision.eq("type B").sum())))
table[columns].to_csv(out / "figure3_supp_table_amelogenin_by_line.tsv", sep="\t")
print(pd.DataFrame(summary).to_string(index=False), flush=True)
pd.DataFrame(summary).to_csv(out / "figure3_agreement_by_source.tsv", sep="\t", index=False)
