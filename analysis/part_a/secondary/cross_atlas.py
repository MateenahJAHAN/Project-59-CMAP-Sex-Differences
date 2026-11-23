"""A2: EmeraldBay vehicle-cell states (task F1) against the frozen full-data Tahoe Part A table.

Frozen rule (analysis_plan/PartA_frozen_rules.md, section B, row A2, commit 7cac9a4):
  "EmeraldBay calls agree with Tahoe's for at least 90% of shared lines of high or medium confidence
   (task F1 compares with the pilot table; repeat against the full-data table)".
Implementation mirrors F1 (src/F1/analyze.py lines 525-530): both confidences in {high, medium}; agreement on state.
Tahoe side: state_final and confidence_final of the frozen consensus table (commit 1ba488a).
"""
import hashlib, json
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(sys.argv[1]).resolve()
EM = ROOT / "results/F1/F1_emerald_states.tsv"
TH = ROOT / "results/consensus/tahoe_partA_states_consensus.tsv"
OUT = ROOT / "results/consensus"
def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()

em = pd.read_csv(EM, sep="\t")
th = pd.read_csv(TH, sep="\t")
assert em.cellosaurus_id.is_unique and th.cellosaurus_id.is_unique
assert len(th) == 50 and th.callers_agree.all()
STATES = {"XY", "XX-like", "Y-lost", "XIST-lost", "discordant", "ambiguous"}
CONF = {"high", "medium", "low"}
assert set(em.state) <= STATES and set(th.state_final) <= STATES
assert set(em.confidence) <= CONF and set(th.confidence_final) <= CONF

m = em.merge(th[["cellosaurus_id", "cell_name", "state_final", "confidence_final"]], on="cellosaurus_id",
             how="inner", suffixes=("_em", "_th"), validate="one_to_one")
# names may be spelled differently between atlases; report, do not join on them
name_diff = m.loc[m.cell_name_em != m.cell_name_th, ["cellosaurus_id", "cell_name_em", "cell_name_th"]]
m["agree"] = m.state == m.state_final
hm = m[m.confidence.isin(["high", "medium"]) & m.confidence_final.isin(["high", "medium"])]

def row(name, sub, rule, thr):
    k, n = int(sub.agree.sum()), len(sub)
    dis = ";".join(f"{r.cell_name_em}: EmeraldBay {r.state} vs Tahoe {r.state_final}" for r in sub[~sub.agree].itertuples())
    res = ("pass" if k / n >= thr else "fail") if thr is not None else "descriptive"
    return dict(comparison="EmeraldBay vs Tahoe full-data frozen table", subset=name, n_compared=n, n_agree=k,
                pct_agree=k / n, rule=rule, result=res, disagreements=dis)

rows = [row("lines with high or medium confidence in both", hm, ">= 90% agreement (pre-registered, A2)", 0.90),
        row("all shared lines", m, "descriptive", None)]
R = pd.DataFrame(rows)
R.to_csv(OUT / "A2_emerald_vs_full_agreement.tsv", sep="\t", index=False)
m.sort_values("cell_name_em").to_csv(OUT / "A2_emerald_vs_full_lines.tsv", sep="\t", index=False)

# what changed relative to F1's pilot comparison (descriptive)
pil = pd.read_csv(ROOT / "results/F1/F1_agreement.tsv", sep="\t")
pil = pil[pil.comparison == "EmeraldBay vs Tahoe pilot state"]
summary = dict(inputs_sha256={"F1_emerald_states.tsv": sha(EM), "tahoe_partA_states_consensus.tsv": sha(TH)},
               n_emerald_lines=len(em), n_tahoe_lines=len(th), n_shared=len(m), n_hm_both=len(hm),
               a2_k=int(hm.agree.sum()), a2_n=len(hm), a2_pct=float(hm.agree.mean()), a2_pass=bool(hm.agree.mean() >= 0.90),
               all_shared_k=int(m.agree.sum()), all_shared_n=len(m),
               pilot_comparison=pil[["subset", "n_compared", "n_agree", "result"]].to_dict("records"),
               name_spelling_differences=name_diff.to_dict("records"))
(OUT / "A2_summary.json").write_text(json.dumps(summary, indent=1))
print(R.to_string(index=False))
print(json.dumps({k: v for k, v in summary.items() if k != "inputs_sha256"}, indent=1))
print("excluded from A2 (low confidence in either):")
print(m.loc[~m.index.isin(hm.index), ["cell_name_em", "state", "confidence", "state_final", "confidence_final"]].to_string(index=False))
