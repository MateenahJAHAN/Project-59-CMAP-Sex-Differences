"""Task D script 3: regulon tables for the secondary TF readout (task G).
Reads CollecTRI (OmniPath), the Lopes-Ramos 2026 supplement (mmc5.xlsx) and the
landmark list. No state and no expression value is read.
Outputs: results/D/D_regulon_collectri_landmark.tsv, results/D/D_lopesramos_tfs.tsv,
results/D/D_regulon_summary.tsv. The PANDA-LCL table is used only if
inputs/D/regulon_pandaLCL_top30_landmark.tsv exists."""
import os, time
import pandas as pd
ROOT = os.environ["P59_LINCS_WORK"]
RES = os.path.join(ROOT, "results", "D")
LOG = os.path.join(ROOT, "logs", "D_runlog.txt")
t0 = time.time()
lm = pd.read_csv(os.path.join(RES, "D_landmarks.tsv"), sep="\t")
lmset = set(lm.gene_symbol)
ct = pd.read_csv(os.path.join(ROOT, "data/D/reg/collectri_omnipath.tsv"), sep="\t")
assert {"source_genesymbol", "target_genesymbol", "is_stimulation", "is_inhibition", "consensus_stimulation", "consensus_inhibition"} <= set(ct.columns)
n_all = len(ct)
ct = ct[ct.target_genesymbol.isin(lmset)].copy()
# sign: consensus direction; +1 stimulation, -1 inhibition, 0 if neither/both
ct["sign"] = ct.consensus_stimulation.astype(int) - ct.consensus_inhibition.astype(int)
ct = ct[["source_genesymbol", "target_genesymbol", "sign", "is_stimulation", "is_inhibition"]].drop_duplicates(["source_genesymbol", "target_genesymbol"])
ct.rename(columns={"source_genesymbol": "tf", "target_genesymbol": "target"}).to_csv(os.path.join(RES, "D_regulon_collectri_landmark.tsv"), sep="\t", index=False)
lr = pd.read_excel(os.path.join(ROOT, "data/D/reg/mmc5.xlsx"), sheet_name="Sheet1")
assert list(lr.columns[:3]) == ["TF", "Number of tissues", "Chromosome"], lr.columns
lr.to_csv(os.path.join(RES, "D_lopesramos_tfs.tsv"), sep="\t", index=False)
rows = [["collectri_all_interactions", n_all], ["collectri_landmark_target_edges", len(ct)],
        ["collectri_tfs_with_landmark_targets", ct.source_genesymbol.nunique() if "source_genesymbol" in ct else ct.tf.nunique()],
        ["lopesramos_tfs", len(lr)], ["lopesramos_tfs_with_collectri_landmark_targets", int(lr.TF.isin(set(ct.source_genesymbol)).sum())]]
pp = os.path.join(os.environ["P59_PUBLIC_REPOSITORY"], "analysis/lincs/inputs/D/regulon_pandaLCL_top30_landmark.tsv")
rows.append(["pandaLCL_input_present", int(os.path.exists(pp))])
pd.DataFrame(rows, columns=["item", "value"]).to_csv(os.path.join(RES, "D_regulon_summary.tsv"), sep="\t", index=False)
with open(LOG, "a") as f:
    f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} [03_regulons +{time.time()-t0:.0f}s] {rows}\n")
print(rows)
