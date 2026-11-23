"""Diagnostic for exploratory X1 (EXPLORATORY): why X1b (level 5 minus the line's DMSO response) differs from X1a (level 3).

For the P1 X term on the 41-line set (frozen design), compute each compound's 978-gene X coefficient vector (gcore.run_term,
no permutations) on: level 5 as is, X1a (level 3, plate DMSO), X1b (level 5 minus line DMSO). Also the DMSO pseudo-compound's
X coefficient vector on level 5. Report, per tensor, the correlation over genes between each compound's X vector and minus
the DMSO X vector, and the mean T_c. If X1b's excess comes from subtracting one shared vector, its compounds' X vectors
will track minus the DMSO vector, and X1a's will not.
Usage (repository root): <task G python> -I src/G_run/diag_X1b_common_vector.py
"""
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
ROOT = os.getcwd()
from study import open_study
import numpy as np
import pandas as pd
import gcore as C
import gsets as GS

DMSO = "DMSO (pseudo-compound)"
cfg, I = open_study(ROOT)
cov = I.covariates()
P1 = I.panel("P1", cov, cfg.run("B_full"))
prim = GS.primary_settings(cfg)
given = pd.read_csv(os.path.join(ROOT, "data/G/inputs/taskD/D_line_table.tsv"), sep="\t")
key_seed = cfg.run("key_seed_base") + cfg.const("panel_numbers")["P1"]
lp, cp = I.file("taskD/tensors/P1_lines.tsv"), I.file("taskD/tensors/P1_compounds.tsv")
tens = {"level5 (confirmatory)": I.d_tensor("P1_resp.npy"),
        "X1a level 3, plate DMSO": os.path.join(ROOT, "data/G/exploratory/X1a_level3_plateDMSO_P1_Dorder.npy"),
        "X1b level 5 minus line DMSO": os.path.join(ROOT, "data/G/exploratory/X1b_level5_minus_lineDMSO_P1_Dorder.npy")}
d = GS.Aligned(P1, given).design("factX", prim)

# DMSO pseudo-compound X vector on level 5 (as in the confirmatory run: extra target)
df, _, recs = C.run_term(P1, d["Xo"], d["k"], d["mask"], 2, extra={DMSO: P1.dmso}, lin=d["lin"], rule=d["rule"], want_perm=False)
bd = [r for r in recs if r["compound"] == DMSO][0]["_b"]
rows = []
for name, path in tens.items():
    Q = C.Panel("P1", path, lp, cp, P1.genes, cov, key_seed, cfg.run("B_full"), None)
    df, _, recs = C.run_term(Q, d["Xo"], d["k"], d["mask"], 2, lin=d["lin"], rule=d["rule"], want_perm=False)
    for r in recs:
        if not r["tested"]:
            continue
        b = r["_b"]
        rows.append(dict(tensor=name, compound=r["compound"], corr_with_minus_DMSO_X=float(np.corrcoef(b, -bd)[0, 1]),
                         corr_with_plus_DMSO_X=float(np.corrcoef(b, bd)[0, 1]), T_obs=float(np.mean(r["_t"] ** 2))))
out = pd.DataFrame(rows)
os.makedirs(os.path.join(ROOT, "results/G_exploratory"), exist_ok=True)
out.assign(label="EXPLORATORY X1 diagnostic").to_csv(os.path.join(ROOT, "results/G_exploratory/X1_common_vector_by_compound.tsv"), sep="\t", index=False)
summ = out.groupby("tensor").agg(n_compounds=("compound", "size"), median_corr_minus_DMSO=("corr_with_minus_DMSO_X", "median"),
                                 share_corr_minus_DMSO_gt_05=("corr_with_minus_DMSO_X", lambda s: float((s > 0.5).mean())),
                                 median_corr_plus_DMSO=("corr_with_plus_DMSO_X", "median"), mean_T_obs=("T_obs", "mean")).reset_index()
summ.assign(label="EXPLORATORY X1 diagnostic").to_csv(os.path.join(ROOT, "results/G_exploratory/X1_common_vector_summary.tsv"), sep="\t", index=False)
print(summ.to_string(index=False))
