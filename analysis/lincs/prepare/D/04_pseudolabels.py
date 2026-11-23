"""Task D script 4: pseudo-labels from per-lineage group counts only.

Reads results/D/D_lineage_counts.tsv and results/D/D_eligible_lines.tsv. Within
each lineage, the eligible lines (sorted by name) receive a uniformly random
arrangement of the lineage's group labels with exactly its group counts. The
true assignment is never read, so it is never evaluated or saved.

Output: results/D/pseudolabels/D_pseudo_<panel>_<contrast>.tsv.gz
(rows: eligible lines; columns s000..s199: group code 1..4 as in D_lineage_counts.tsv)
"""
import os, time
import numpy as np, pandas as pd

ROOT = os.environ["P59_LINCS_WORK"]
RES = os.path.join(ROOT, "results", "D")
OUT = os.path.join(RES, "pseudolabels")
os.makedirs(OUT, exist_ok=True)
S = 200
CONTRASTS = {"factorial": 0, "strict": 1, "annotated": 2}
PANEL_OFF = {"P1": 0, "P1b": 10, "P2": 20}
t0 = time.time()
cnt = pd.read_csv(os.path.join(RES, "D_lineage_counts.tsv"), sep="\t")
el = pd.read_csv(os.path.join(RES, "D_eligible_lines.tsv"), sep="\t")
assert set(el.columns) == {"panel", "contrast", "cell_iname"}
cl = pd.read_csv(os.path.join(RES, "D_line_covariates.tsv"), sep="\t", usecols=["cell_iname", "lineage"])
assert cl.cell_iname.is_unique
for pn in PANEL_OFF:
    for cn, ci in CONTRASTS.items():
        e = el[(el.panel == pn) & (el.contrast == cn)].merge(cl, on="cell_iname", how="left", validate="1:1")
        assert e.lineage.notna().all()
        c = cnt[(cnt.panel == pn) & (cnt.contrast == cn)].set_index("lineage")
        assert set(c.index) == set(e.lineage), (pn, cn)
        rng = np.random.default_rng(20261117 + ci + PANEL_OFF[pn])
        e = e.sort_values(["lineage", "cell_iname"]).reset_index(drop=True)
        lab = np.zeros((len(e), S), dtype=np.int8)
        for lin, idx in e.groupby("lineage").groups.items():
            r = c.loc[lin]
            codes = np.concatenate([np.full(int(r[f"n_g{g}"]), g) for g in (1, 2, 3, 4) if pd.notna(r[f"n_g{g}"])])
            assert len(codes) == len(idx) == r.n_eligible, (pn, cn, lin)
            for s in range(S):
                lab[np.asarray(idx), s] = rng.permutation(codes)
        df = pd.DataFrame(lab, columns=[f"s{s:03d}" for s in range(S)])
        df.insert(0, "cell_iname", e.cell_iname.values)
        df.insert(1, "lineage", e.lineage.values)
        df.to_csv(os.path.join(OUT, f"D_pseudo_{pn}_{cn}.tsv.gz"), sep="\t", index=False)
with open(os.path.join(ROOT, "logs", "D_runlog.txt"), "a") as f:
    f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} [04_pseudolabels +{time.time()-t0:.0f}s] {S} pseudo-labellings per panel x contrast written\n")
print("done", time.time() - t0)
