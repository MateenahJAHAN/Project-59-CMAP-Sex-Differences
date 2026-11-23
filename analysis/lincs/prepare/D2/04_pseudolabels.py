# Copied for task D2 from origin/task-D commit a84f566 (a84f566f5a3aa97e443f23488c1ac0c3395f4110), src/D/04_pseudolabels.py.
# D2 changes: the line sets, counts, numbers of pseudo-labellings and seeds of runs R1 and R2
# (analysis_plan/D2_plan.md, sections 2 and 3). The drawing code (`draw`) is task D's loop, unchanged.
"""Task D2 script 4: pseudo-labels from per-lineage group counts only.

R1 (seed 20261120, 200 pseudo-labellings, factorial codes g1..g4):
  F41 and F44: lines of results/D2/D2_sets_lines.tsv, counts of results/D2/D2_sets_counts.tsv.
  A new generator default_rng(20261120) is made for each set.
R2 (seed 20261121, 100 pseudo-labellings, annotated codes g1 = M, g2 = F):
  P1b and P2: task D's eligible lines (D_eligible_lines.tsv) and per-lineage counts (D_lineage_counts.tsv)
  for the annotated contrast. A new generator default_rng(20261121) is made for each panel.
Within each lineage, the lines (sorted by lineage and name) receive a uniformly random arrangement of the
lineage's group codes with exactly its counts. The true assignment is never read.

Output: results/D2/pseudolabels/D2_pseudo_<tag>.tsv.gz (rows: lines; columns s000..: group code)
"""
import os, sys, time
import numpy as np, pandas as pd

ROOT = os.environ["P59_LINCS_WORK"]
RES = os.path.join(ROOT, "results", "D2")
IN_D = os.path.join(ROOT, "data", "D2", "inputs", "taskD")
OUT = os.path.join(RES, "pseudolabels")
os.makedirs(OUT, exist_ok=True)
t0 = time.time()
cl = pd.read_csv(os.path.join(IN_D, "D_line_covariates.tsv"), sep="\t", usecols=["cell_iname", "lineage"])
assert cl.cell_iname.is_unique


def draw(e, c, S, rng, tag):
    """Task D's loop: e = lines with lineage; c = counts indexed by lineage (n_g1..n_g4, n_eligible)."""
    assert set(c.index) == set(e.lineage), tag
    e = e.sort_values(["lineage", "cell_iname"]).reset_index(drop=True)
    lab = np.zeros((len(e), S), dtype=np.int8)
    for lin, idx in e.groupby("lineage").groups.items():
        r = c.loc[lin]
        codes = np.concatenate([np.full(int(r[f"n_g{g}"]), g) for g in (1, 2, 3, 4) if pd.notna(r[f"n_g{g}"])])
        assert len(codes) == len(idx) == r.n_eligible, (tag, lin)
        for s in range(S):
            lab[np.asarray(idx), s] = rng.permutation(codes)
    df = pd.DataFrame(lab, columns=[f"s{s:03d}" for s in range(S)])
    df.insert(0, "cell_iname", e.cell_iname.values)
    df.insert(1, "lineage", e.lineage.values)
    # every pseudo-labelling keeps the per-lineage counts
    for s in (0, S - 1):
        got = df.groupby("lineage")[f"s{s:03d}"].apply(lambda v: tuple(np.bincount(v, minlength=5)[1:]))
        for lin, tup in got.items():
            want = tuple(int(c.loc[lin, f"n_g{g}"]) if pd.notna(c.loc[lin, f"n_g{g}"]) else 0 for g in (1, 2, 3, 4))
            assert tup == want, (tag, lin, tup, want)
    df.to_csv(os.path.join(OUT, f"D2_pseudo_{tag}.tsv.gz"), sep="\t", index=False)
    return df


written = []
# ---- R1: fixed sets F41 and F44 ----
sl = pd.read_csv(os.path.join(RES, "D2_sets_lines.tsv"), sep="\t")
sc = pd.read_csv(os.path.join(RES, "D2_sets_counts.tsv"), sep="\t")
assert set(sl.columns) == {"set", "cell_iname", "lineage", "TP53_cov"}
for name in ("F41", "F44"):
    e = sl[sl.set == name][["cell_iname", "lineage"]]
    e2 = e.merge(cl, on="cell_iname", how="left", suffixes=("", "_cov"), validate="1:1")
    assert len(e2) == len(e) and (e2.lineage == e2.lineage_cov).all()
    c = sc[sc.set == name].set_index("lineage").rename(columns={"n_lines": "n_eligible"})
    df = draw(e, c, 200, np.random.default_rng(20261120), f"R1_{name}")
    written.append((f"R1_{name}", len(df), 200))
# ---- R2: task D's eligible lines and counts, annotated contrast, P1b and P2 ----
cnt = pd.read_csv(os.path.join(IN_D, "D_lineage_counts.tsv"), sep="\t")
el = pd.read_csv(os.path.join(IN_D, "D_eligible_lines.tsv"), sep="\t")
assert set(el.columns) == {"panel", "contrast", "cell_iname"}
for pn in ("P1b", "P2"):
    e = el[(el.panel == pn) & (el.contrast == "annotated")].merge(cl, on="cell_iname", how="left", validate="1:1")
    assert e.lineage.notna().all()
    c = cnt[(cnt.panel == pn) & (cnt.contrast == "annotated")].set_index("lineage")
    df = draw(e[["cell_iname", "lineage"]], c, 100, np.random.default_rng(20261121), f"R2_{pn}_annotated")
    written.append((f"R2_{pn}_annotated", len(df), 100))
with open(os.path.join(ROOT, "logs", "D2_runlog.txt"), "a") as f:
    f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} [04_pseudolabels +{time.time()-t0:.0f}s] "
            + "; ".join(f"{t}: {n} lines x {S} pseudo-labellings" for t, n, S in written) + "\n")
print(written)
