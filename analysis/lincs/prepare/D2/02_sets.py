"""Task D2 step 2: the line sets of runs R1, R3, R4 and R5 and their per-lineage group counts.

New for D2. This is the only D2 script that opens task D's line table, which holds state columns.
It reads no response value (no tensor, no response-size file), and it writes no state or label
next to a line name: the line-list output has no group column, and the counts output has only
per-lineage counts (as task D's D_lineage_counts.tsv does).

Sets (task D's analysis-set rule, src/D/common.py analysis_set, applied to the real eligibility):
  F44  P1 factorial: factorial-eligible P1 lines with a TP53_cov call, in lineages with >= 2 such
       lines in which Y or X varies (task D's 44-line set).
  F41  F44 without BEN, HCC1588 and MCF7 (tasks/D2_blind_followups.md, R1), with the same rule
       re-applied (a lineage left with < 2 lines or no variation would drop out).
  A25  P1 annotated: annotated-eligible P1 lines with a TP53_cov call, in lineages with >= 2 such
       lines in which annotated M vs F varies (task D's P1 annotated analysis set).
Checks: each set must equal the set derived without any state column, from D_eligible_lines.tsv,
D_lineage_counts.tsv and TP53_cov (possible because every line of the varying lineages has a
TP53_cov call; asserted), and F44 / A25 counts must equal D_lineage_counts.tsv for those lineages.

Outputs: results/D2/D2_sets_lines.tsv (set, cell_iname, lineage, TP53_cov)
         results/D2/D2_sets_counts.tsv (set, lineage, n_lines, n_g1..n_g4, group_definitions)
         results/D2/D2_sets_summary.tsv
Usage: python 02_sets.py <repo_root>
"""
import sys, time
import numpy as np, pandas as pd

ROOT = sys.argv[1]
IN = f"{ROOT}/data/D2/inputs/taskD"
OUT = f"{ROOT}/results/D2"
DROP = ["BEN", "HCC1588", "MCF7"]
t0 = time.time()

# Columns read from the line table. The state-derived columns (elig_*, fact_Y, fact_X, annot_M) are
# used only to apply the analysis-set rule and to count groups per lineage.
USE = ["cell_iname", "lineage", "in_P1", "TP53_cov", "elig_factorial", "fact_Y", "fact_X", "elig_annotated", "annot_M",
       "BROAD_ID", "ccle_name"]  # identifiers: for the R4 check that the state-free covariate table carries the same IDs
lt = pd.read_csv(f"{IN}/D_line_table.tsv", sep="\t", usecols=USE)
assert lt.cell_iname.is_unique and len(lt) == 129, len(lt)
cov = pd.read_csv(f"{IN}/D_line_covariates.tsv", sep="\t", usecols=["cell_iname", "lineage", "TP53_cov", "BROAD_ID", "ccle_name"])
m = lt.merge(cov, on="cell_iname", how="inner", suffixes=("", "_cov"), validate="1:1")
assert len(m) == len(lt)
assert (m.lineage == m.lineage_cov).all()
assert ((m.TP53_cov == m.TP53_cov_cov) | (m.TP53_cov.isna() & m.TP53_cov_cov.isna())).all()
ids_equal = all(((m[c] == m[f"{c}_cov"]) | (m[c].isna() & m[f"{c}_cov"].isna())).all() for c in ("BROAD_ID", "ccle_name"))
cov = cov[["cell_iname", "lineage", "TP53_cov"]]
assert set(m.TP53_cov.dropna()) <= {"mutant", "wild_type", "unknown"}, set(m.TP53_cov.dropna())
KNOWN = ["mutant", "wild_type"]  # as task D: "unknown" lines are excluded from the models
p1 = m[m.in_P1.astype(bool)].copy()
assert len(p1) == 85

GROUP_DEF = {"factorial": "g1=Y+X-;g2=Y-X+;g3=Y-X-;g4=Y+X+", "annotated": "g1=annotated M;g2=annotated F"}


def codes(df, contrast):
    if contrast == "factorial":
        y, x = df.fact_Y.astype(float).values, df.fact_X.astype(float).values
        assert set(np.unique(y)) <= {0.0, 1.0} and set(np.unique(x)) <= {0.0, 1.0}
        return np.where((y == 1) & (x == 0), 1, np.where((y == 0) & (x == 1), 2, np.where((y == 0) & (x == 0), 3, 4)))
    a = df.annot_M.astype(float).values
    assert set(np.unique(a)) <= {0.0, 1.0}
    return np.where(a == 1, 1, 2)


def analysis_set(df, contrast):
    """task D's rule: model columns known; lineages with >= 2 lines in which a tested term varies."""
    ok = df[df.TP53_cov.isin(KNOWN)].copy()
    ok["_g"] = codes(ok, contrast)
    keep = []
    for lin, g in ok.groupby("lineage"):
        if contrast == "factorial":
            y = np.isin(g._g, [1, 4]); x = np.isin(g._g, [2, 4])
            varies = (np.unique(y).size > 1) or (np.unique(x).size > 1)
        else:
            varies = np.unique(g._g).size > 1
        if len(g) >= 2 and varies:
            keep.append(lin)
    return ok[ok.lineage.isin(keep)]


def counts(s, contrast, name):
    rows = []
    for lin, g in s.groupby("lineage"):
        c = np.bincount(g._g, minlength=5)
        r = dict(set=name, lineage=lin, n_lines=len(g), n_g1=int(c[1]), n_g2=int(c[2]),
                 n_g3=int(c[3]) if contrast == "factorial" else np.nan, n_g4=int(c[4]) if contrast == "factorial" else np.nan,
                 group_definitions=GROUP_DEF[contrast])
        rows.append(r)
    return pd.DataFrame(rows)


fe = p1[p1.elig_factorial.astype(bool)]
ae = p1[p1.elig_annotated.astype(bool)]
F44 = analysis_set(fe, "factorial")
assert set(DROP) <= set(F44.cell_iname), "BEN, HCC1588 and MCF7 must be in the 44-line set"
F41 = analysis_set(fe[~fe.cell_iname.isin(DROP)], "factorial")
A25 = analysis_set(ae, "annotated")
sets = {"F44": (F44, "factorial"), "F41": (F41, "factorial"), "A25": (A25, "annotated")}

# ---- check: same sets from state-free files (eligible lines, per-lineage counts, TP53_cov) ----
el = pd.read_csv(f"{IN}/D_eligible_lines.tsv", sep="\t")
dc = pd.read_csv(f"{IN}/D_lineage_counts.tsv", sep="\t")
checks = []
for name, contrast in [("F44", "factorial"), ("A25", "annotated")]:
    c = dc[(dc.panel == "P1") & (dc.contrast == contrast)].set_index("lineage")
    gcols = ["n_g1", "n_g2", "n_g3", "n_g4"] if contrast == "factorial" else ["n_g1", "n_g2"]
    if contrast == "factorial":
        yv = (c.n_g1 + c.n_g4 > 0) & (c.n_g2 + c.n_g3 > 0)
        xv = (c.n_g2 + c.n_g4 > 0) & (c.n_g1 + c.n_g3 > 0)
        vary = yv | xv
    else:
        vary = (c.n_g1 > 0) & (c.n_g2 > 0)
    lins = sorted(c.index[(c.n_eligible >= 2) & vary])
    e = el[(el.panel == "P1") & (el.contrast == contrast)].merge(cov, on="cell_iname", how="left", validate="1:1")
    e = e[e.lineage.isin(lins)]
    assert e.TP53_cov.isin(KNOWN).all(), f"{name}: a line of a varying lineage lacks a TP53_cov call"
    s, _ = sets[name]
    same = set(e.cell_iname) == set(s.cell_iname)
    cs = counts(s.assign(_g=codes(s, contrast)), contrast, name).set_index("lineage")
    same_counts = all((cs.loc[l, g] == c.loc[l, g]) for l in lins for g in gcols) and sorted(cs.index) == lins
    checks.append(dict(set=name, check="line set equals the state-free derivation", passed=bool(same)))
    checks.append(dict(set=name, check="per-lineage counts equal D_lineage_counts.tsv", passed=bool(same_counts)))
    assert same and same_counts, name

lines_out, counts_out, summ = [], [], []
for name, (s, contrast) in sets.items():
    s = s.assign(_g=codes(s, contrast))
    lines_out.append(s[["cell_iname", "lineage", "TP53_cov"]].assign(set=name))
    cnt = counts(s, contrast, name)
    counts_out.append(cnt)
    r = dict(set=name, contrast=contrast, panel="P1", n_lines=len(s), n_lineages=s.lineage.nunique(),
             n_g1=int(cnt.n_g1.sum()), n_g2=int(cnt.n_g2.sum()))
    if contrast == "factorial":
        r.update(n_g3=int(cnt.n_g3.sum()), n_g4=int(cnt.n_g4.sum()),
                 n_Y_plus=int(cnt.n_g1.sum() + cnt.n_g4.sum()), n_X_plus=int(cnt.n_g2.sum() + cnt.n_g4.sum()))
    r["n_TP53_mutant"] = int((s.TP53_cov == "mutant").sum()); r["n_TP53_wild_type"] = int((s.TP53_cov == "wild_type").sum())
    r["lines_removed_from_F44"] = ";".join(DROP) if name == "F41" else ""
    summ.append(r)
    del s
lines_df = pd.concat(lines_out)[["set", "cell_iname", "lineage", "TP53_cov"]].sort_values(["set", "lineage", "cell_iname"])
assert set(lines_df.columns) == {"set", "cell_iname", "lineage", "TP53_cov"}  # no group column
lines_df.to_csv(f"{OUT}/D2_sets_lines.tsv", sep="\t", index=False)
pd.concat(counts_out).to_csv(f"{OUT}/D2_sets_counts.tsv", sep="\t", index=False)
sm = pd.DataFrame(summ)
sm.to_csv(f"{OUT}/D2_sets_summary.tsv", sep="\t", index=False)
checks.append(dict(set="all 129 lines", check="BROAD_ID and ccle_name equal in D_line_table.tsv and D_line_covariates.tsv (used by R4)",
                   passed=bool(ids_equal)))
assert ids_equal
pd.DataFrame(checks).to_csv(f"{OUT}/D2_sets_checks.tsv", sep="\t", index=False)
lost = sorted(set(F44.lineage) - set(F41.lineage))
with open(f"{ROOT}/logs/D2_runlog.txt", "a") as f:
    f.write(f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} [02_sets +{time.time()-t0:.0f}s] "
            + "; ".join(f"{r['set']}: {r['n_lines']} lines in {r['n_lineages']} lineages" for r in summ)
            + f"; lineages lost from F44 to F41: {lost if lost else 'none'}\n")
print(sm.to_string())
