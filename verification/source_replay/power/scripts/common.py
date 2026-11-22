"""Shared code for the Part B sensitivity analyses (post hoc, added at revision v2; see ../PLAN.md).

Imports the frozen task G modules from the repository unchanged (src/G, checked against src/G/SHA256SUMS.txt) and
reads every input through ginputs.Inputs, which checks each file's sha256. Nothing is written to the repository.
"""
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import hashlib
import time
import types

REPO = os.environ["P59_G_REPO"]
RUN_DIR = os.environ["P59_G_RUN"]
INPUT_ROOT = os.path.join(RUN_DIR, "inputs", "data", "G")
OVERRIDES = {"inputs": os.path.join(INPUT_ROOT, "inputs"), "tensors_D": os.path.join(INPUT_ROOT, "tensors", "D"),
             "tensors_E": os.path.join(INPUT_ROOT, "tensors", "E")}
CAPTURE = os.path.join(RUN_DIR, "capture")
HERE = os.path.dirname(os.path.abspath(__file__))
ADIR = os.path.dirname(HERE)
RES, RAW, LOGS = (os.path.join(ADIR, x) for x in ("results", "raw", "logs"))
GRES = os.path.join(REPO, "results", "G")
sys.path.insert(0, os.path.join(REPO, "src", "G"))

import numpy as np
import pandas as pd
import gconfig
import ginputs
import gcore as C
import ggates as GG
import gsets as GS
import ganalyses as GA

assert os.environ.get("PYTHONHASHSEED") == "0", "run with PYTHONHASHSEED=0, as the real run"
for _m in (gconfig, ginputs, C, GG, GS, GA):
    assert os.path.dirname(os.path.abspath(_m.__file__)) == os.path.join(REPO, "src", "G"), _m.__file__

TERMS = ("factY", "factX")
TERM_CODE = {"factY": 0, "factX": 1}
SEED_PSEUDO = 20261131      # pseudo-labels of the power run (PLAN.md)
SEED_SCEN = 20261132        # planted compounds, gene modules and unaffected lines (PLAN.md)
SEED_TIMING = 20261199      # throwaway timing replicates (PLAN.md)
MODULE_SIZE = 49            # 5% of 978 landmark genes
SCENARIOS = (["N0"] + [f"M{int(round(d * 100))}" for d in (0.10, 0.20, 0.30, 0.50)]
             + [f"C{k}" for k in (1, 6, 23)] + [f"G{k}" for k in (1, 6, 23)] + ["A1", "A2"])
SCEN_CODE = {s: i for i, s in enumerate(SCENARIOS)}
M_DELTA = {"M10": 0.10, "M20": 0.20, "M30": 0.30, "M50": 0.50}
C_K = {"C1": 1, "C6": 6, "C23": 23}
G_K = {"G1": 1, "G6": 6, "G23": 23}
A_U = {"A1": 1, "A2": 2}
C_FACTOR = 1.5               # per-compound magnitude: a 50% larger response
A_DELTA = 0.20               # attenuation scenarios use M(0.20)
E2_SCEN = ("N0", "M10", "M20", "M30", "M50")   # scenarios whose E2 is kept for every compound


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def log(path, msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}"
    print(line, flush=True)
    with open(path, "a") as f:
        f.write(line + "\n")


def check_frozen_code():
    """Every file in src/G/SHA256SUMS.txt against its listed sha256."""
    rows = []
    for ln in open(os.path.join(REPO, "src", "G", "SHA256SUMS.txt")):
        if ln.startswith("#") or not ln.strip():
            continue
        sha, rel = ln.split(None, 1)
        rel = rel.strip()
        got = sha256_file(os.path.join(REPO, rel))
        rows.append(dict(check="frozen code", item=rel, expected=sha, observed=got, passes=got == sha))
    return rows


class Data:
    """P1 (level 5, task D) and its level-3 tensor R3_P1_d111 (task E) in task D's axis order, via the frozen loaders."""

    def __init__(self, with_level3=True):
        self.cfg = gconfig.load(os.path.join(REPO, "config", "G_defaults.json"))
        self.I = ginputs.Inputs(REPO, self.cfg, OVERRIDES, strict=True)
        self.cov = self.I.covariates()
        self.P = self.I.panel("P1", self.cov, self.cfg.run("B_full"))
        self.R3 = self.I.level3("R3_P1_d111", self.P, order="D")[0] if with_level3 else None
        self.tpX = GS.tp53_columns(self.P.cov, GS.tp53_code(self.cfg))[0]
        self.rule = GS.rule_for(self.cfg)
        self.ctx = types.SimpleNamespace(cfg=self.cfg)
        self.avail5 = self.P.avail
        self.avail3 = ~np.isnan(self.R3[:, :, 0]) if with_level3 else None
        self._s0 = None

    def s0(self):
        """Unplanted response sizes (label-free): s on level 5 and level 3, lines x compounds."""
        if self._s0 is None:
            self._s0 = (C.response_size(self.P.R, self.avail5), C.response_size(self.R3, self.avail3))
        return self._s0


def factorial_set(P):
    """The 41-line P1 factorial set and its per-lineage group counts (codes 1 Y+ only, 2 X+ only, 3 neither)."""
    sl = pd.read_csv(os.path.join(GRES, "G_sets_lines.tsv"), sep="\t")
    sl = sl[(sl.panel == "P1") & (sl.family == "factorial")]
    cnt = pd.read_csv(os.path.join(GRES, "G_sets_counts.tsv"), sep="\t")
    cnt = cnt[(cnt.panel == "P1") & (cnt.family == "factorial")].drop(columns=["panel", "family"]).reset_index(drop=True)
    lin_of = dict(zip(P.lines, P.lin))
    assert all(lin_of[a] == b for a, b in zip(sl.cell_iname, sl.lineage))
    pairs = list(zip(sl.cell_iname, sl.lineage))
    assert int(cnt.n_lines.sum()) == len(pairs) and int(cnt.n_g4.sum()) == 0
    return pairs, cnt


def draw_codes(P, pairs, cnt, R, seed):
    """Task D's pseudo-label scheme (ggates.fresh_codes) with the real per-lineage counts; codes on the panel's lines."""
    flines, _, codes = GG.fresh_codes(pairs, cnt, R, seed)
    full, elig = GG.codes_full(P, flines, codes)
    return full, elig


def term_design(D, full_col, elig, term):
    """The engine design of one pseudo-labelling, as ggates.calib_rows builds it."""
    P = D.P
    Xo, k, tcols = GG.xo_from_codes(P, D.tpX, full_col, "factorial", term)
    mask = C.analysis_set(P.lin, Xo, tcols, elig)
    g1 = np.flatnonzero(mask & (Xo[:, k] == 1))
    return dict(term=term, lin=P.lin, Xo=Xo, k=k, tcols=tcols, mask=mask, rule=D.rule, g1=g1)


def run_tests(D, d, B, Y_override=None, comps=None):
    """The frozen per-compound tests and panel-level test (gcore.run_term + gcore.finish_table)."""
    df, permT, recs = C.run_term(D.P, d["Xo"], d["k"], d["mask"], B, comps=comps, Y_override=Y_override, rule=d["rule"])
    return df, permT, recs


def finish(df, permT):
    return C.finish_table(df, permT)


def splice(df0, permT0, recs0, df1, permT1, recs1):
    """N0's table with the rows and permutation columns of the recomputed (planted) compounds replaced."""
    df = df0.copy()
    permT = permT0.copy()
    recs = list(recs0)
    pos = {c: i for i, c in enumerate(df0.compound)}
    for i1, c in enumerate(df1.compound):
        i0 = pos[c]
        for col in df1.columns:
            df.at[i0, col] = df1.at[i1, col]
        permT[:, i0] = permT1[:, i1]
        recs[i0] = recs1[i1]
    return df, permT, recs


def detect(D, dfF, planted_names):
    """Detection readouts from a finished table (rules 1 to 3 of H2; rule 4 is not applied)."""
    cfg = D.cfg
    t = dfF[dfF.in_ranking]
    cut = int(np.floor(cfg.const("H2_top_fraction") * len(t)))
    pl = t.compound.isin(planted_names).values
    q, p, rk = t.q.values, t.p.values, t["rank"].values
    n_p05 = int((p < cfg.const("E4_count_p")).sum())
    out = dict(n_tested=int(len(t)), n_q_lt_01=int((q < cfg.const("H2_bh_q")).sum()), n_p_lt_005=n_p05,
               n_p_lt_0025=int((p < 0.025).sum()), rule3_E4=bool(n_p05 > cfg.const("E4_threshold_7kA")),
               n_planted=int(pl.sum()), n_planted_q_lt_01=int(((q < cfg.const("H2_bh_q")) & pl).sum()),
               n_planted_p_lt_005=int(((p < 0.05) & pl).sum()),
               n_planted_rule12=int(((q < cfg.const("H2_bh_q")) & (rk <= cut) & pl).sum()),
               n_other_q_lt_01=int(((q < cfg.const("H2_bh_q")) & ~pl).sum()), top5pct_cut=cut)
    return out


def e3_cols(D, d, S, scale, prefix):
    """E3 exactly as ganalyses.e3_row computes it, flattened into columns with a prefix."""
    r = GA.e3_row(D.ctx, "P1", d["term"], d, S, scale, "", "sensitivity")
    keys = ("estimate", "lo90", "hi90", "lo95", "hi95", "n_lines", "df", "mean_S", "difference_declared")
    return {f"{prefix}_{k}": r.get(k, np.nan) for k in keys}


def e2_for(D, d, recs, s5, js=None):
    """E2 per tested compound on level 5, as ganalyses.run_main_term computes it, at 90% and 95%.
    Returns {j: (est, lo90, hi90, lo95, hi95)}."""
    out = {}
    for r in recs:
        if not r["tested"] or (js is not None and r["_j"] not in js):
            continue
        f = C.scalar_fit2(d["lin"], d["Xo"], d["k"], r["_rows"], s5[:, r["_j"]], rule=d["rule"])
        est, lo90, hi90 = GA.express(f, D.cfg, 0.90)
        _, lo95, hi95 = GA.express(f, D.cfg, 0.95)
        out[r["_j"]] = (est, lo90, hi90, lo95, hi95)
    return out


def template_classes(R):
    """Label-free template strength per compound: |mean response over the available panel lines|; returns
    (norms, weak indices = the 29 smallest, strong indices = the 29 largest) over all compounds of the panel."""
    with np.errstate(invalid="ignore"):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            m = np.nanmean(R, axis=0)
    norms = np.linalg.norm(m, axis=1)
    nq = int(round(0.25 * len(norms)))
    o = np.argsort(norms, kind="stable")
    return norms, np.sort(o[:nq]), np.sort(o[-nq:])


def mean_over(s, cols):
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(s[:, cols], axis=1)


def frames_equal(a, b):
    """Exact equality of two tables: same columns and rows; numbers equal bit for bit (NaN equals NaN)."""
    a, b = a.reset_index(drop=True), b.reset_index(drop=True)
    if list(a.columns) != list(b.columns) or len(a) != len(b):
        return False
    for c in a.columns:
        if a[c].dtype.kind in "fiub" and b[c].dtype.kind in "fiub":
            if not np.array_equal(a[c].values.astype(float), b[c].values.astype(float), equal_nan=True):
                return False
        else:
            sa = a[c].astype(object).where(a[c].notna(), "<NA>").astype(str).values
            sb = b[c].astype(object).where(b[c].notna(), "<NA>").astype(str).values
            if not (sa == sb).all():
                return False
    return True


def wilson(x, n, z=1.959963984540054):
    if n == 0:
        return np.nan, np.nan
    p = x / n
    den = 1 + z * z / n
    c = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, c - h), min(1.0, c + h)


def clopper_pearson(x, n, a=0.05):
    from scipy import stats
    lo = 0.0 if x == 0 else float(stats.beta.ppf(a / 2, x, n - x + 1))
    hi = 1.0 if x == n else float(stats.beta.ppf(1 - a / 2, x + 1, n - x))
    return lo, hi
