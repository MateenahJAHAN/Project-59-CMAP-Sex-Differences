"""Analysis 1c (../PLAN.md, deviation 5): why the panel-level test is conservative under pseudo-labels.

On the first N pseudo-labellings of the power run (seed 20261131, the same draw), the frozen engine runs the panel-level
test for the Y and X terms with (b) lineage + Y + X (no TP53) and (d) lineage + the tested term alone. With lineage
blocks only, the Freedman-Lane permutation of reduced-model residuals within lineage is a relabelling of the lines, so
(d) is an exact test. The frozen design's values on the same pseudo-labellings are read from the N0 rows of
raw/power_runs.tsv.gz. Writes results/size_diagnostic.tsv and raw/size_diagnostic_rows.tsv.gz.
Usage: <python> -I size_diagnostic.py --N 200 --B 2000 --workers 8
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as K  # noqa: E402
import argparse
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

C, GG = K.C, K.GG
_W = {}


def _init():
    _W["D"] = K.Data(with_level3=False)


def variants(D, full_col, elig, term):
    P = D.P
    none = np.zeros((len(P.lines), 0))
    Xb, kb, tb = GG.xo_from_codes(P, none, full_col, "factorial", term)       # lineage + Y + X
    Xd = Xb[:, [kb]]                                                          # lineage + tested term
    return (("b: lineage + Y + X (no TP53)", Xb, kb, tb), ("d: lineage + tested term only", Xd, 0, [0]))


def _task(args):
    rep, full_col, elig, B = args
    D = _W["D"]
    out = []
    for term in K.TERMS:
        for name, Xo, k, tcols in variants(D, full_col, elig, term):
            mask = C.analysis_set(D.P.lin, Xo, tcols, elig)
            df, permT, recs = C.run_term(D.P, Xo, k, mask, B, rule=D.rule)
            dfF, obs, pp, _ = C.finish_table(df, permT)
            t = dfF[dfF.in_ranking]
            out.append(dict(rep=rep, term=term, variant=name, B=B, n_lines=int(mask.sum()), n_tested=len(t), panel_stat=obs,
                            panel_p=pp, n_p_lt_005=int((t.p < 0.05).sum()), n_p_lt_0025=int((t.p < 0.025).sum())))
    return out


def summarise(g, term, name, B):
    p = g.panel_p.values
    n = len(p)
    r = dict(term=term, variant=name, B=B, n_pseudo=n, mean_panel_p=float(p.mean()), mean_panel_p_mc_se=float(p.std(ddof=1) / np.sqrt(n)),
             median_panel_p=float(np.median(p)), compound_share_p_lt_005=float((g.n_p_lt_005 / g.n_tested).mean()))
    for a, lab in ((0.05, "005"), (0.025, "0025")):
        x = int((p < a).sum())
        lo, hi = K.clopper_pearson(x, n)
        r.update({f"n_p_lt_{lab}": x, f"share_p_lt_{lab}": x / n, f"cp95_lo_{lab}": lo, f"cp95_hi_{lab}": hi,
                  f"reading_{lab}": "conservative" if hi < a else ("liberal" if lo > a else "consistent with nominal")})
    lo, hi = r["mean_panel_p"] - 1.959964 * r["mean_panel_p_mc_se"], r["mean_panel_p"] + 1.959964 * r["mean_panel_p_mc_se"]
    r["mean_p_reading"] = "consistent with 0.5" if lo <= 0.5 <= hi else ("above 0.5" if lo > 0.5 else "below 0.5")
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=200)
    ap.add_argument("--B", type=int, default=2000)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    assert a.workers <= 8
    logp = os.path.join(K.LOGS, "size_diagnostic.log")
    t0 = time.time()
    K.log(logp, f"start: N {a.N}, B {a.B}, workers {a.workers}")
    runs = pd.read_csv(os.path.join(K.RAW, "power_runs.tsv.gz"), sep="\t", float_precision="round_trip")
    R = int(runs.rep.nunique())
    _init()
    D = _W["D"]
    pairs, cnt = K.factorial_set(D.P)
    full, elig = K.draw_codes(D.P, pairs, cnt, R, K.SEED_PSEUDO)     # the power run's draw (R columns)
    tasks = [(r, full[:, r].copy(), elig, a.B) for r in range(a.N)]
    with ProcessPoolExecutor(max_workers=a.workers, initializer=_init) as ex:
        rows = [x for out in ex.map(_task, tasks, chunksize=1) for x in out]
    raw = pd.DataFrame(rows).sort_values(["term", "variant", "rep"])
    raw.to_csv(os.path.join(K.RAW, "size_diagnostic_rows.tsv.gz"), sep="\t", index=False)
    write_summary(runs, raw, a.N, a.B)
    K.log(logp, f"done in {time.time() - t0:.0f} s")


def write_summary(runs, raw, N, B):
    out = []
    for term in K.TERMS:
        g0 = runs[(runs.term == term) & (runs.scenario == "N0") & (runs.rep < N)].sort_values("rep")
        assert int(g0.B.iloc[0]) == B
        r0 = summarise(g0, term, "a: frozen design (lineage + TP53 + Y + X)", B)
        r0.update(paired_diff_vs_frozen=0.0, paired_diff_se=0.0)
        out.append(r0)
        for name in raw.variant.unique():
            g = raw[(raw.term == term) & (raw.variant == name)].sort_values("rep")
            assert list(g.rep) == list(g0.rep)
            r = summarise(g, term, name, B)
            dd = g0.panel_p.values - g.panel_p.values            # PLAN.md deviation 9: same pseudo-labellings
            r.update(paired_diff_vs_frozen=float(dd.mean()), paired_diff_se=float(dd.std(ddof=1) / np.sqrt(len(dd))))
            out.append(r)
    pd.DataFrame(out).to_csv(os.path.join(K.RES, "size_diagnostic.tsv"), sep="\t", index=False)


def summary_only():
    """Rebuilds results/size_diagnostic.tsv from the saved rows (no new computation)."""
    runs = pd.read_csv(os.path.join(K.RAW, "power_runs.tsv.gz"), sep="\t", float_precision="round_trip")
    raw = pd.read_csv(os.path.join(K.RAW, "size_diagnostic_rows.tsv.gz"), sep="\t", float_precision="round_trip")
    write_summary(runs, raw, int(raw.rep.nunique()), int(raw.B.iloc[0]))


if __name__ == "__main__":
    if "--summary-only" in sys.argv:
        summary_only()
    else:
        main()
