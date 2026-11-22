"""Planted-effect power and interval coverage on the final Part B design (post hoc; ../PLAN.md, analysis 2).

One replicate = one pseudo-labelling (task D's scheme, the real per-lineage counts) with every scenario of the plan,
for the factorial Y and X terms of P1. Effects are planted into the real response tensors of the pseudo-group-1 lines;
the frozen engine (gcore.run_term, gcore.finish_table) runs the tests, and E3 and E2 are scored with ganalyses.e3_row
and gcore.scalar_fit2 + ganalyses.express, from response sizes rebuilt by gcore.response_size in every replicate.

C(k) and G(k) change only k compounds. Their runs recompute those compounds with gcore.run_term(comps=...) and reuse
the replicate's N0 rows for every other compound (each compound's statistic and permutation column depend only on its
own data, design and the shared keys); harness check 6 shows this equals a full run.

Usage: <python> -I power.py --B <B> --R <R> --workers <W> --tag <name>
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as K  # noqa: E402  (sets the thread variables before numpy loads)
import argparse
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

C, GG = K.C, K.GG
_W = {}


def _init():
    D = K.Data(with_level3=True)
    D.s0()
    D.norm5, D.weak5, D.strong5 = K.template_classes(D.P.R)
    D.norm3, D.weak3, D.strong3 = K.template_classes(D.R3)
    _W["D"] = D


def e3_block(D, d, s5, s3):
    out = {}
    for sc, s, weak, strong in (("l3", s3, D.weak3, D.strong3), ("l5", s5, D.weak5, D.strong5)):
        scale = "level3" if sc == "l3" else "level5"
        allc = np.arange(s.shape[1])
        out.update(K.e3_cols(D, d, K.mean_over(s, allc), scale, f"e3{sc}"))
        out.update(K.e3_cols(D, d, K.mean_over(s, weak), scale, f"e3{sc}w"))
        out.update(K.e3_cols(D, d, K.mean_over(s, strong), scale, f"e3{sc}s"))
    return out


def e2_array(e2d, n_comp):
    a = np.full((n_comp, 5), np.nan)
    for j, v in e2d.items():
        a[j] = v
    return a


def replicate(D, rep, full_col, elig, B):
    """All scenarios of one replicate for both terms. Returns (rows, planted rows, {(term, scen): E2 array})."""
    P = D.P
    nC = len(P.comps)
    s5_0, s3_0 = D.s0()
    rows, prow, e2arr = [], [], {}
    for term in K.TERMS:
        d = K.term_design(D, full_col, elig, term)
        g1 = d["g1"]
        base = dict(rep=rep, term=term, B=B, n_lines=int(d["mask"].sum()), n_g1=len(g1))
        # ---------------- N0
        t0 = time.time()
        df0, permT0, recs0 = K.run_tests(D, d, B)
        dfF0, obs0, pp0, _ = K.finish(df0, permT0)
        tested = [r["_j"] for r in recs0 if r["tested"]]
        sig0 = {r["_j"]: r["_sig"] for r in recs0 if r["tested"]}
        e20 = K.e2_for(D, d, recs0, s5_0)
        e2arr[(term, "N0")] = e2_array(e20, nC)
        rows.append(dict(base, scenario="N0", panel_stat=obs0, panel_p=pp0, **K.detect(D, dfF0, []),
                         **e3_block(D, d, s5_0, s3_0), time_s=time.time() - t0))
        # ---------------- M(d): pooled magnitude, every compound, both scales
        for scen, delta in K.M_DELTA.items():
            t0 = time.time()
            R5p = P.R.copy(); R5p[g1] *= 1.0 + delta
            R3p = D.R3.copy(); R3p[g1] *= 1.0 + delta
            Yov = {j: R5p[:, j, :] for j in range(nC)}
            df, permT, recs = K.run_tests(D, d, B, Y_override=Yov)
            dfF, obs, pp, _ = K.finish(df, permT)
            s5 = C.response_size(R5p, D.avail5)
            s3 = C.response_size(R3p, D.avail3)
            e2d = K.e2_for(D, d, recs, s5)
            e2arr[(term, scen)] = e2_array(e2d, nC)
            rows.append(dict(base, scenario=scen, delta=delta, panel_stat=obs, panel_p=pp,
                             **K.detect(D, dfF, list(P.comps)), **e3_block(D, d, s5, s3), time_s=time.time() - t0))
            del R5p, R3p, Yov
        # ---------------- C(k): 50% larger responses in k compounds, both scales; G(k): +1 residual SD in a module
        for scen in list(K.C_K) + list(K.G_K):
            t0 = time.time()
            rng = np.random.default_rng([K.SEED_SCEN, rep, K.TERM_CODE[term], K.SCEN_CODE[scen]])
            k = K.C_K.get(scen) or K.G_K.get(scen)
            js = sorted(int(x) for x in rng.choice(tested, size=k, replace=False))
            module = np.sort(rng.choice(len(P.genes), size=K.MODULE_SIZE, replace=False)) if scen in K.G_K else None
            Y5 = P.R[:, js, :].copy()
            Y3 = D.R3[:, js, :].copy()
            if scen in K.C_K:
                Y5[g1] *= K.C_FACTOR
                Y3[g1] *= K.C_FACTOR
            else:
                for a, j in enumerate(js):
                    Y5[np.ix_(g1, [a], module)] += sig0[j][module][None, None, :]
            Yov = {j: Y5[:, a, :] for a, j in enumerate(js)}
            df1, permT1, recs1 = K.run_tests(D, d, B, Y_override=Yov, comps=js)
            df, permT, recs = K.splice(df0, permT0, recs0, df1, permT1, recs1)
            dfF, obs, pp, _ = K.finish(df, permT)
            s5 = s5_0.copy()
            s5[:, js] = C.response_size(Y5, D.avail5[:, js])
            s3 = s3_0.copy()
            if scen in K.C_K:
                s3[:, js] = C.response_size(Y3, D.avail3[:, js])
            e2d = K.e2_for(D, d, recs, s5, js=set(js))
            names = [P.comps[j] for j in js]
            rows.append(dict(base, scenario=scen, k=k, panel_stat=obs, panel_p=pp, **K.detect(D, dfF, names),
                             **e3_block(D, d, s5, s3), time_s=time.time() - t0))
            t = dfF.set_index("compound")
            rmap = {r["_j"]: r for r in recs if r["tested"]}
            for j, nm in zip(js, names):
                e1m = np.nan
                if module is not None and j in rmap:
                    r = rmap[j]
                    with np.errstate(invalid="ignore", divide="ignore"):
                        e1m = float(np.mean(r["_b"][module] / r["_sig"][module]))
                v, v0 = e2d.get(j, (np.nan,) * 5), e20.get(j, (np.nan,) * 5)
                prow.append(dict(rep=rep, term=term, scenario=scen, compound=nm, j=j, T_obs=t.at[nm, "T_obs"],
                                 p=t.at[nm, "p"], q=t.at[nm, "q"], rank=t.at[nm, "rank"], n_ranked=t.at[nm, "n_ranked"],
                                 e1_module_mean=e1m, e2_est=v[0], e2_lo90=v[1], e2_hi90=v[2], e2_lo95=v[3], e2_hi95=v[4],
                                 e2_N0_est=v0[0],
                                 module=";".join(str(int(x)) for x in module) if module is not None else ""))
        # ---------------- A(u): M(0.20) with u labelled group-1 lines carrying no effect (E3 only)
        for scen, u in K.A_U.items():
            t0 = time.time()
            rng = np.random.default_rng([K.SEED_SCEN, rep, K.TERM_CODE[term], K.SCEN_CODE[scen]])
            unaff = np.sort(rng.choice(g1, size=u, replace=False))
            g1e = np.setdiff1d(g1, unaff)
            R5p = P.R.copy(); R5p[g1e] *= 1.0 + K.A_DELTA
            R3p = D.R3.copy(); R3p[g1e] *= 1.0 + K.A_DELTA
            s5 = C.response_size(R5p, D.avail5)
            s3 = C.response_size(R3p, D.avail3)
            rows.append(dict(base, scenario=scen, u=u, delta=K.A_DELTA,
                             unaffected=";".join(P.lines[i] for i in unaff), **e3_block(D, d, s5, s3),
                             time_s=time.time() - t0))
            del R5p, R3p
    return rows, prow, e2arr


def _task(args):
    rep, full_col, elig, B = args
    t0 = time.time()
    out = replicate(_W["D"], rep, full_col, elig, B)
    return rep, out, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--B", type=int, required=True)
    ap.add_argument("--R", type=int, required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--tag", default="power")
    a = ap.parse_args()
    assert a.workers <= 8
    logp = os.path.join(K.LOGS, f"{a.tag}.log")
    K.log(logp, f"start: B {a.B}, R {a.R}, workers {a.workers}, pseudo-label seed {K.SEED_PSEUDO}")
    t_all = time.time()
    _init()
    D = _W["D"]
    pairs, cnt = K.factorial_set(D.P)
    full, elig = K.draw_codes(D.P, pairs, cnt, a.R, K.SEED_PSEUDO)
    tasks = [(r, full[:, r].copy(), elig, a.B) for r in range(a.R)]
    res = {}
    with ProcessPoolExecutor(max_workers=a.workers, initializer=_init) as ex:
        for n, (rep, out, dt) in enumerate(ex.map(_task, tasks, chunksize=1)):
            res[rep] = out
            if (n + 1) % 20 == 0 or n + 1 == a.R:
                K.log(logp, f"{n + 1} of {a.R} replicates done; last {dt:.1f} s; elapsed {time.time() - t_all:.0f} s")
    rows = [x for r in range(a.R) for x in res[r][0]]
    prow = [x for r in range(a.R) for x in res[r][1]]
    runs = pd.DataFrame(rows)
    planted = pd.DataFrame(prow)
    runs.to_csv(os.path.join(K.RAW, f"{a.tag}_runs.tsv.gz"), sep="\t", index=False)
    planted.to_csv(os.path.join(K.RAW, f"{a.tag}_planted.tsv.gz"), sep="\t", index=False)
    nC = len(D.P.comps)
    E2 = np.full((a.R, len(K.TERMS), len(K.E2_SCEN), nC, 5), np.nan, dtype=np.float64)
    for r in range(a.R):
        for ti, term in enumerate(K.TERMS):
            for si, sc in enumerate(K.E2_SCEN):
                E2[r, ti, si] = res[r][2][(term, sc)]
    np.savez_compressed(os.path.join(K.RAW, f"{a.tag}_e2.npz"), E2=E2, terms=np.array(K.TERMS),
                        scenarios=np.array(K.E2_SCEN), compounds=np.array(D.P.comps),
                        values=np.array(["est", "lo90", "hi90", "lo95", "hi95"]),
                        norm5=D.norm5, weak5=D.weak5, strong5=D.strong5, norm3=D.norm3, weak3=D.weak3, strong3=D.strong3)
    K.log(logp, f"pool done in {time.time() - t_all:.0f} s; {len(runs)} run rows, {len(planted)} planted rows")
    # check 5: replicate 0 again in this process must equal the pool's rows
    t0 = time.time()
    r0 = replicate(D, 0, full[:, 0].copy(), elig, a.B)
    a0 = pd.DataFrame(r0[0]).drop(columns=["time_s"])
    b0 = runs[runs.rep == 0].drop(columns=["time_s"])
    same_runs = K.frames_equal(a0, b0[a0.columns]) and list(a0.columns) == [c for c in b0.columns]
    same_pl = K.frames_equal(pd.DataFrame(r0[1]), planted[planted.rep == 0])
    same_e2 = all(np.array_equal(r0[2][(t, s)], E2[0, ti, si], equal_nan=True)
                  for ti, t in enumerate(K.TERMS) for si, s in enumerate(K.E2_SCEN))
    chk = pd.DataFrame([dict(check="5 determinism: replicate 0 recomputed in the main process", item=f"{a.tag}",
                             expected="equal to the pool's rows", observed=f"runs {same_runs}; planted {same_pl}; E2 {same_e2}",
                             passes=bool(same_runs and same_pl and same_e2))])
    chk.to_csv(os.path.join(K.RES, f"harness_check5_{a.tag}.tsv"), sep="\t", index=False)
    K.log(logp, f"check 5 (replicate 0 recomputed, {time.time() - t0:.0f} s): runs {same_runs}, planted {same_pl}, E2 {same_e2}")
    K.log(logp, f"end; total {time.time() - t_all:.0f} s")


if __name__ == "__main__":
    main()
