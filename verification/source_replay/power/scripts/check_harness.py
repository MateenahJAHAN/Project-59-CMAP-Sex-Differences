"""Harness checks 1 to 4 and 6 of ../PLAN.md, run before the power run. Writes results/harness_checks.tsv.

1. Frozen code equals src/G/SHA256SUMS.txt.
2. Every input opened through ginputs.Inputs equals its frozen sha256.
3. This harness reproduces the frozen calibration rows (results/G/G_calibration.tsv, P1 factorial Y and X, configured
   runs, pseudo-labellings 0 to 3 of taskD2/pseudolabels/D2_pseudo_R1_F41.tsv.gz, B = 2,000).
4. E3 as scored here reproduces results/G/G_H3.tsv for P1 (factorial Y, factorial X, annotated; level 3 and level 5).
6. A C(6) and a G(6) run that recomputes only the planted compounds and reuses N0 for the others (power.py) equals a
   full run of gcore.run_term over every compound with the same planted data.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as K  # noqa: E402
import time

import numpy as np
import pandas as pd

C, GG, GS, GA = K.C, K.GG, K.GS, K.GA


def main():
    logp = os.path.join(K.LOGS, "check_harness.log")
    t_all = time.time()
    rows = K.check_frozen_code()
    K.log(logp, f"check 1: {sum(r['passes'] for r in rows)} of {len(rows)} frozen files match SHA256SUMS.txt")
    D = K.Data(with_level3=True)
    P = D.P
    # ---------------- check 3: calibration rows
    fl, _, cd = GG.read_pseudo(D.I.file("taskD2/pseudolabels/D2_pseudo_R1_F41.tsv.gz"))
    full, elig = GG.codes_full(P, fl, cd)
    prior = pd.read_csv(os.path.join(K.GRES, "G_calibration.tsv"), sep="\t", float_precision="round_trip")
    for term in K.TERMS:
        new = []
        for s in range(4):
            d = K.term_design(D, full[:, s], elig, term)
            df, permT, recs = K.run_tests(D, d, 2000)
            dfF, obs, pp, _ = K.finish(df, permT)
            t = dfF[dfF.tested]
            new.append(dict(pseudo=s, n_lines=int(d["mask"].sum()), n_g1=int(np.nansum(d["Xo"][d["mask"], d["k"]] == 1)),
                            n_g0=int(np.nansum(d["Xo"][d["mask"], d["k"]] == 0)), n_lineages=len(set(P.lin[d["mask"]])),
                            n_tested=len(t), n_p_lt_005=int((t.p < 0.05).sum()), n_q_lt_01=int((t.q < 0.1).sum()),
                            min_p=float(t.p.min()), panel_stat=obs, panel_p=pp))
        pr = prior[(prior.panel == "P1") & (prior.term == term) & (prior.run == "configured")]
        ok = GG.compare_rows(new, pr)
        rows.append(dict(check="3 calibration reproduced", item=f"P1 {term}, pseudo-labellings 0 to 3, B 2000",
                         expected="equal to results/G/G_calibration.tsv", observed=f"{int(ok.sum())} of {len(ok)} rows equal",
                         passes=bool(ok.all())))
        K.log(logp, f"check 3 {term}: {int(ok.sum())} of {len(ok)} calibration rows reproduced")
    # ---------------- check 4: E3 with the real labels
    lt = pd.read_csv(D.I.real_line_table_path(), sep="\t")
    al = GS.Aligned(P, lt)
    prim = GS.primary_settings(D.cfg)
    s5, s3 = D.s0()
    S = {"level5": np.nanmean(s5, axis=1), "level3": np.nanmean(s3, axis=1)}
    h3 = pd.read_csv(os.path.join(K.GRES, "G_H3.tsv"), sep="\t", float_precision="round_trip")
    for term in ("factY", "factX", "annotated"):
        d = al.design(term, prim)
        for scale, role, tn in (("level3", "primary", "R3_P1_d111"), ("level5", "S12", "P1_resp")):
            r = GA.e3_row(D.ctx, "P1", term, d, S[scale], scale, tn, role)
            ref = h3[(h3.panel == "P1") & (h3.term == term) & (h3.scale == scale) & (h3.role == role)].iloc[0]
            diffs = [abs(r[c] - ref[c]) / max(1e-300, abs(ref[c])) for c in ("estimate", "lo90", "hi90", "lo95", "hi95", "mean_S")]
            same_n = int(r["n_lines"]) == int(ref.n_lines) and int(r["df"]) == int(ref["df"])
            ok = bool(max(diffs) <= 1e-9 and same_n)
            rows.append(dict(check="4 E3 reproduced (real labels)", item=f"P1 {term} {scale} ({role})",
                             expected="equal to results/G/G_H3.tsv (relative 1e-9)", observed=f"max relative difference {max(diffs):.3g}; n and df equal {same_n}",
                             passes=ok))
    K.log(logp, f"check 4: {sum(1 for r in rows if r['check'].startswith('4') and r['passes'])} of 6 E3 rows reproduced")
    # ---------------- check 6: recomputing only the planted compounds equals a full run
    d = K.term_design(D, full[:, 0], elig, "factY")
    g1 = d["g1"]
    df0, permT0, recs0 = K.run_tests(D, d, 2000)
    tested = [r["_j"] for r in recs0 if r["tested"]]
    sig0 = {r["_j"]: r["_sig"] for r in recs0 if r["tested"]}
    rng = np.random.default_rng(12345)
    for scen in ("C6", "G6"):
        js = sorted(int(x) for x in rng.choice(tested, size=6, replace=False))
        module = np.sort(rng.choice(len(P.genes), size=K.MODULE_SIZE, replace=False))
        Y5 = P.R[:, js, :].copy()
        if scen == "C6":
            Y5[g1] *= K.C_FACTOR
        else:
            for a, j in enumerate(js):
                Y5[np.ix_(g1, [a], module)] += sig0[j][module][None, None, :]
        Yov = {j: Y5[:, a, :] for a, j in enumerate(js)}
        dfa, pa, ra = K.run_tests(D, d, 2000, Y_override=Yov)                 # full run over every compound
        A = K.finish(dfa, pa)
        df1, p1, r1 = K.run_tests(D, d, 2000, Y_override=Yov, comps=js)     # planted compounds only
        dfs, ps, rs = K.splice(df0, permT0, recs0, df1, p1, r1)
        Bf = K.finish(dfs, ps)
        cols = ["compound", "tested", "T_obs", "p", "q", "rank", "n_ranked", "in_ranking"]
        same_tab = K.frames_equal(A[0][cols], Bf[0][cols])
        same_perm = bool(np.array_equal(pa, ps, equal_nan=True))
        same_panel = bool(A[1] == Bf[1] and A[2] == Bf[2])
        rows.append(dict(check="6 recomputing only the planted compounds equals a full run", item=f"{scen}, P1 factY, pseudo-labelling 0, B 2000",
                         expected="identical tables, permutation statistics and panel p",
                         observed=f"table {same_tab}; permutation statistics {same_perm}; panel statistic and p {same_panel}",
                         passes=bool(same_tab and same_perm and same_panel)))
        K.log(logp, f"check 6 {scen}: table {same_tab}, permutations {same_perm}, panel {same_panel}")
    # ---------------- check 2: inputs
    ht = D.I.hash_table()
    accepted = ht.status.eq("equal")
    if "P59_LINCS_WORK" in os.environ:
        # ginputs verifies these separately pinned raw-rebuild hashes before assigning this status.
        accepted |= ht.status.eq("equal (alternative)") & ht.object.isin(
            ["taskD/D_line_table.tsv", "taskD/D_line_covariates.tsv"])
    rows.append(dict(check="2 input hashes", item=f"{len(ht)} files opened through ginputs.Inputs",
                     expected="registered or independently pinned reconstruction digest",
                     observed=f"{int(accepted.sum())} of {len(ht)} matched; see harness_inputs_sha256.tsv",
                     passes=bool(accepted.all())))
    ht.to_csv(os.path.join(K.RES, "harness_inputs_sha256.tsv"), sep="\t", index=False)
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(K.RES, "harness_checks.tsv"), sep="\t", index=False)
    K.log(logp, f"done in {time.time() - t_all:.0f} s: {int(out.passes.sum())} of {len(out)} checks pass")
    if not out.passes.all():
        raise SystemExit("a harness check failed: " + "; ".join(out[~out.passes].item))


if __name__ == "__main__":
    main()
