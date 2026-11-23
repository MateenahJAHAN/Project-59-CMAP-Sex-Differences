#!/usr/bin/env python3
"""Task G entry point: the single Part B unblinding run, in four steps.

  1. guard   compares the given line table with the frozen real one (guard.py). Real labels need --mode real,
             --plan with P59_UNBLIND_PLAN_SHA256 equal to its sha256, code equal to src/G/SHA256SUMS.txt (7m),
             --plan-tag (8b) and full settings. Scrambled labels need --mode test.
  2. sets    per panel and contrast, the analysis-set membership (no group column) and per-lineage counts.
  3. gates   label-free: PC1 (level 5), K4 (level 3), pseudo-label calibration with replay of tasks D and D2,
             label-free E3 precision, DMSO gate (7b C). A real run stops if a gate fails or a frozen value is not
             reproduced.
  4. merge   joins states and labels to responses (ganalyses.py and the modules it calls) and writes results.

Usage (test, scrambled labels):
  python src/G/run_g.py --mode test --line-table <scrambled.tsv> --out <dir> [--B 200 --calib-n 4 ...]
Usage (real, at unblinding only; see src/G/README.md):
  P59_UNBLIND_PLAN_SHA256=<sha256 of the frozen plan> python src/G/run_g.py --mode real \\
      --line-table data/G/inputs/taskD/D_line_table.tsv --plan analysis_plan/G_plan.md --plan-tag <tag> --out results/G
No step prints a statistic of the merge; tables are written to --out.
"""
import os, sys
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
import argparse, time, platform, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
ROOT = os.path.dirname(os.path.dirname(HERE))
import numpy as np, pandas as pd
import gconfig, guard, ginputs, gcore as C, gsets as GS, ggates as GG, ganalyses as GA
import gsecondary as GSC, gsensitivity as GSN, gviability as GV, gdecisions as GD, gfigures as GF

TAG = "run_g"


class Ctx:
    pass


class PanelMeta:
    """Lines, lineages and covariates of a panel, without any response (used while labels are read for sets)."""

    def __init__(self, pn, I, cov):
        self.pn = pn
        self.lines = pd.read_csv(I.file(f"taskD/tensors/{pn}_lines.tsv"), sep="\t").cell_iname.tolist()
        c = pd.DataFrame({"cell_iname": self.lines}).merge(cov, on="cell_iname", how="left", validate="1:1")
        assert len(c) == len(self.lines) and c.lineage.notna().all()
        self.cov = c
        self.lin = np.array(c.lineage.astype(str).tolist(), dtype="U64")


def parse():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default="config/G_defaults.json")
    ap.add_argument("--mode", choices=["test", "real"], required=True)
    ap.add_argument("--line-table", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--plan"); ap.add_argument("--plan-tag"); ap.add_argument("--osf-id")
    for n in ("B", "calib-n", "calib-B", "replay-n", "relabel-n", "random-n"):
        ap.add_argument(f"--{n}", type=int, help="reduced setting (test runs only)")
    ap.add_argument("--workers", type=int)
    for n in ("inputs", "tensors-D", "tensors-E", "persig-E", "siginfo", "prism"):
        ap.add_argument(f"--{n}")
    ap.add_argument("--log", default="logs/G_runlog.txt")
    ap.add_argument("--stop-after", choices=["guard", "sets", "gates"])
    return ap.parse_args()


def w(out, name, df):
    df.to_csv(os.path.join(out, name), sep="\t", index=False)


def main(*, args=None, replay_check=None):
    a = parse() if args is None else args
    t0 = time.time()
    cfg = gconfig.load(gconfig.path(ROOT, a.config))
    os.makedirs(os.path.dirname(gconfig.path(ROOT, a.log)), exist_ok=True)
    C.set_log(gconfig.path(ROOT, a.log))
    out = gconfig.path(ROOT, a.out)
    if os.path.exists(out) and os.listdir(out):
        raise SystemExit(f"output folder {out} is not empty")
    reduced ={n: getattr(a, n) for n in ("B", "calib_n", "calib_B", "replay_n", "relabel_n", "random_n") if getattr(a, n) is not None}
    ov = {"inputs": a.inputs, "tensors_D": a.tensors_D, "tensors_E": a.tensors_E, "persig_E": a.persig_E, "siginfo": a.siginfo, "prism": a.prism}
    I = ginputs.Inputs(ROOT, cfg, {k: v for k, v in ov.items() if v}, strict=(a.mode == "real"))
    C.log(TAG, f"start: mode {a.mode}; config sha256 {cfg.sha256[:16]}; reduced settings {sorted(reduced)}")

    # ======================= 1. guard =======================
    real = pd.read_csv(I.real_line_table_path(), sep="\t")
    given = pd.read_csv(a.line_table, sep="\t")
    info = guard.classify(given, real, cfg.run("test_min_share_label_rows_changed"))
    del real
    if replay_check is None:
        gi = guard.check(info, a.mode, cfg, ROOT, a.plan, a.plan_tag, a.osf_id, list(reduced))
    else:
        gi = replay_check(info, a, cfg, reduced)
    C.log(TAG, f"guard: labels {info['status']} ({info['n_rows_changed']} of {info['n_label_rows']} label rows differ from the frozen table); run allowed as {a.mode}")
    os.makedirs(os.path.join(out, "figures"), exist_ok=True)
    manifest = [dict(item="mode", value=a.mode), dict(item="label_status", value=info["status"]),
                dict(item="label_rows_changed", value=f"{info['n_rows_changed']} of {info['n_label_rows']}"),
                dict(item="line_table", value=os.path.abspath(a.line_table)),
                dict(item="line_table_sha256", value=guard.sha256_file(a.line_table)),
                dict(item="config", value=os.path.relpath(cfg.file_path, ROOT)), dict(item="config_sha256", value=cfg.sha256),
                dict(item="python", value=platform.python_version()), dict(item="numpy", value=np.__version__),
                dict(item="pandas", value=pd.__version__), dict(item="platform", value=platform.platform()),
                dict(item="reduced_settings", value=";".join(f"{k}={v}" for k, v in sorted(reduced.items())))]
    manifest += [dict(item=f"guard_{k}", value=str(v)) for k, v in gi.items()]
    if a.stop_after == "guard":
        finish(out, cfg, manifest, I, t0); return

    B = a.B or cfg.run("B_full")
    calB = a.calib_B or cfg.run("calibration_B")
    relabel_n = a.relabel_n or cfg.run("relabel_n_viability")
    random_n = a.random_n or cfg.run("random_set_n")
    workers = a.workers or cfg.run("workers")
    panels = cfg.const("panels")
    terms = ["factY", "factX", "annotated"] + (["strict"] if cfg.decide("7c") in ("A", "C") else [])
    families = ["factorial", "annotated"] + (["strict"] if "strict" in terms else [])
    fam_of = {"factY": "factorial", "factX": "factorial", "annotated": "annotated", "strict": "strict"}
    prim = GS.primary_settings(cfg)
    tcode = GS.tp53_code(cfg)

    # ======================= 2. sets (labels -> membership and counts only) =======================
    cov = I.covariates()
    meta = {pn: PanelMeta(pn, I, cov) for pn in panels}
    sets, sl_rows, sc_rows = {}, [], []
    for pn in panels:
        al = GS.Aligned(meta[pn], given)
        for fam in families:
            lines, cnt, d = al.set_and_counts(fam, prim)
            sets[(pn, fam)] = (lines, cnt)
            lin_of = dict(zip(meta[pn].lines, meta[pn].lin))
            sl_rows += [dict(panel=pn, family=fam, cell_iname=l, lineage=lin_of[l]) for l in lines]
            sc_rows += [dict(panel=pn, family=fam, **r) for r in cnt.to_dict("records")]
    del given, al
    w(out, "G_sets_lines.tsv", pd.DataFrame(sl_rows)); w(out, "G_sets_counts.tsv", pd.DataFrame(sc_rows))
    # comparison with task D2's frozen sets (membership and counts)
    d2l, d2c = I.table("taskD2/D2_sets_lines.tsv"), I.table("taskD2/D2_sets_counts.tsv")
    frozen_name = {("P1", "factorial"): {"A": "F41", "B": "F44"}.get(cfg.decide("7f")), ("P1", "annotated"): "A25",
                   ("P1b", "factorial"): {"A": "F41", "B": "F44"}.get(cfg.decide("7f")), ("P1b", "annotated"): "A25"}
    chk = []
    for (pn, fam), nm in frozen_name.items():
        if nm is None or (pn, fam) not in sets:
            continue
        lines, cnt = sets[(pn, fam)]
        fl = d2l[d2l.set == nm]
        chk.append(dict(panel=pn, family=fam, frozen_set=f"task D2 {nm}", same_lines=sorted(lines) == sorted(fl.cell_iname),
                        same_counts=GG.same_counts(cnt, d2c[d2c.set == nm])))
    sets_chk = pd.DataFrame(chk)
    w(out, "G_sets_check.tsv", sets_chk)
    stops = []
    if a.mode == "real" and len(sets_chk) and not (sets_chk.same_lines & sets_chk.same_counts).all():
        stops.append("an analysis set differs from task D2's frozen set")
    C.log(TAG, f"sets: {len(sets)} analysis sets written (membership and counts only)")
    if a.stop_after == "sets":
        finish(out, cfg, manifest, I, t0); return

    # ======================= 3. gates (label-free) =======================
    P = {pn: I.panel(pn, cov, cfg.run("B_full")) for pn in panels}
    genes_E = I.e_genes()
    gates, calib_rows, calib_summ, replay_rows, prec_rows = [], [], [], [], []
    full_B = (B == cfg.run("B_full"))
    # PC1
    row, pc1_tab = GG.pc1(P["P1"], cfg, B)
    pr = I.table("taskD/D_pc1_summary.tsv")
    prior = pr[(pr.panel == "P1") & (pr.test == "TP53 (CMP)")].iloc[0].to_dict()
    row["reproduces_frozen"] = GG.pc1_reproduces(row, prior) if full_B else np.nan
    row["frozen_source"] = "task D results/D/D_pc1_summary.tsv, P1 / TP53 (CMP)"
    gates.append(row)
    C.log(TAG, f"gate PC1: passes {row['passes']}; reproduces task D {row['reproduces_frozen']}", t0)
    if a.mode == "real" and not row["passes"]:
        stops.append("PC1 fails (K2)")
    if a.mode == "real" and row["reproduces_frozen"] is not True:
        stops.append("PC1 does not reproduce task D")
    # K4 on the level-3 tensor of 7d (A or C: R3_P1_d111; B: R3_P1)
    l3obj = "R3_P1" if cfg.decide("7d") == "B" else "R3_P1_d111"
    z5obj = {"R3_P1": "Z5_P1", "R3_P1_d111": "Z5_P1_d111"}[l3obj]
    k4row = None
    R3E, Z5E = I.level3(l3obj, P["P1"], order="E"), I.level3(z5obj, P["P1"], order="E")
    if R3E is not None and Z5E is not None:
        k4row = GG.k4(R3E[0], Z5E[0], R3E[1], R3E[2], R3E[3], I.table("taskE/E_pc1_lines.tsv", dtype=str), cfg, B, l3obj)
        if l3obj == "R3_P1_d111":
            ks = I.table("taskD2/R3/D2_K4_summary.tsv")
            l3, l5 = ks[ks.analysis == "level3"].iloc[0], ks[ks.analysis == "level5"].iloc[0]
            p3 = dict(rank=l3["rank"], n_ranked=l3.n_compounds_ranked, p=l3.p_perm, T=l3.T_c, beta=l3.CDKN1A_coef_mut_minus_wt)
            p5 = dict(beta=l5.CDKN1A_coef_mut_minus_wt)
            k4row["frozen_source"] = "task D2 results/D2/R3/D2_K4_summary.tsv"
        else:
            es = I.table("taskE/E_pc1_summary.tsv")
            l3, l5 = es[es.analysis == "level3"].iloc[0], es[es.analysis == "level5"].iloc[0]
            p3 = dict(rank=l3["rank"], n_ranked=l3.n_compounds_ranked, p=l3.p_perm, T=l3.T_c, beta=l3.CDKN1A_coef_mut_minus_wt)
            p5 = dict(beta=l5.CDKN1A_coef_mut_minus_wt)
            k4row["frozen_source"] = "task E results/E/E_pc1_summary.tsv"
        k4row["reproduces_frozen"] = GG.k4_reproduces(k4row, p3, p5) if full_B else np.nan
        gates.append(k4row)
        C.log(TAG, f"gate K4 ({l3obj}): passes {k4row['passes']}; reproduces {k4row['reproduces_frozen']}", t0)
        if a.mode == "real" and k4row["reproduces_frozen"] is not True:
            stops.append("K4 does not reproduce its frozen value")
    else:
        C.log(TAG, f"gate K4: tensor {l3obj} not available in this run (test); K4 not run")
    k4_pass = bool(k4row is not None and k4row["passes"])
    # DMSO gate (7b C)
    if cfg.decide("7b") == "C":
        dp = row.get("dmso_p", np.nan)
        gates.append(dict(gate="DMSO", panel="P1", p=dp, passes=bool(np.isfinite(dp) and dp >= 0.05)))
        if a.mode == "real" and not gates[-1]["passes"]:
            stops.append("DMSO gate fails (7b C)")
    # calibration: configured runs, prior comparison and replays
    sources = cfg.raw["pseudo_sources"]["list"]
    prule = dict(min_group=cfg.const("min_group_7iB"), min_relab=0)
    for pn in panels:
        Pp = P[pn]
        panel_args = (pn, I.d_tensor(f"{pn}_resp.npy"), I.file(f"taskD/tensors/{pn}_lines.tsv"), I.file(f"taskD/tensors/{pn}_compounds.tsv"),
                      Pp.genes, cov, cfg.run("key_seed_base") + cfg.const("panel_numbers")[pn], cfg.run("B_full"), None)
        tpX, _ = GS.tp53_columns(Pp.cov, tcode)
        tpB, _ = GS.tp53_columns(Pp.cov, "7gB")
        for term in terms:
            fam = fam_of[term]
            lines, cnt = sets[(pn, fam)]
            n_want = a.calib_n or cfg.run("calibration_n_pseudo")[pn]
            src, codes, flines = None, None, None
            cands = sorted([s for s in sources if s["family"] == fam and pn in s["panels"]], key=lambda s: term not in s["prior"].get(pn, {}))
            for s in cands:
                fl, flin, cd = GG.read_pseudo(I.file(s["file"]))
                il, ic = GG.implied_set(meta[pn], tpB, fl, cd, fam)
                if sorted(il) == sorted(lines) and GG.same_counts(ic, cnt):
                    src, codes, flines = s, cd, fl
                    break
            if src is None:
                lin_of = dict(zip(meta[pn].lines, meta[pn].lin))
                seed = cfg.run("calibration_fresh_seed") + 10 * cfg.const("panel_numbers")[pn] + families.index(fam)
                flines, _, codes = GG.fresh_codes([(l, lin_of[l]) for l in lines], cnt, n_want, seed)
                src_name = f"fresh draws (seed {seed})"
            else:
                src_name = src["name"]
            n_use = min(n_want, codes.shape[1])
            full, elig = GG.codes_full(Pp, flines, codes)
            rows = GG.run_calibration(panel_args, Pp, tpX, full, elig, fam, term, calB, prim["rule"], range(n_use), workers)
            for r in rows:
                r.update(panel=pn, source=src_name, run="configured")
            prior_file = src["prior"].get(pn, {}).get(term) if src else None
            rep_note, reproduced = "no prior run with these lines", np.nan
            if prior_file and calB == cfg.raw["pseudo_sources"]["prior_B"]:
                prior = pd.read_csv(I.file(prior_file), sep="\t", float_precision="round_trip")
                # rows whose computation can differ from the prior run's: all rows if the TP53 coding differs from
                # task D's, else the rows in which the 7i A rule left a compound untested
                need = list(range(n_use)) if tcode != "7gB" else [r["pseudo"] for r in rows if r.get("n_untested_by_7iA_rule", 0) > 0]
                direct = [r for r in rows if r["pseudo"] not in need]
                ok_d = GG.compare_rows(direct, prior) if direct else pd.Series(dtype=bool)
                ok_r = pd.Series(dtype=bool)
                if need:
                    rr = GG.run_calibration(panel_args, Pp, tpB, full, elig, fam, term, calB, prule, need, workers)
                    for r in rr:
                        r.update(panel=pn, source=src_name, run="replay with the prior's settings")
                    replay_rows += rr
                    ok_r = GG.compare_rows(rr, prior)
                reproduced = bool(ok_d.all() and ok_r.all())
                rep_note = f"{len(direct)} rows compared directly, {len(need)} replayed with task D's rule; prior {prior_file}"
                if a.mode == "real" and not reproduced:
                    stops.append(f"calibration {pn} {term} does not reproduce {prior_file}")
            sm = GG.summarize(rows, cfg.const("calibration_gate"))
            calib_rows += rows
            calib_summ.append(dict(panel=pn, term=term, source=src_name, B=calB, **sm, reproduces_prior=reproduced, reproduction_note=rep_note))
            gates.append(dict(gate="calibration", panel=pn, term=term, passes=sm["gate_pass"], mean_share_p_lt_005=sm.get("mean_share_p_lt_005", np.nan),
                              n_pseudo=sm["n_pseudo"], reproduces_frozen=reproduced, frozen_source=prior_file or ""))
            C.log(TAG, f"gate calibration {pn} {term}: {src_name}, {n_use} pseudo-labellings; gate pass {sm['gate_pass']}; reproduces prior {reproduced}", t0)
            if a.mode == "real" and not sm["gate_pass"]:
                stops.append(f"calibration gate fails for {pn} {term} (K6)")
    # test runs: replay the first rows of every prior run (checks the reproduction path)
    if a.mode == "test":
        nrep = a.replay_n if a.replay_n is not None else 2
        for s in sources:
            for pn, tt in s["prior"].items():
                Pp = P[pn]
                panel_args = (pn, I.d_tensor(f"{pn}_resp.npy"), I.file(f"taskD/tensors/{pn}_lines.tsv"), I.file(f"taskD/tensors/{pn}_compounds.tsv"),
                              Pp.genes, cov, cfg.run("key_seed_base") + cfg.const("panel_numbers")[pn], cfg.run("B_full"), None)
                tpB, _ = GS.tp53_columns(Pp.cov, "7gB")
                fl, flin, cd = GG.read_pseudo(I.file(s["file"]))
                full, elig = GG.codes_full(Pp, fl, cd)
                for term, pf in tt.items():
                    if term not in terms or nrep <= 0:
                        continue
                    rr = GG.run_calibration(panel_args, Pp, tpB, full, elig, s["family"], term, cfg.raw["pseudo_sources"]["prior_B"], prule, range(nrep), workers)
                    prior = pd.read_csv(I.file(pf), sep="\t", float_precision="round_trip")
                    ok = GG.compare_rows(rr, prior)
                    for r in rr:
                        r.update(panel=pn, source=s["name"], run="test replay of the prior run", reproduces_prior=bool(ok.loc[r["pseudo"]]))
                    replay_rows += rr
                    C.log(TAG, f"test replay {s['name']} {pn} {term}: {int(ok.sum())} of {len(ok)} rows reproduce {pf}", t0)
    # label-free E3 precision (task D2's scheme) for the primary panel
    scale_E3 = "level5" if (cfg.decide("7d") == "C" or not k4_pass) else "level3"
    S_lab = {}
    if scale_E3 == "level3" and R3E is not None:
        R3 = R3E[0]
        S_lab["P1"] = (np.nanmean(C.response_size(R3, ~np.isnan(R3[:, :, 0])), axis=1), R3E[1], l3obj)
    else:
        S_lab["P1"] = (np.nanmean(C.response_size(P["P1"].R, P["P1"].avail), axis=1), P["P1"].lines, "P1_resp")
    S_all, plines, tname = S_lab["P1"]
    tps = dict(zip(P["P1"].cov.cell_iname, GS.tp53_strings(P["P1"].cov, tcode)))
    d2p = I.table("taskD2/R3/D2_R3_precision.tsv")
    pl_free = {}
    for fam in ("factorial", "annotated"):
        lines, cnt = sets[("P1", fam)]
        lin_of = dict(zip(meta["P1"].lines, meta["P1"].lin))
        rows = GG.precision_e3(S_all, plines, lines, [lin_of[l] for l in lines], [tps[l] for l in lines], cnt, fam == "factorial",
                               cfg.run("precision_draws"), cfg.run("precision_seed"), cfg.const("E3_margin_pct"))
        nm = frozen_name.get(("P1", fam))
        match = sets_chk[(sets_chk.panel == "P1") & (sets_chk.family == fam)]
        same_set = bool(len(match) and match.same_lines.iloc[0] and match.same_counts.iloc[0])
        for r in rows:
            q = d2p[(d2p.tensor == tname.replace("P1_resp", "Z5_P1_d111")) & (d2p.set == (nm or "")) & (d2p.term == r["term"])]
            rep = np.nan
            if same_set and len(q) and tcode == "7gB" and tname != "P1_resp":
                q = q.iloc[0]
                rep = bool(all(r[c] == q[c] for c in ("P_ci_inside_pm15pct", "P_ci_inside_pm20pct", "P_ci_inside_pm25pct",
                                                      "P_ci_inside_pm15pct_hc3", "P_ci_inside_pm20pct_hc3", "P_ci_inside_pm25pct_hc3")))
            col = "P_ci_inside_pm20pct_hc3" if cfg.decide("7j")["interval"] == "A" else "P_ci_inside_pm20pct"
            pl_free[r["term"]] = r[col]
            prec_rows.append(dict(panel="P1", family=fam, scale=scale_E3, tensor=tname, interval_used=col, **r, reproduces_D2=rep))
            if a.mode == "real" and rep is False:
                stops.append(f"label-free E3 precision {r['term']} does not reproduce task D2")
    gt = pd.DataFrame(gates)
    w(out, "G_gates.tsv", gt); w(out, "G_calibration.tsv", pd.DataFrame(calib_rows)); w(out, "G_calibration_summary.tsv", pd.DataFrame(calib_summ))
    w(out, "G_calibration_replay.tsv", pd.DataFrame(replay_rows)); w(out, "G_E3_precision_labelfree.tsv", pd.DataFrame(prec_rows))
    if stops:
        C.log(TAG, f"STOP before unblinding: {stops}")
        w(out, "G_decisions.tsv", pd.DataFrame([dict(section="5", item="run status", rule="gates must pass and frozen values must be reproduced",
                                                     outcome="stopped before any state was joined", consequence="; ".join(stops))]))
        finish(out, cfg, manifest, I, t0); return
    if a.stop_after == "gates":
        finish(out, cfg, manifest, I, t0); return

    # ======================= 4. merge =======================
    given = pd.read_csv(a.line_table, sep="\t")
    al = {pn: GS.Aligned(P[pn], given) for pn in panels}
    del given
    ctx = Ctx()
    ctx.cfg, ctx.P = cfg, P
    ctx.s_level5 = {pn: C.response_size(P[pn].R, P[pn].avail) for pn in panels}
    act = I.activity()
    ctx.activity = {pn: dict(zip(act[act.panel == pn].cmap_name, act[act.panel == pn].above_panel_median.astype(bool))) for pn in panels}
    pdef = I.panel_definitions()
    ctx.auc, prism_cov = {}, []
    pa = I.prism_auc(P["P1"], pdef)
    if pa is not None:
        ctx.auc["P1"] = pa[0]; prism_cov.append(pa[1])
    old = set(I.old_list().cmap_name)
    C.log(TAG, "merge: designs built; response sizes computed", t0)
    res, h1_rows, h2_parts = {}, [], []
    for pn in panels:
        for term in terms:
            d = al[pn].design(term, prim)
            if term == "strict" and cfg.decide("7c") == "C":
                r = GA.run_main_term(ctx, pn, term, d, 2, keep_perm=False)
                if r["per_cpd"] is not None:
                    for c in ("p", "q", "rank", "rule1_q", "rule2_top5pct"):
                        r["per_cpd"][c] = np.nan
                    r["panel_row"].update(panel_p=np.nan, note="7c C: estimates only, no p-value; null not calibrated")
            else:
                r = GA.run_main_term(ctx, pn, term, d, B, keep_perm=(pn == "P1"))
            res[(pn, term)] = r
            h1_rows.append(r["panel_row"])
            C.log(TAG, f"merge: {pn} {term} done", t0)
    # H2 rule 4 and findings (P1, Y and X); other panels and terms carry the rule columns without rule 4
    for (pn, term), r in res.items():
        if r["per_cpd"] is None:
            continue
        if pn == cfg.const("primary_panel") and term in ("factY", "factX"):
            r4 = GA.rule4(ctx, pn, r)
        else:
            r4 = pd.DataFrame(columns=["panel", "term", "compound", "rule4_met"])
        r["per_cpd"] = GA.apply_findings(ctx, pn, r["per_cpd"], r4, old)
        h2_parts.append(r["per_cpd"])
    h1 = pd.DataFrame(h1_rows)
    h2 = pd.concat(h2_parts, ignore_index=True)
    # gene level: findings and the top-ranked compound of each P1 Y and X term
    gl = []
    for term in ("factY", "factX"):
        r = res[("P1", term)]
        if r["per_cpd"] is None:
            continue
        t = r["per_cpd"]
        pick = set(t[t.finding | t.rule1_q.fillna(False)].compound) | set(t[t["rank"] == 1].compound)
        pick = [c for c in pick if c in P["P1"].comps]
        if pick:
            gl.append(GA.gene_level(ctx, "P1", r, sorted(pick), B))
    # ---------------- H3 (E3) and S12; H4 rows on the same scales
    l3_by_panel = {"P1": l3obj, "P1b": cfg.const("level3_tensors")["P1b"], "P2": cfg.const("level3_tensors")["P2"]}
    S_line = {}
    for pn in panels:
        S_line[(pn, "level5")] = (np.nanmean(ctx.s_level5[pn], axis=1), f"{pn}_resp")
        L3 = I.level3(l3_by_panel[pn], P[pn], order="D")
        if L3 is not None:
            S_line[(pn, "level3")] = (np.nanmean(C.response_size(L3[0], ~np.isnan(L3[0][:, :, 0])), axis=1), l3_by_panel[pn])
    other = {"level3": "level5", "level5": "level3"}
    h3_rows = []
    for pn in panels:
        for term in terms:
            d = res[(pn, term)].get("design") or al[pn].design(term, prim)
            for scale, role in ((scale_E3, "primary"), (other[scale_E3], "S12")):
                if (pn, scale) not in S_line:
                    h3_rows.append(dict(panel=pn, term=term, scale=scale, role=role, note="tensor not available in this run")); continue
                S, tn = S_line[(pn, scale)]
                row = GA.e3_row(ctx, pn, term, d, S, scale, tn, role)
                if scale == "level5" and role == "primary":
                    row["note"] = "level 5 (K4 failed or 7d C): biased toward zero, no equivalence claim" if cfg.decide("7d") != "C" else "level 5 (7d C)"
                key = {"factY": "Y_plus", "factX": "X_plus", "annotated": "g1_vs_g2"}.get(term)
                pf = pl_free.get(key, np.nan) if (pn == "P1" and role == "primary" and scale == scale_E3) else np.nan
                row["labelfree_P_inside_margin"] = pf
                row["equivalence"] = GA.equivalence(row, pf, cfg) if (pn == "P1" and role == "primary" and term != "strict") else "not applicable (secondary panel, S12 or strict)"
                if scale == "level5" and role == "primary" and cfg.decide("7d") != "C":
                    row["equivalence"] = "no equivalence claim (K4: E3 on level 5)"
                h3_rows.append(row)
    # E3 with the line's PRISM viability score as a covariate (draft section 6, mixed results): the score is the line's
    # mean over P1 compounds of the per-compound z-scored AUC; lines without any AUC leave this model
    if "P1" in ctx.auc and ("P1", scale_E3) in S_line:
        A = ctx.auc["P1"]
        with np.errstate(invalid="ignore", divide="ignore"):
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                Zs = (A - np.nanmean(A, axis=0)) / np.nanstd(A, axis=0, ddof=1)
                vscore = np.nanmean(np.where(np.isfinite(Zs), Zs, np.nan), axis=1)
        S, tn = S_line[("P1", scale_E3)]
        for term in terms:
            if term == "strict":
                continue
            d = res[("P1", term)].get("design") or al["P1"].design(term, prim)
            dv = dict(d, Xo=np.c_[vscore, d["Xo"]], k=d["k"] + 1, tcols=[c + 1 for c in d["tcols"]], mask=d["mask"] & np.isfinite(vscore))
            row = GA.e3_row(ctx, "P1", term, dv, S, scale_E3, tn, "viability-adjusted")
            prow = [r for r in h3_rows if r["panel"] == "P1" and r["term"] == term and r["role"] == "primary"]
            row["equivalence"] = "not applicable (viability-adjusted E3)"
            row["difference_disappears_after_viability_adjustment"] = bool(len(prow) and prow[0].get("difference_declared", False)
                                                                          and not row.get("difference_declared", False))
            h3_rows.append(row)
    h3 = pd.DataFrame(h3_rows)
    # ---------------- H4 attenuation
    h4 = []
    for pn in panels:
        dA = res[(pn, "annotated")].get("design") or al[pn].design("annotated", prim)
        att = GA.h4_attenuation(al[pn], dA)
        h4.append(dict(panel=pn, **att))
    h4 = pd.DataFrame(h4)
    C.log(TAG, "merge: H1, H2, H3, H4 done", t0)
    # ---------------- comparison exposures, NC5 and E5
    ctrl, e45 = [], []
    csumm = pd.DataFrame(calib_summ)
    for pn in panels:
        d = res[(pn, "factY")].get("design") or al[pn].design("factY", prim)
        tpX, _ = GS.tp53_columns(P[pn].cov, tcode)
        m = d["mask"]
        lin = d["lin"]
        tp1 = GS.tp53_columns(P[pn].cov, "7gB" if tcode == "7gB" else "cmp")[0]
        tp_row = GA.exposure_run(ctx, pn, "TP53 status (comparison exposure)", tp1, 0, [0], m, lin, prim["rule"], B)
        hi = GA.ploidy_split(P[pn], m, lin)
        pl_row = GA.exposure_run(ctx, pn, "ploidy within-lineage median split", np.c_[tpX, hi], tpX.shape[1], [tpX.shape[1]], m, lin, prim["rule"], B)
        q = P[pn].cov[f"qc_pass_frac_{pn}"].to_numpy(dtype=float)
        qlow = np.where(np.isfinite(q), (q < cfg.const("qc_low_threshold")).astype(float), np.nan)
        qc_row = GA.exposure_run(ctx, pn, "NC5 line QC status (qc_pass fraction < 0.2)", np.c_[tpX, qlow], tpX.shape[1], [tpX.shape[1]], m, lin, prim["rule"], B)
        ctrl += [tp_row, pl_row, qc_row]
        for term in terms:
            pr_ = res[(pn, term)]["panel_row"]
            cs = csumm[(csumm.panel == pn) & (csumm.term == term)]
            base = dict(panel=pn, term=term, panel_stat=pr_.get("panel_stat", np.nan), panel_p=pr_.get("panel_p", np.nan),
                        E5_excess=pr_.get("E5_excess", np.nan), n_p_lt_005=pr_.get("n_p_lt_005", np.nan), E4_threshold=pr_.get("E4_threshold", np.nan),
                        E4_met=pr_.get("E4_met", np.nan), n_q_lt_01=pr_.get("n_q_lt_01", np.nan))
            for c in (tp_row, pl_row):
                e45.append(dict(base, comparison=c["exposure"], comparison_panel_stat=c.get("panel_stat", np.nan),
                                comparison_panel_p=c.get("panel_p", np.nan), E5_comparison=c.get("E5_excess", np.nan)))
            e45.append(dict(base, comparison="NC2 random partitions (calibration pseudo-labels)",
                            comparison_panel_stat=np.nan, comparison_panel_p=np.nan, E5_comparison=np.nan,
                            random_partition_panel_stat_p05=cs.panel_stat_p05.iloc[0] if len(cs) and "panel_stat_p05" in cs else np.nan,
                            random_partition_panel_stat_p95=cs.panel_stat_p95.iloc[0] if len(cs) and "panel_stat_p95" in cs else np.nan))
            ctrl.append(dict(panel=pn, exposure=f"NC3 DMSO pseudo-compound under {term}", dmso_T=pr_.get("dmso_T", np.nan),
                             dmso_p=pr_.get("dmso_p", np.nan), dmso_rank_if_inserted=pr_.get("dmso_rank_if_inserted", np.nan),
                             reading_rule_7b=cfg.decide("7b")))
    ctrl = pd.DataFrame(ctrl); e45 = pd.DataFrame(e45)
    C.log(TAG, "merge: controls and comparison exposures done", t0)
    # ---------------- P2 agreement (7a)
    rep_c, rep_p = [], []
    for term in ("factY", "factX"):
        r1, r2 = res[("P1", term)], res[("P2", term)]
        p1p, p2p = r1["panel_row"].get("panel_p", np.nan), r2["panel_row"].get("panel_p", np.nan)
        seen = bool(p1p < cfg.const("H1_alpha_per_term") and p2p < cfg.const("replication_p"))
        rep_p.append(dict(term=term, P1_panel_p=p1p, P2_panel_p=p2p, P1_rejected=bool(p1p < cfg.const("H1_alpha_per_term")),
                          P2_p_lt_005=bool(p2p < cfg.const("replication_p")), label=GA.p2_labels(cfg, seen) if p1p < cfg.const("H1_alpha_per_term") else "P1 not rejected: no P2 label"))
        if r1["recs"] is not None and r2["recs"] is not None:
            ag = GA.p2_agreement(ctx, r1, r2, B)
            if len(ag):
                ag["label"] = [GA.p2_labels(cfg, x) if cfg.decide("7a") == "A" else
                               ("compound-level agreement reported, not called replication (7a C)" if cfg.decide("7a") == "C" else "description only (7a B)")
                               for x in ag.shape_and_E2_agree]
                rep_c.append(ag)
    rep_p = pd.DataFrame(rep_p); rep_c = pd.concat(rep_c, ignore_index=True) if rep_c else pd.DataFrame()
    # ---------------- leads (7l)
    lrows = []
    opt7l = cfg.decide("7l")
    for lead in cfg.const("leads"):
        for pn in lead["panels"]:
            d = res[(pn, lead["term"])].get("design") or al[pn].design(lead["term"], prim)
            if opt7l == "B" and (pn, "level3") in S_line:
                L3 = I.level3(l3_by_panel[pn], P[pn], order="D")
                s_scale, sc = C.response_size(L3[0], ~np.isnan(L3[0][:, :, 0])), "E3 scale (level 3)"
            else:
                s_scale, sc = ctx.s_level5[pn], "level 5"
            shape = (opt7l == "C" and lead["id"] == "L4")
            r = GA.lead_test(ctx, pn, lead, d, B, s_scale, use_shape=shape, h2_table=res[(pn, lead["term"])]["per_cpd"])
            r["scale"] = "level 5 (shape)" if shape else sc
            lrows.append(r)
    leads = pd.DataFrame(lrows)
    prim_leads = leads.panel == cfg.const("primary_panel")
    leads["p_holm"] = np.nan
    hp = prim_leads & leads["p"].notna() if "p" in leads else prim_leads & False
    if hp.any():
        leads.loc[hp, "p_holm"] = GA.holm(leads.loc[hp, "p"].values)
    C.log(TAG, "merge: P2 agreement and leads done", t0)
    # ---------------- secondary family
    reg = I.collectri()
    sec_rows, sec_cpd, xl_rows = [], [], []
    opt7n, opt7e = cfg.decide("7n"), cfg.decide("7e")
    fam_main = {"A": "all", "B": "single", "C": "single"}[opt7n]
    Wc, tfc, infc = GSC.collectri_weights(reg, P["P1"].genes, cfg.const("collectri_min_landmark_targets"), fam_main)
    Wx, tfx, infx = GSC.collectri_weights(reg, P["P1"].genes, cfg.const("collectri_min_landmark_targets"), "complex")
    lr = I.lopesramos()
    if opt7e == "A":
        pan = I.panda()
        tf_all = sorted(pan.TF.unique())
        Wp, tfp, _ = GSC.weights_from_edges(pan, P["P1"].genes, tf_all, "TF", "gene_symbol", "weight")
        focus = [t for t in lr.TF if t in tfp]
        not_testable = [t for t in lr.TF if t not in tfp]
    elif opt7e == "B":
        n_all = reg.groupby("tf").size()
        tf_all = sorted(t for t in n_all.index if "_" not in t)
        Wp, tfp, _ = GSC.weights_from_edges(reg[reg.tf.isin(tf_all)], P["P1"].genes, tf_all, "tf", "target", "sign")
        focus = [t for t in lr.TF if t in tfp]
        not_testable = [t for t in lr.TF if t not in tfp]
    xlinked = I.chromosomes().set_index("gene_symbol").loc[P["P1"].genes].is_X.values.astype(float)
    for pn in panels:
        for term in ("factY", "factX"):
            r = res[(pn, term)]
            if r["per_cpd"] is None:
                continue
            d = r["design"]
            h1set = set(r["per_cpd"][r["per_cpd"].in_ranking & (r["per_cpd"].complete_over_analysis_set if cfg.decide("7h") == "A" else True)].compound)
            h1_sub = (lambda df, s=h1set: df.compound.isin(s).values & df.in_ranking.values)
            gene_snr = r["panel_row"]["panel_stat"] / r["panel_row"]["null_median"] if r["panel_row"].get("null_median") else np.nan
            for fam_lab, W_, tfs_, perm in (("CollecTRI TF layer", Wc, tfc, True), ("CollecTRI complexes", Wx, tfx, opt7n == "C")):
                if opt7n == "A" and fam_lab == "CollecTRI complexes":
                    continue
                tl = GSC.tf_layer(P[pn], d, W_, B, h1_sub, with_perm=perm)
                if tl is None:
                    continue
                tdf = tl["df"]
                if perm:
                    tdf["q_bh_readout"] = np.nan
                    tdf.loc[tdf.in_ranking, "q_bh_readout"] = C.bh(tdf.loc[tdf.in_ranking, "p"].values)
                sec_rows.append(dict(panel=pn, term=term, readout=fam_lab + ("" if perm else " (descriptive, 7n B)"), n_tfs=len(tfs_),
                                     n_compounds=tl["n"], panel_stat=tl["panel_stat"], panel_p=tl["panel_p"],
                                     n_q_lt_01=int((tdf.get("q_bh_readout", pd.Series(dtype=float)) < cfg.const("secondary_bh_q")).sum()) if perm else np.nan,
                                     signal_to_null=tl["panel_stat"] / tl["null_median"] if perm and tl["null_median"] else np.nan,
                                     gene_signal_to_null=gene_snr, label="secondary"))
                keep = ["compound", "tested", "n_lines", "T_obs"] + (["p", "q_bh_readout", "rank"] if perm else [])
                sec_cpd.append(tdf[keep].assign(panel=pn, term=term, readout=fam_lab))
            if opt7e in ("A", "B") and focus:
                fz = GSC.focused(P[pn], d, Wp, tfp, focus, B, random_n, cfg.run("random_set_seed"), h1_sub)
                if fz is not None:
                    sec_rows.append(dict(panel=pn, term=term, readout=f"focused Lopes-Ramos TFs ({'PANDA-LCL' if opt7e == 'A' else 'CollecTRI'}, 7e {opt7e})",
                                         n_tfs=len(focus), n_compounds=fz["n"], panel_stat=fz["panel_stat"], panel_p=fz["panel_p"],
                                         competitive_p=fz["competitive_p"], n_random_sets=fz["n_random_sets"],
                                         not_testable=";".join(not_testable), label="secondary, one-sided"))
                    fdf = fz["df"]
                    sec_cpd.append(fdf[["compound", "tested", "n_lines", "T_obs", "p", "rank"]].assign(panel=pn, term=term, readout="focused Lopes-Ramos TFs"))
            if term == "factX":
                xdf, xp = GSC.xlinked(P[pn], d, xlinked, B, random_n, cfg.run("random_set_seed") + 1, h1set)
                xdf.insert(0, "panel", pn)
                xl_rows.append(xdf)
                sec_rows.append(dict(panel=pn, term=term, readout="X-linked landmark score (centred)", **xp, label="secondary, two-sided"))
    sec = pd.DataFrame(sec_rows)
    sec_cpd = pd.concat(sec_cpd, ignore_index=True) if sec_cpd else pd.DataFrame()
    xl = pd.concat(xl_rows, ignore_index=True) if xl_rows else pd.DataFrame()
    C.log(TAG, "merge: secondary family done", t0)
    # ---------------- sensitivity analyses
    sens, sens_cpd = [], []
    pm = res[("P1", "factY")]["design"]["mask"]
    for sid, desc, st in GSN.variants(cfg, prim, pm):
        pns = ["P1", "P2"] if sid == "S2" else ["P1"]
        for pn in pns:
            for term in (["factY"] if st.get("y_only") else ["factY", "factX"]):
                st2 = dict(st)
                if st2.get("same_lines") is not None and pn != "P1":
                    st2["same_lines"] = res[(pn, "factY")]["design"]["mask"]
                d = al[pn].design(term, st2)
                pdz = res[(pn, term)].get("design") or al[pn].design(term, prim)
                if sid == "S2" and (d["mask"] == pdz["mask"]).all():
                    sens.append(dict(sensitivity=sid, description=desc, panel=pn, term=term, n_lines=d["n_lines"], note="no line dropped: identical to the primary analysis"))
                    continue
                r = GA.run_main_term(ctx, pn, term, d, B, label=sid)
                pr_ = r["panel_row"]
                sens.append(dict(sensitivity=sid, description=desc, panel=pn, term=term, **{k: pr_.get(k) for k in
                            ("n_lines", "n_g1", "n_g0", "n_lineages", "contrast_ok", "n_tested", "panel_stat", "panel_p", "n_p_lt_005", "n_q_lt_01")}, note=""))
                if r["per_cpd"] is not None:
                    sens_cpd.append(r["per_cpd"][["panel", "term", "compound", "tested", "T_obs", "p", "q", "rank", "E2_est"]].assign(sensitivity=sid))
        C.log(TAG, f"merge: sensitivity {sid} done", t0)
    ps = I.persig_level5()
    si = I.siginfo(["sig_id", "cell_iname", "cmap_name", "pert_idose", "nsample"]) if ps is not None else None
    if ps is not None and si is not None:
        Q, s11info = GSN.build_s11(P["P1"], ps, si)
        ctxQ = Ctx(); ctxQ.cfg, ctxQ.P = cfg, dict(P, P1=Q)
        ctxQ.s_level5 = dict(ctx.s_level5, P1=C.response_size(Q.R, Q.avail)); ctxQ.activity = ctx.activity
        for term in ("factY", "factX"):
            r = GA.run_main_term(ctxQ, "P1", term, res[("P1", term)]["design"], B, label="S11")
            pr_ = r["panel_row"]
            sens.append(dict(sensitivity="S11", description="signatures with nsample >= 2", panel="P1", term=term, **{k: pr_.get(k) for k in
                        ("n_lines", "n_g1", "n_g0", "n_lineages", "contrast_ok", "n_tested", "panel_stat", "panel_p", "n_p_lt_005", "n_q_lt_01")},
                        note="; ".join(f"{k}={v}" for k, v in s11info.items())))
            if r["per_cpd"] is not None:
                sens_cpd.append(r["per_cpd"][["panel", "term", "compound", "tested", "T_obs", "p", "q", "rank", "E2_est"]].assign(sensitivity="S11"))
    else:
        sens.append(dict(sensitivity="S11", description="signatures with nsample >= 2", panel="P1", note="per-signature level 5 or siginfo not available in this run"))
    for term in ("factY", "factX"):
        pr_ = res[("P1", term)]["panel_row"]
        sens.append(dict(sensitivity="S4", description="P1b (0.1 uM)", panel="P1b", term=term, **{k: res[("P1b", term)]["panel_row"].get(k) for k in ("n_lines", "n_tested", "panel_stat", "panel_p", "n_p_lt_005", "n_q_lt_01")}, note="rows of G_H1.tsv"))
        sens.append(dict(sensitivity="7h", description=f"H1 over {pr_.get('sens_compound_set')}", panel="P1", term=term, n_tested=pr_.get("sens_n"),
                         panel_stat=pr_.get("sens_panel_stat"), panel_p=pr_.get("sens_panel_p"), note="other compound set of 7h"))
        sens.append(dict(sensitivity="activity", description="H1 over compounds above the panel median activity", panel="P1", term=term,
                         n_tested=pr_.get("activity_subset_n"), panel_stat=pr_.get("activity_subset_panel_stat"), panel_p=pr_.get("activity_subset_panel_p"), note=""))
        s12 = h3[(h3.panel == "P1") & (h3.term == term) & (h3.role == "S12")]
        if len(s12):
            sens.append(dict(sensitivity="S12", description=f"E3 on {s12.scale.iloc[0]}", panel="P1", term=term, note="rows of G_H3.tsv with role S12"))
    sens = pd.DataFrame(sens); sens_cpd = pd.concat(sens_cpd, ignore_index=True) if sens_cpd else pd.DataFrame()
    # ---------------- viability step
    vt, vs = [], []
    for term in ("factY", "factX"):
        if res[("P1", term)]["recs"] is None:
            continue
        tab, summ = GV.run(ctx, res[("P1", term)], relabel_n)
        vt.append(tab); vs.append(summ)
    vt = pd.concat(vt, ignore_index=True) if vt else pd.DataFrame(); vs = pd.DataFrame(vs)
    C.log(TAG, "merge: sensitivity and viability done", t0)
    # ---------------- decisions, figures, outputs
    dec = GD.build(cfg, gt, h1, h2, h3, e45, rep_p, leads, sec, [])
    w(out, "G_H1.tsv", h1); w(out, "G_H2_per_compound.tsv", h2)
    pd.concat(gl, ignore_index=True).to_csv(os.path.join(out, "G_H2_gene_level.tsv.gz"), sep="\t", index=False) if gl else \
        pd.DataFrame(columns=["panel", "term", "compound", "gene"]).to_csv(os.path.join(out, "G_H2_gene_level.tsv.gz"), sep="\t", index=False)
    w(out, "G_H3.tsv", h3); w(out, "G_H4.tsv", h4); w(out, "G_E4_E5.tsv", e45); w(out, "G_controls.tsv", ctrl)
    w(out, "G_replication_P2.tsv", rep_p); w(out, "G_replication_P2_compounds.tsv", rep_c); w(out, "G_leads.tsv", leads)
    w(out, "G_secondary_tf.tsv", sec); w(out, "G_secondary_tf_per_compound.tsv", sec_cpd); w(out, "G_secondary_xlinked.tsv", xl)
    w(out, "G_sensitivity.tsv", sens); w(out, "G_sensitivity_per_compound.tsv", sens_cpd)
    w(out, "G_viability.tsv", vs); w(out, "G_viability_per_compound.tsv", vt)
    w(out, "G_prism_coverage.tsv", pd.concat(prism_cov, ignore_index=True) if prism_cov else pd.DataFrame(columns=["panel", "compound"]))
    w(out, "G_decisions.tsv", dec)
    fo = os.path.join(out, "figures")
    perm_by = {t: (res[("P1", t)]["permT"], np.flatnonzero(res[("P1", t)]["per_cpd"].in_ranking.values))
               for t in ("factY", "factX") if res[("P1", t)]["permT"] is not None}
    if perm_by:
        GF.qq(perm_by, h2[h2.panel == "P1"], fo)
        h1cols = {}
        for t in perm_by:
            pc = res[("P1", t)]["per_cpd"]
            sel = pc.in_ranking.values & (pc.complete_over_analysis_set.values if cfg.decide("7h") == "A" else True)
            h1cols[t] = np.nanmean(res[("P1", t)]["permT"][:, np.flatnonzero(sel)], axis=1)
        GF.panel_null(h1cols, h1, fo)
    GF.e3(h3, cfg.const("E3_margin_pct"), fo)
    if len(rep_c):
        GF.p2_shape(rep_c, fo)
    if len(vt) and "E2_auc_adjusted" in vt:
        GF.viability(vt, fo)
    C.log(TAG, "merge: all tables and figures written", t0)
    finish(out, cfg, manifest, I, t0)


def finish(out, cfg, manifest, I, t0):
    w(out, "G_hash_check.tsv", I.hash_table())
    w(out, "G_decide_options.tsv", pd.DataFrame(cfg.decide_table()))
    sums = os.path.join(ROOT, cfg.raw["paths"]["code_sums"])
    code = []
    for f in sorted(os.listdir(HERE)):
        if f.endswith((".py", ".sh")):
            code.append(dict(item=f"code_sha256 src/G/{f}", value=guard.sha256_file(os.path.join(HERE, f))))
    if os.path.exists(sums):
        code.append(dict(item="code_sha256 src/G/SHA256SUMS.txt", value=guard.sha256_file(sums)))
    manifest = manifest + code + [dict(item="decide_ids_read", value=";".join(cfg.read_ids)), dict(item="runtime_s", value=f"{time.time()-t0:.0f}")]
    w(out, "G_run_manifest.tsv", pd.DataFrame(manifest))
    C.log(TAG, f"done in {time.time()-t0:.0f}s; outputs in {out}")


if __name__ == "__main__":
    main()
