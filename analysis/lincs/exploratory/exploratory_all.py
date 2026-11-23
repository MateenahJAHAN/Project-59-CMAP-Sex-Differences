"""Exploratory runs on all panels and all idarubicin conditions (extends X1 and X2 at the coordinator's request).

EXPLORATORY. Frozen functions in src/G, unchanged; only response tensors differ. Every result is written, whatever it shows.
X1 on P1b and P2: the H1 test on vehicle-referenced responses (level 3: each well minus the same plate's DMSO mean;
  level 5 minus the line's DMSO pseudo-compound response; P1b uses P1's DMSO response, same lines, project and time).
  For each panel: the calibration code is first validated on the original level-5 tensor against results/G/G_calibration.tsv
  (must be equal), then the pseudo-label calibration runs on each new tensor, then the real-label test.
X2 on all idarubicin conditions with a testable design: the L4 statistic (ganalyses.lead_test, E2, Y term, two-sided) and
  the shape test, at "0.08 uM" 6 h (CPC; P2 design; read by range reads), at "0.12 uM" 24 h (P1b tensor; the condition of
  an earlier look) and, for the "12 uM" 24 h condition, the count of eligible lines in every P1 line with data.
Usage (repository root): <task G python> -I src/G_run/exploratory_all.py
"""
import hashlib
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
ROOT = os.getcwd()
from study import open_study
import numpy as np
import pandas as pd
import fsspec
import h5py
import gcore as C
import gsets as GS
import ggates as GG
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from calibration import run_calibration
GG.run_calibration = run_calibration
import ganalyses as GA

URL = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/level5_beta_trt_cp_n720216x12328.gctx"
OUT = os.path.join(ROOT, "results", "G_exploratory")
DATA = os.path.join(ROOT, "data", "G", "exploratory")
WORKERS = int(os.environ.get("P59_CAL_WORKERS", 4))
CAL_SOURCE = {"P1b": "D2_R1_F41", "P2": "D_P2_factorial"}


class Ctx:
    pass


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def log(msg):
    line = f"{pd.Timestamp.now(tz='UTC').strftime('%Y-%m-%dT%H:%M:%SZ')} [all] {msg}"
    print(line, flush=True)
    with open(os.path.join(ROOT, "logs", "G_exploratory_runlog.txt"), "a") as f:
        f.write(line + "\n")


def ctx_for(cfg, I, pn, Q):
    act = I.activity()
    ctx = Ctx()
    ctx.cfg, ctx.P = cfg, {pn: Q}
    ctx.s_level5 = {pn: C.response_size(Q.R, Q.avail)}
    ctx.activity = {pn: dict(zip(act[act.panel == pn].cmap_name, act[act.panel == pn].above_panel_median.astype(bool)))}
    return ctx


def x1_panel(cfg, I, cov, pn, given, conf, P1dmso_by_line):
    P = I.panel(pn, cov, cfg.run("B_full"))
    prim, tcode = GS.primary_settings(cfg), GS.tp53_code(cfg)
    key_seed = cfg.run("key_seed_base") + cfg.const("panel_numbers")[pn]
    lp, cp = I.file(f"taskD/tensors/{pn}_lines.tsv"), I.file(f"taskD/tensors/{pn}_compounds.tsv")
    l3obj = cfg.const("level3_tensors")[pn]
    L3 = I.level3(l3obj, P, order="D")[0]
    if P.dmso is not None:
        dm, dm_src = P.dmso, f"{pn}_DMSO_resp.npy"
    else:
        dm, dm_src = np.stack([P1dmso_by_line[l] for l in P.lines]), "P1_DMSO_resp.npy (same lines, project and time)"
    tens = {"level5_original": I.d_tensor(f"{pn}_resp.npy")}
    for name, arr in ((f"X1a_level3_plateDMSO", L3), (f"X1b_level5_minus_lineDMSO", P.R - dm[:, None, :])):
        p = os.path.join(DATA, f"{name}_{pn}_Dorder.npy")
        np.save(p, arr.astype(np.float64))
        tens[name] = p
    trows = [dict(panel=pn, tensor=k, path=os.path.relpath(v, ROOT), sha256=sha(v), note=(dm_src if k.startswith("X1b") else (l3obj if k.startswith("X1a") else ""))) for k, v in tens.items()]
    src = [s for s in cfg.raw["pseudo_sources"]["list"] if s["name"] == CAL_SOURCE[pn]][0]
    fl, flin, cd = GG.read_pseudo(I.file(src["file"]))
    n_use = min(cfg.run("calibration_n_pseudo")[pn], cd.shape[1])
    cal_rows, cal_summ, panels = [], [], {}
    for name, path in tens.items():
        pa = (pn, path, lp, cp, P.genes, cov, key_seed, cfg.run("B_full"), None)
        Q = C.Panel(*pa)
        assert list(Q.lines) == list(P.lines) and list(Q.comps) == list(P.comps)
        panels[name] = Q
        tpX, _ = GS.tp53_columns(Q.cov, tcode)
        full, elig = GG.codes_full(Q, fl, cd)
        for term in ("factY", "factX"):
            rows = GG.run_calibration(pa, Q, tpX, full, elig, "factorial", term, cfg.run("calibration_B"), prim["rule"], range(n_use), WORKERS)
            sm = GG.summarize(rows, cfg.const("calibration_gate"))
            note = ""
            if name == "level5_original":
                a = pd.DataFrame(rows).set_index("pseudo").sort_index()
                b = conf[(conf.panel == pn) & (conf.term == term)].set_index("pseudo").sort_index().loc[a.index]
                eq = all(np.allclose(a[c].astype(float), b[c].astype(float), rtol=1e-9, atol=0, equal_nan=True)
                         for c in ("n_tested", "n_p_lt_005", "n_q_lt_01", "min_p", "panel_stat", "panel_p"))
                note = "validation against results/G/G_calibration.tsv: " + ("all rows equal" if eq else "ROWS DIFFER")
                assert eq, f"{pn} {term}: calibration code does not reproduce the confirmatory rows"
            for r in rows:
                r.update(panel=pn, tensor=name, label="EXPLORATORY X1")
            cal_rows += rows
            cal_summ.append(dict(panel=pn, tensor=name, term=term, label="EXPLORATORY X1", source=CAL_SOURCE[pn], **sm, note=note))
            log(f"X1 {pn} calibration {name} {term}: mean share {sm.get('mean_share_p_lt_005', float('nan')):.4f}, p95 {sm.get('p95_n_p_lt_005')}, "
                f"in range {sm['gate_pass']} {note}")
    cs = pd.DataFrame(cal_summ)
    h1_rows, cpd = [], []
    for name in ("X1a_level3_plateDMSO", "X1b_level5_minus_lineDMSO"):
        Q = panels[name]
        al = GS.Aligned(Q, given)
        ctx = ctx_for(cfg, I, pn, Q)
        for term in ("factY", "factX"):
            r = GA.run_main_term(ctx, pn, term, al.design(term, prim), cfg.run("B_full"), label=f"EXPLORATORY X1 {name}")
            pr = dict(r["panel_row"]); pr.pop("null", None)
            c = cs[(cs.tensor == name) & (cs.term == term)].iloc[0]
            pr.update(tensor=name, label="EXPLORATORY X1", calib_mean_share=c.mean_share_p_lt_005, calib_in_range=bool(c.gate_pass),
                      calib_p95_count=c.p95_n_p_lt_005)
            if not pr["calib_in_range"]:
                pr["reading"] = "uncalibrated: not read"
            elif pr.get("n_tested", 0) and pr["n_p_lt_005"] <= pr["calib_p95_count"] and pr["panel_p"] >= 0.05:
                pr["reading"] = "no excess after vehicle referencing"
            else:
                pr["reading"] = "excess or panel p < 0.05 after vehicle referencing"
            h1_rows.append(pr)
            if r["per_cpd"] is not None:
                cpd.append(r["per_cpd"].drop(columns=[x for x in r["per_cpd"].columns if x.startswith("_")], errors="ignore").assign(tensor=name))
            log(f"X1 {pn} real labels {name} {term}: panel p {pr.get('panel_p', float('nan')):.4f}, n p<0.05 {pr.get('n_p_lt_005')}, p95 {pr['calib_p95_count']}")
    return pd.DataFrame(cal_rows), cs, pd.DataFrame(h1_rows), (pd.concat(cpd, ignore_index=True) if cpd else pd.DataFrame()), pd.DataFrame(trows)


def read_signatures(sig_ids, lm):
    fs = fsspec.filesystem("https", block_size=4 * 1024 * 1024)
    with fs.open(URL, "rb") as fo, h5py.File(fo, "r") as h:
        col = [x.decode() for x in h["0/META/COL/id"][:]]
        row = [x.decode() for x in h["0/META/ROW/id"][:]]
        ds = h["0/DATA/0/matrix"]
        assert ds.id.get_create_plist().get_layout() == h5py.h5d.CONTIGUOUS
        offset, dtype = ds.id.get_offset(), ds.dtype
        pos = {s: i for i, s in enumerate(col)}
        order = sorted(sig_ids, key=lambda s: pos[s])
        M = np.stack([ds[pos[s], :] for s in order]).astype(np.float32)
    assert [row[i] for i in lm.gctx_row_index] == [str(g) for g in lm.gene_id]
    nb = len(row) * np.dtype(dtype).itemsize
    for s in order[:5]:
        st = offset + pos[s] * nb
        raw = np.frombuffer(fs.cat_file(URL, start=st, end=st + nb), dtype=dtype)
        assert np.array_equal(raw, M[order.index(s)]), "byte-range read differs"
    return order, M[:, lm.gctx_row_index.to_numpy()], nb


def lead_on(cfg, I, cov, pn, tensor_path, given, label):
    prim = GS.primary_settings(cfg)
    key_seed = cfg.run("key_seed_base") + cfg.const("panel_numbers")[pn]
    P = I.panel(pn, cov, cfg.run("B_full"))
    pa = (pn, tensor_path, I.file(f"taskD/tensors/{pn}_lines.tsv"), I.file(f"taskD/tensors/{pn}_compounds.tsv"), P.genes, cov, key_seed,
          cfg.run("B_full"), None)
    Q = C.Panel(*pa)
    j = list(Q.comps).index("idarubicin")
    s = C.response_size(Q.R, Q.avail)
    d = GS.Aligned(Q, given).design("factY", prim)
    rows_all = np.flatnonzero(d["mask"])
    rows = rows_all[Q.avail[rows_all, j]]
    kept, dd, why = C.compound_design(d["lin"], d["Xo"], d["k"], rows, int(d["rule"]["min_group"]), int(d["rule"]["min_relab"]))
    tv = d["Xo"][kept, d["k"]] if len(kept) else np.array([])
    out = dict(label=label, panel_design=pn, n_lines_with_data=int(Q.avail[:, j].sum()), n_set_lines_with_data=len(rows),
               n_lines_after_lineage_rule=len(kept), n_Y1=int(np.nansum(tv == 1)), n_Y0=int(np.nansum(tv == 0)), testable=dd is not None,
               reason="" if dd is not None else why)
    if dd is not None:
        ctx = Ctx(); ctx.cfg, ctx.P, ctx.s_level5 = cfg, {pn: Q}, {pn: s}
        lr = GA.lead_test(ctx, pn, dict(id="X2", compound="idarubicin", term="factY", alternative="two-sided", panels=[pn]), d, cfg.run("B_full"), s)
        out.update({f"lead_{k}": v for k, v in lr.items() if k not in ("lead", "compound", "panel", "term", "alternative")})
        df, _, _ = C.run_term(Q, d["Xo"], d["k"], d["mask"], cfg.run("B_full"), comps=[j], lin=d["lin"], rule=d["rule"])
        rr = df[df.compound == "idarubicin"].iloc[0]
        out.update(shape_tested=bool(rr.tested), shape_T_obs=float(rr.T_obs) if "T_obs" in df else np.nan, shape_p=float(rr.p) if "p" in df else np.nan)
    return out


def main():
    os.makedirs(OUT, exist_ok=True); os.makedirs(DATA, exist_ok=True)
    cfg, I = open_study(ROOT)
    cov = I.covariates()
    given = pd.read_csv(os.path.join(ROOT, "data/G/inputs/taskD/D_line_table.tsv"), sep="\t")
    conf = pd.read_csv(os.path.join(ROOT, "results/G/G_calibration.tsv"), sep="\t", float_precision="round_trip")
    P1 = I.panel("P1", cov, cfg.run("B_full"))
    P1dmso_by_line = {l: P1.dmso[i] for i, l in enumerate(P1.lines)}

    # ---- X1 on P1b and P2
    parts = [x1_panel(cfg, I, cov, pn, given, conf, P1dmso_by_line) for pn in ("P1b", "P2")]
    for k, name in enumerate(("X1all_calibration_rows", "X1all_calibration_summary", "X1all_H1", "X1all_per_compound", "X1all_tensors")):
        pd.concat([p[k] for p in parts], ignore_index=True).to_csv(os.path.join(OUT, f"{name}.tsv"), sep="\t", index=False)
    log("X1 on P1b and P2 written")

    # ---- X2 on all idarubicin conditions
    lm = I.landmarks()
    si = pd.read_csv(os.path.join(ROOT, "data/G/meta/siginfo_beta.txt"), sep="\t", low_memory=False,
                     usecols=["sig_id", "cmap_name", "pert_idose", "pert_itime", "pert_type", "project_code", "cell_iname"])
    res = []
    # (a) 0.08 uM, 6 h, CPC: P2 design, new range reads
    sel = si[(si.cmap_name == "idarubicin") & (si.pert_idose == "0.08 uM") & (si.pert_itime == "6 h") & (si.pert_type == "trt_cp") & (si.project_code == "CPC")]
    P2 = I.panel("P2", cov, cfg.run("B_full"))
    assert len(sel) > 0 and sel.cell_iname.isin(set(P2.lines)).all()
    t0 = time.time()
    order, Z, nb = read_signatures(list(sel.sig_id), lm)
    zp = os.path.join(DATA, "X2_idarubicin_0.08uM_6h_signatures_landmarks.npy"); np.save(zp, Z)
    with open(os.path.join(ROOT, "provenance/downloads.tsv"), "a") as f:
        f.write("\t".join(["G-exploratory", URL, os.path.relpath(zp, ROOT), str(len(order) * nb), sha(zp), pd.Timestamp.now(tz="UTC").strftime("%Y-%m-%dT%H:%M:%SZ"),
                           f"range reads of the {len(order)} idarubicin 0.08 uM 6 h CPC signatures; byte-range check 5 of 5 equal; sha256 of the saved landmark array"]) + "\n")
    log(f"X2 read {len(order)} signatures (0.08 uM, 6 h) in {time.time() - t0:.0f} s")
    if "idarubicin" in list(P2.comps):
        R = np.full(P2.R.shape, np.nan); j = list(P2.comps).index("idarubicin")
        lines_of = sel.set_index("sig_id").loc[order].cell_iname.values
        for l in sorted(set(lines_of)):
            R[list(P2.lines).index(l), j, :] = Z[lines_of == l].astype(np.float64).mean(axis=0)
        tp = os.path.join(DATA, "X2_idarubicin_0.08uM_6h_P2shape_Dorder.npy"); np.save(tp, R)
        res.append(dict(condition="0.08 uM, 6 h (CPC)", new_data=True, **lead_on(cfg, I, cov, "P2", tp, given, "EXPLORATORY X2")))
    else:
        # idarubicin is not a P2 compound: place it in a P2-shaped tensor by replacing a column is not allowed; build a one-compound panel instead
        comps = list(P2.comps)
        lp, cpth = I.file("taskD/tensors/P2_lines.tsv"), os.path.join(DATA, "X2_P2_idarubicin_compounds.tsv")
        pd.DataFrame({"cmap_name": ["idarubicin"]}).to_csv(cpth, sep="\t", index=False)
        R = np.full((len(P2.lines), 1, 978), np.nan)
        lines_of = sel.set_index("sig_id").loc[order].cell_iname.values
        for l in sorted(set(lines_of)):
            R[list(P2.lines).index(l), 0, :] = Z[lines_of == l].astype(np.float64).mean(axis=0)
        tp = os.path.join(DATA, "X2_idarubicin_0.08uM_6h_P2lines_Dorder.npy"); np.save(tp, R)
        prim = GS.primary_settings(cfg)
        pa = ("P2", tp, lp, cpth, P2.genes, cov, cfg.run("key_seed_base") + cfg.const("panel_numbers")["P2"], cfg.run("B_full"), None)
        Q = C.Panel(*pa)
        s = C.response_size(Q.R, Q.avail)
        d = GS.Aligned(Q, given).design("factY", prim)
        rows_all = np.flatnonzero(d["mask"]); rows = rows_all[Q.avail[rows_all, 0]]
        kept, dd, why = C.compound_design(d["lin"], d["Xo"], d["k"], rows, int(d["rule"]["min_group"]), int(d["rule"]["min_relab"]))
        tv = d["Xo"][kept, d["k"]] if len(kept) else np.array([])
        out = dict(condition="0.08 uM, 6 h (CPC)", new_data=True, label="EXPLORATORY X2", panel_design="P2", n_lines_with_data=int(Q.avail[:, 0].sum()),
                   n_set_lines_with_data=len(rows), n_lines_after_lineage_rule=len(kept), n_Y1=int(np.nansum(tv == 1)), n_Y0=int(np.nansum(tv == 0)),
                   testable=dd is not None, reason="" if dd is not None else why)
        if dd is not None:
            ctx = Ctx(); ctx.cfg, ctx.P, ctx.s_level5 = cfg, {"P2": Q}, {"P2": s}
            lr = GA.lead_test(ctx, "P2", dict(id="X2", compound="idarubicin", term="factY", alternative="two-sided", panels=["P2"]), d, cfg.run("B_full"), s)
            out.update({f"lead_{k}": v for k, v in lr.items() if k not in ("lead", "compound", "panel", "term", "alternative")})
            df, _, _ = C.run_term(Q, d["Xo"], d["k"], d["mask"], cfg.run("B_full"), comps=[0], lin=d["lin"], rule=d["rule"])
            rr = df.iloc[0]
            out.update(shape_tested=bool(rr.tested), shape_T_obs=float(rr.T_obs) if "T_obs" in df else np.nan, shape_p=float(rr.p) if "p" in df else np.nan)
        res.append(out)
    # (b) 0.12 uM, 24 h: P1b tensor as it is (condition of an earlier look)
    res.append(dict(condition="0.12 uM, 24 h (ASG; earlier look)", new_data=False, **lead_on(cfg, I, cov, "P1b", I.d_tensor("P1b_resp.npy"), given, "EXPLORATORY X2")))
    # (c) 12 uM, 24 h: the analysis set rebuilt for the P1 lines with data at this label (frozen set function, same eligibility)
    prim = GS.primary_settings(cfg)
    tp12 = os.path.join(DATA, "X2_idarubicin_12uM_P1shape_Dorder.npy")
    pa = ("P1", tp12, I.file("taskD/tensors/P1_lines.tsv"), I.file("taskD/tensors/P1_compounds.tsv"), P1.genes, cov,
          cfg.run("key_seed_base") + cfg.const("panel_numbers")["P1"], cfg.run("B_full"), None)
    Q12 = C.Panel(*pa)
    j = list(Q12.comps).index("idarubicin")
    d = GS.Aligned(Q12, given).design("factY", prim)
    av = Q12.avail[:, j]
    mask12 = C.analysis_set(d["lin"], d["Xo"], d["tcols"], d["base"] & av)
    d12 = dict(d, mask=mask12)
    rows = np.flatnonzero(mask12 & av)
    kept, dd, why = C.compound_design(d12["lin"], d12["Xo"], d12["k"], rows, int(d12["rule"]["min_group"]), int(d12["rule"]["min_relab"]))
    tv = d12["Xo"][kept, d12["k"]] if len(kept) else np.array([])
    out = dict(condition="12 uM, 24 h (ASG); analysis set rebuilt for the lines with data", new_data=True, label="EXPLORATORY X2",
               panel_design="P1 eligibility, set rebuilt", n_lines_with_data=int(av.sum()), n_eligible_with_data=int((d["base"] & av).sum()),
               n_set_lines_with_data=len(rows), n_lines_after_lineage_rule=len(kept), n_Y1=int(np.nansum(tv == 1)), n_Y0=int(np.nansum(tv == 0)),
               testable=dd is not None, reason="" if dd is not None else why)
    if dd is not None:
        s12 = C.response_size(Q12.R, Q12.avail)
        ctx = Ctx(); ctx.cfg, ctx.P, ctx.s_level5 = cfg, {"P1": Q12}, {"P1": s12}
        lr = GA.lead_test(ctx, "P1", dict(id="X2", compound="idarubicin", term="factY", alternative="two-sided", panels=["P1"]), d12, cfg.run("B_full"), s12)
        out.update({f"lead_{k}": v for k, v in lr.items() if k not in ("lead", "compound", "panel", "term", "alternative")})
        df, _, _ = C.run_term(Q12, d12["Xo"], d12["k"], d12["mask"], cfg.run("B_full"), comps=[j], lin=d12["lin"], rule=d12["rule"])
        rr = df[df.compound == "idarubicin"].iloc[0]
        out.update(shape_tested=bool(rr.tested), shape_T_obs=float(rr.T_obs) if "T_obs" in df else np.nan, shape_p=float(rr.p) if "p" in df else np.nan)
    res.append(out)
    pd.DataFrame(res).to_csv(os.path.join(OUT, "X2all_idarubicin_conditions.tsv"), sep="\t", index=False)
    log("X2 on all idarubicin conditions written")


if __name__ == "__main__":
    main()
