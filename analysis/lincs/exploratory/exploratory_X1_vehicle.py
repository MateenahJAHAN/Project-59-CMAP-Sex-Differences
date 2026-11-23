"""Exploratory X1 (analysis_plan/G_exploratory_X1_vehicle.md): the confirmatory H1 test on vehicle-referenced P1 responses.

EXPLORATORY. Uses the frozen functions in src/G unchanged; only the response tensor differs.
  X1a  task E's level-3 tensor R3_P1_d111 (each well minus the same plate's DMSO mean), in task D's order
  X1b  task D's level-5 P1 tensor minus each line's level-5 DMSO pseudo-compound response
Order: (0) validation: the calibration code below, run on the original level-5 tensor, must reproduce the confirmatory
calibration rows exactly; (1) blind pseudo-label calibration on X1a and X1b; (2) the real-label test.
Usage (repository root): <task G python> -I src/G_run/exploratory_X1_vehicle.py
"""
import hashlib
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
import ggates as GG
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from calibration import run_calibration
GG.run_calibration = run_calibration
import ganalyses as GA

OUT = os.path.join(ROOT, "results", "G_exploratory")
DATA = os.path.join(ROOT, "data", "G", "exploratory")
LABEL = "EXPLORATORY X1"
WORKERS = int(os.environ.get("P59_CAL_WORKERS", 4))


class Ctx:
    pass


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def log(msg):
    line = f"{pd.Timestamp.now(tz='UTC').strftime('%Y-%m-%dT%H:%M:%SZ')} [X1] {msg}"
    print(line, flush=True)
    with open(os.path.join(ROOT, "logs", "G_exploratory_runlog.txt"), "a") as f:
        f.write(line + "\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(DATA, exist_ok=True)
    cfg, I = open_study(ROOT)
    cov = I.covariates()
    P1 = I.panel("P1", cov, cfg.run("B_full"))
    prim = GS.primary_settings(cfg)
    tcode = GS.tp53_code(cfg)
    key_seed = cfg.run("key_seed_base") + cfg.const("panel_numbers")["P1"]
    lines_p, comps_p = I.file("taskD/tensors/P1_lines.tsv"), I.file("taskD/tensors/P1_compounds.tsv")

    # ---- tensors, written once in task D's order
    l3obj = cfg.const("level3_tensors")["P1"]["7dA"]
    L3 = I.level3(l3obj, P1, order="D")[0]
    assert P1.dmso is not None and P1.dmso.shape == (P1.R.shape[0], P1.R.shape[2]), None if P1.dmso is None else P1.dmso.shape
    R5adj = P1.R - P1.dmso[:, None, :]
    tens = {"level5_original": I.d_tensor("P1_resp.npy")}
    for name, arr in (("X1a_level3_plateDMSO", L3), ("X1b_level5_minus_lineDMSO", R5adj)):
        p = os.path.join(DATA, f"{name}_P1_Dorder.npy")
        np.save(p, arr.astype(np.float64))
        tens[name] = p
    pd.DataFrame([dict(tensor=k, path=os.path.relpath(v, ROOT), sha256=sha(v), bytes=os.path.getsize(v)) for k, v in tens.items()]).to_csv(
        os.path.join(OUT, "X1_tensors.tsv"), sep="\t", index=False)
    log(f"tensors written: {list(tens)}; level-3 object {l3obj}")

    # ---- pseudo-labels: the confirmatory run's source for P1 factorial
    src = [s for s in cfg.raw["pseudo_sources"]["list"] if s["name"] == "D2_R1_F41"][0]
    fl, flin, cd = GG.read_pseudo(I.file(src["file"]))
    n_use = min(cfg.run("calibration_n_pseudo")["P1"], cd.shape[1])
    calB = cfg.run("calibration_B")
    conf = pd.read_csv(os.path.join(ROOT, "results/G/G_calibration.tsv"), sep="\t", float_precision="round_trip")
    conf = conf[(conf.panel == "P1") & conf.term.isin(["factY", "factX"])]

    cal_rows, cal_summ = [], []
    panels = {}
    for name, path in tens.items():
        panel_args = ("P1", path, lines_p, comps_p, P1.genes, cov, key_seed, cfg.run("B_full"), None)
        Q = C.Panel(*panel_args)
        assert list(Q.lines) == list(P1.lines) and list(Q.comps) == list(P1.comps) and list(Q.genes) == list(P1.genes)
        panels[name] = (Q, panel_args)
        tpX, _ = GS.tp53_columns(Q.cov, tcode)
        full, elig = GG.codes_full(Q, fl, cd)
        for term in ("factY", "factX"):
            rows = GG.run_calibration(panel_args, Q, tpX, full, elig, "factorial", term, calB, prim["rule"], range(n_use), WORKERS)
            sm = GG.summarize(rows, cfg.const("calibration_gate"))
            note = ""
            if name == "level5_original":
                a = pd.DataFrame(rows).set_index("pseudo").sort_index()
                b = conf[conf.term == term].set_index("pseudo").sort_index().loc[a.index]
                eq = all(np.allclose(a[c].astype(float), b[c].astype(float), rtol=1e-9, atol=0, equal_nan=True)
                         for c in ("n_tested", "n_p_lt_005", "n_q_lt_01", "min_p", "panel_stat", "panel_p"))
                note = f"validation against results/G/G_calibration.tsv: {'all rows equal' if eq else 'ROWS DIFFER'}"
                assert eq, "the calibration code does not reproduce the confirmatory rows; stop"
            for r in rows:
                r.update(tensor=name, label=LABEL)
            cal_rows += rows
            cal_summ.append(dict(tensor=name, term=term, label=LABEL, **sm, note=note))
            log(f"calibration {name} {term}: mean share {sm.get('mean_share_p_lt_005', float('nan')):.4f}, "
                f"p95 count {sm.get('p95_n_p_lt_005', float('nan'))}, in range {sm['gate_pass']}; {note}")
    pd.DataFrame(cal_rows).to_csv(os.path.join(OUT, "X1_calibration_rows.tsv"), sep="\t", index=False)
    cs = pd.DataFrame(cal_summ)
    cs.to_csv(os.path.join(OUT, "X1_calibration_summary.tsv"), sep="\t", index=False)
    log("calibration written before any real label was read in X1")

    # ---- real labels (after calibration)
    given = pd.read_csv(os.path.join(ROOT, "data/G/inputs/taskD/D_line_table.tsv"), sep="\t")
    act = I.activity()
    h1_rows, cpd = [], []
    for name in ("X1a_level3_plateDMSO", "X1b_level5_minus_lineDMSO"):
        Q, _ = panels[name]
        al = GS.Aligned(Q, given)
        ctx = Ctx()
        ctx.cfg, ctx.P = cfg, {"P1": Q}
        ctx.s_level5 = {"P1": C.response_size(Q.R, Q.avail)}
        ctx.activity = {"P1": dict(zip(act[act.panel == "P1"].cmap_name, act[act.panel == "P1"].above_panel_median.astype(bool)))}
        for term in ("factY", "factX"):
            d = al.design(term, prim)
            r = GA.run_main_term(ctx, "P1", term, d, cfg.run("B_full"), label=f"{LABEL} {name}")
            pr = dict(r["panel_row"])
            pr.pop("null", None)
            c = cs[(cs.tensor == name) & (cs.term == term)].iloc[0]
            pr.update(tensor=name, label=LABEL, calib_mean_share=c.mean_share_p_lt_005, calib_in_range=bool(c.gate_pass),
                      calib_p95_count=c.p95_n_p_lt_005)
            if term == "factX":
                if not pr["calib_in_range"]:
                    reading = "uncalibrated: not read"
                elif pr["n_p_lt_005"] <= pr["calib_p95_count"] and pr["panel_p"] >= 0.05:
                    reading = "the level-5 X-term excess does not persist when responses are referenced to vehicle wells"
                else:
                    reading = "the X-term pattern persists after vehicle referencing"
                pr["reading_fixed_in_plan"] = reading
            h1_rows.append(pr)
            if r["per_cpd"] is not None:
                cpd.append(r["per_cpd"].drop(columns=[x for x in r["per_cpd"].columns if x.startswith("_")], errors="ignore")
                           .assign(tensor=name, label=LABEL))
            log(f"real labels {name} {term}: panel p {pr.get('panel_p', float('nan')):.4f}, n p<0.05 {pr.get('n_p_lt_005')}, "
                f"calibration p95 {pr['calib_p95_count']}")
    h1 = pd.DataFrame(h1_rows)
    h1.to_csv(os.path.join(OUT, "X1_H1.tsv"), sep="\t", index=False)
    cp = pd.concat(cpd, ignore_index=True)
    cp.to_csv(os.path.join(OUT, "X1_per_compound.tsv"), sep="\t", index=False)
    # the 11 X-term compounds at p < 0.05 in the confirmatory run
    h2 = pd.read_csv(os.path.join(ROOT, "results/G/G_H2_per_compound.tsv"), sep="\t")
    c11 = h2[(h2.panel == "P1") & (h2.term == "factX") & (h2.p < 0.05) & (h2.compound != "DMSO (pseudo-compound)")][["compound", "p"]]
    c11 = c11.rename(columns={"p": "p_confirmatory_level5"})
    for name in ("X1a_level3_plateDMSO", "X1b_level5_minus_lineDMSO"):
        x = cp[(cp.tensor == name) & (cp.term == "factX")][["compound", "p"]].rename(columns={"p": f"p_{name}"})
        c11 = c11.merge(x, on="compound", how="left", validate="1:1")
    assert len(c11) == int(h2[(h2.panel == "P1") & (h2.term == "factX") & (h2.compound != "DMSO (pseudo-compound)")].p.lt(0.05).sum())
    c11.assign(label=LABEL).to_csv(os.path.join(OUT, "X1_confirmatory_X_compounds.tsv"), sep="\t", index=False)
    log("X1 done")


if __name__ == "__main__":
    main()
