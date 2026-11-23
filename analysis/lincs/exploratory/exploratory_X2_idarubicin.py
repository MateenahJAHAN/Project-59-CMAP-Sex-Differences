"""Exploratory X2 (analysis_plan/G_exploratory_X2_idarubicin.md): the L4 statistic for idarubicin at "12 uM", 24 h.

EXPLORATORY. Reads the 51 idarubicin signatures at this label from the LINCS2020 level-5 gctx by HTTP range reads
(task D's method: fsspec + h5py; landmark genes by task D's gctx_row_index), checks 5 of them against a direct byte-range
read of the contiguous matrix, averages per line, and applies the frozen lead test (ganalyses.lead_test, 7l A) and the
shape test (gcore.run_term) with the P1 design and keys.
Usage (repository root): <task G python> -I src/G_run/exploratory_X2_idarubicin.py
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
import ganalyses as GA

URL = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/level5_beta_trt_cp_n720216x12328.gctx"
OUT = os.path.join(ROOT, "results", "G_exploratory")
DATA = os.path.join(ROOT, "data", "G", "exploratory")
LABEL = "EXPLORATORY X2"


class Ctx:
    pass


def sha_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def log(msg):
    line = f"{pd.Timestamp.now(tz='UTC').strftime('%Y-%m-%dT%H:%M:%SZ')} [X2] {msg}"
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
    key_seed = cfg.run("key_seed_base") + cfg.const("panel_numbers")["P1"]
    lm = I.landmarks()
    assert list(lm.gene_symbol) == list(P1.genes), "landmark order differs from the P1 tensor's"

    si = pd.read_csv(os.path.join(ROOT, "data/G/meta/siginfo_beta.txt"), sep="\t", low_memory=False,
                     usecols=["sig_id", "cmap_name", "pert_idose", "pert_itime", "pert_type", "project_code", "cell_iname", "pert_id"])
    sel = si[(si.cmap_name == "idarubicin") & (si.pert_idose == "12 uM") & (si.pert_itime == "24 h") & (si.pert_type == "trt_cp")
             & (si.project_code == "ASG")].copy()
    assert len(sel) == 51 and sel.cell_iname.nunique() == 49 and sel.cell_iname.isin(set(P1.lines)).all(), (len(sel), sel.cell_iname.nunique())
    log(f"selected {len(sel)} signatures in {sel.cell_iname.nunique()} lines (pert_id {';'.join(sorted(sel.pert_id.unique()))})")

    # ---- range reads (task D's method), then a byte-range check of 5 signatures
    t0 = time.time()
    fs = fsspec.filesystem("https", block_size=4 * 1024 * 1024)
    with fs.open(URL, "rb") as fo, h5py.File(fo, "r") as h:
        col = [x.decode() for x in h["0/META/COL/id"][:]]
        row = [x.decode() for x in h["0/META/ROW/id"][:]]
        ds = h["0/DATA/0/matrix"]
        assert ds.shape == (len(col), len(row)) and ds.id.get_create_plist().get_layout() == h5py.h5d.CONTIGUOUS
        offset, dtype = ds.id.get_offset(), ds.dtype
        pos = {s: i for i, s in enumerate(col)}
        sig_order = sorted(sel.sig_id, key=lambda s: pos[s])
        M = np.stack([ds[pos[s], :] for s in sig_order]).astype(np.float32)
    assert [row[i] for i in lm.gctx_row_index] == [str(g) for g in lm.gene_id], "gctx row ids do not match task D's landmarks"
    Z = M[:, lm.gctx_row_index.to_numpy()]
    nbytes_row = len(row) * np.dtype(dtype).itemsize
    chk = []
    for s in sig_order[:5]:
        start = offset + pos[s] * nbytes_row
        raw = np.frombuffer(fs.cat_file(URL, start=start, end=start + nbytes_row), dtype=dtype)
        chk.append(dict(sig_id=s, equal=bool(np.array_equal(raw, M[sig_order.index(s)]))))
    assert all(c["equal"] for c in chk), "byte-range read differs from the h5py read"
    pd.DataFrame(chk).assign(label=LABEL).to_csv(os.path.join(OUT, "X2_rangeread_check.tsv"), sep="\t", index=False)
    accessed = pd.Timestamp.now(tz="UTC").strftime("%Y-%m-%dT%H:%M:%SZ")
    zp = os.path.join(DATA, "X2_idarubicin_12uM_24h_signatures_landmarks.npy")
    np.save(zp, Z)
    pd.DataFrame({"sig_id": sig_order}).merge(sel, on="sig_id", how="left", validate="1:1").to_csv(
        os.path.join(OUT, "X2_signatures.tsv"), sep="\t", index=False)
    log(f"read {len(sig_order)} signatures in {time.time() - t0:.0f} s; byte-range check 5 of 5 equal")
    with open(os.path.join(ROOT, "provenance/downloads.tsv"), "a") as f:
        f.write("\t".join(["G-exploratory", URL, os.path.relpath(zp, ROOT), str(len(sig_order) * nbytes_row), sha_file(zp), accessed,
                           "range reads of the 51 idarubicin 12 uM 24 h ASG signatures (all 12328 genes read; 978 landmarks kept); "
                           "sha256 is of the saved landmark array"]) + "\n")

    # ---- P1-shaped tensor with only idarubicin filled (task D's rule: mean over a line's signatures)
    j = list(P1.comps).index("idarubicin")
    R = np.full(P1.R.shape, np.nan)
    sl = pd.Series(sel.set_index("sig_id").loc[sig_order].cell_iname.values)
    for l in sorted(sl.unique()):
        R[list(P1.lines).index(l), j, :] = Z[(sl == l).values].astype(np.float64).mean(axis=0)
    tp = os.path.join(DATA, "X2_idarubicin_12uM_P1shape_Dorder.npy")
    np.save(tp, R)
    pd.DataFrame([dict(tensor="X2_idarubicin_12uM_P1shape_Dorder", path=os.path.relpath(tp, ROOT), sha256=sha_file(tp)),
                  dict(tensor="X2_signatures_landmarks", path=os.path.relpath(zp, ROOT), sha256=sha_file(zp))]).to_csv(
        os.path.join(OUT, "X2_tensors.tsv"), sep="\t", index=False)

    # ---- the test (frozen functions)
    panel_args = ("P1", tp, I.file("taskD/tensors/P1_lines.tsv"), I.file("taskD/tensors/P1_compounds.tsv"), P1.genes, cov, key_seed,
                  cfg.run("B_full"), None)
    Q = C.Panel(*panel_args)
    assert int(Q.avail[:, j].sum()) == 49 and int(Q.avail.sum()) == 49
    s = C.response_size(Q.R, Q.avail)
    given = pd.read_csv(os.path.join(ROOT, "data/G/inputs/taskD/D_line_table.tsv"), sep="\t")
    al = GS.Aligned(Q, given)
    d = al.design("factY", prim)
    rows_all = np.flatnonzero(d["mask"])
    rows = rows_all[Q.avail[rows_all, j]]
    kept, dd, why = C.compound_design(d["lin"], d["Xo"], d["k"], rows, int(d["rule"]["min_group"]), int(d["rule"]["min_relab"]))
    tv = d["Xo"][kept, d["k"]] if len(kept) else np.array([])
    counts = dict(n_set_lines_with_data=len(rows), n_lines_after_lineage_rule=len(kept), n_Y1=int(np.nansum(tv == 1)),
                  n_Y0=int(np.nansum(tv == 0)), testable=dd is not None, reason="" if dd is not None else why)
    log(f"design: {counts}")
    out = dict(label=LABEL, compound="idarubicin", dose_label="12 uM", time="24 h", term="factY", alternative="two-sided", **counts)
    if dd is None:
        out.update(reading_fixed_in_plan="not testable")
    else:
        ctx = Ctx()
        ctx.cfg, ctx.P, ctx.s_level5 = cfg, {"P1": Q}, {"P1": s}
        lead = dict(id="X2", compound="idarubicin", term="factY", alternative="two-sided", panels=["P1"])
        lr = GA.lead_test(ctx, "P1", lead, d, cfg.run("B_full"), s)
        out.update({f"lead_{k}": v for k, v in lr.items() if k not in ("lead", "compound", "panel", "term", "alternative")})
        df, permT, recs = C.run_term(Q, d["Xo"], d["k"], d["mask"], cfg.run("B_full"), comps=[j], lin=d["lin"], rule=d["rule"])
        rr = df[df.compound == "idarubicin"].iloc[0]
        out.update(shape_T_obs=float(rr.T_obs) if "T_obs" in df else np.nan, shape_p=float(rr.p) if "p" in df else np.nan,
                   shape_tested=bool(rr.tested))
        e2, p = out.get("lead_E2_est", np.nan), out.get("lead_p", np.nan)
        out["reading_fixed_in_plan"] = ("consistent with L4 at a dose not used before" if (np.isfinite(e2) and e2 > 0 and p < 0.05)
                                        else "L4 is not supported at 12 uM")
    ld = pd.read_csv(os.path.join(ROOT, "results/G/G_leads.tsv"), sep="\t")
    l4 = ld[(ld.lead == "L4") & (ld.panel == "P1")].iloc[0]
    out.update(P1_L4_p=float(l4.p), P1_L4_E2_est=float(l4.E2_est), P1_L4_E2_lo90=float(l4.E2_lo90), P1_L4_E2_hi90=float(l4.E2_hi90),
               P1_L4_n_lines=int(l4.n_lines))
    pd.DataFrame([out]).to_csv(os.path.join(OUT, "X2_result.tsv"), sep="\t", index=False)
    log(f"X2 result: E2 {out.get('lead_E2_est')}, p {out.get('lead_p')}, shape p {out.get('shape_p')}; reading: {out['reading_fixed_in_plan']}")


if __name__ == "__main__":
    main()
