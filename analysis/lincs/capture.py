"""Capture numerical fits and panel nulls from the registered model without changing its returns."""
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = "1"
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / "src/G"))
import numpy as np
import pandas as pd
import gcore
import ganalyses

CTX, STATE, MANI = [], {"n": 0, "dir": None}, []


def _label():
    return CTX[-1] if CTX else None


def _fname(i, lab, kind):
    safe = lab.replace("|", "__").replace(" ", "_").replace("/", "-").replace("(", "").replace(")", "").replace(",", "")
    return f"{i:04d}_{safe}_{kind}"


_orig_run_term = gcore.run_term
_orig_finish = gcore.finish_table
_orig_rmt = ganalyses.run_main_term
_orig_exp = ganalyses.exposure_run


def run_term(P, Xo, k, mask, B, *a, **kw):
    df, permT, recs = _orig_run_term(P, Xo, k, mask, B, *a, **kw)
    lab = _label()
    if lab is None:
        return df, permT, recs
    STATE["n"] += 1
    i = STATE["n"]
    G = len(P.genes)
    arr = {n: np.full((len(recs), G), np.nan, np.float32) for n in ("beta", "se_hc3", "t", "resid_sd")}
    for j, r in enumerate(recs):
        if r.get("tested") and "_b" in r:
            arr["beta"][j], arr["se_hc3"][j], arr["t"][j], arr["resid_sd"][j] = r["_b"], r["_se"], r["_t"], r["_sig"]
    fn = _fname(i, lab, "run_term.npz")
    np.savez_compressed(os.path.join(STATE["dir"], fn), compounds=np.array([r["compound"] for r in recs]),
                        tested=np.array([bool(r.get("tested")) for r in recs]), genes=np.array(list(P.genes)),
                        permT=(np.asarray(permT, dtype=np.float32) if permT is not None else np.empty((0, 0), np.float32)),
                        B=np.int64(B), **arr)
    df.drop(columns=[c for c in df.columns if c.startswith("_")], errors="ignore").to_csv(
        os.path.join(STATE["dir"], _fname(i, lab, "run_term_table.tsv")), sep="\t", index=False)
    MANI.append(dict(i=i, context=lab, kind="run_term", file=fn, n_records=len(recs),
                     n_tested=int(sum(bool(r.get("tested")) for r in recs)),
                     permT_shape="" if permT is None else "x".join(map(str, np.shape(permT)))))
    return df, permT, recs


def finish_table(df, permT, *a, **kw):
    out = _orig_finish(df, permT, *a, **kw)
    lab = _label()
    if lab is not None:
        STATE["n"] += 1
        i = STATE["n"]
        fn = _fname(i, lab, "panel_null.npz")
        np.savez_compressed(os.path.join(STATE["dir"], fn), panel_stat=np.float64(out[1]), panel_p=np.float64(out[2]),
                            null=np.asarray(out[3], dtype=np.float64))
        MANI.append(dict(i=i, context=lab, kind="finish_table", file=fn, n_records=len(out[0]), n_tested=int(out[0].in_ranking.sum()),
                         permT_shape=str(len(np.asarray(out[3])))))
    return out


def run_main_term(ctx, pn, term, d, B, keep_perm=False, label="primary"):
    CTX.append(f"main|{label}|{pn}|{term}")
    try:
        return _orig_rmt(ctx, pn, term, d, B, keep_perm=keep_perm, label=label)
    finally:
        CTX.pop()


def exposure_run(ctx, pn, name, *a, **kw):
    CTX.append(f"exposure|{name}|{pn}")
    try:
        return _orig_exp(ctx, pn, name, *a, **kw)
    finally:
        CTX.pop()


def record(directory):
    """Record fit arrays and permutation nulls using the same calls that produced the tables."""
    STATE["dir"] = str(directory)
    directory.mkdir(parents=True, exist_ok=False)
    gcore.run_term, gcore.finish_table = run_term, finish_table
    ganalyses.run_main_term, ganalyses.exposure_run = run_main_term, exposure_run


def finish():
    pd.DataFrame(MANI).to_csv(os.path.join(STATE["dir"], "capture_manifest.tsv"), sep="\t", index=False)
