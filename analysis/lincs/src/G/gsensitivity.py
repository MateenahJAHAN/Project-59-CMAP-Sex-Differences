"""Sensitivity analyses S1 to S12 of task G (REC 4g, FROZEN; draft section 4), for the primary panel P1 and the
factorial Y and X terms (S2 also in P2, where its low-QC lines are). Each variant changes one thing:

S1  other line sets for the factorial contrast: task D's 44 lines (7p default for 7f A), REC's exclusion
    (41 lines), and REC's literal S1 (source-conflict lines back in under their Task C call)
S2  drop lines with qc_pass fraction below 0.2
S3  Y-only coding (lineage + TP53 + Y)
S4  P1b (0.1 uM): the P1b rows of G_H1.tsv
S5  log2 doubling time as a covariate          S6  ploidy as a covariate
S7  no TP53 term                               S8  no lineage term (one block, same lines)
S9  CMP cancer type instead of LINCS lineage   S10 TP53 from the CCLE 2018 MAF
S11 signatures with two or more replicates (siginfo nsample >= 2), rebuilt from per-signature level 5
S12 E3 on the other scale (in G_H3.tsv)
S5 to S10 keep the primary analysis set, minus lines that lack the new column.
"""
import copy
import numpy as np, pandas as pd
import gcore as C


def variants(cfg, prim, primary_mask):
    f = cfg.decide("7f")
    v = []
    for lv in ("D44", "F41", "taskC_conflict"):
        if (lv == "D44" and f == "B") or (lv == "F41" and f == "A"):
            continue
        v.append(("S1", f"line set {lv}", dict(prim, line_variant=lv, label=f"S1 {lv}")))
    v.append(("S2", "drop lines with qc_pass fraction < 0.2", dict(prim, qc_drop=True, qc_threshold=cfg.const("qc_low_threshold"), label="S2")))
    v.append(("S3", "Y-only coding", dict(prim, y_only=True, label="S3")))
    v.append(("S5", "log2 doubling time covariate", dict(prim, extra_cov=("log2_doubling_time",), same_lines=primary_mask, label="S5")))
    v.append(("S6", "ploidy covariate", dict(prim, extra_cov=("ploidy",), same_lines=primary_mask, label="S6")))
    v.append(("S7", "no TP53 term", dict(prim, tp53="none", same_lines=primary_mask, label="S7")))
    v.append(("S8", "no lineage term", dict(prim, block="none", same_lines=primary_mask, label="S8")))
    v.append(("S9", "CMP cancer type as block", dict(prim, block="cmp_cancer_type", same_lines=primary_mask, label="S9")))
    v.append(("S10", "TP53 from the CCLE 2018 MAF", dict(prim, tp53="ccle", same_lines=primary_mask, label="S10")))
    return v


def build_s11(P, persig, siginfo, label="1.11 uM", min_nsample=2):
    """P1 tensor from per-signature level 5. Check: with every signature the result must equal task D's P1 tensor
    (max abs difference < 1e-5); then the mean over signatures with nsample >= min_nsample."""
    Z, sig_ids, st, eg = persig
    gp = {g: i for i, g in enumerate(eg)}
    Z = Z[:, [gp[g] for g in P.genes]].astype(np.float64)
    si = siginfo.set_index("sig_id").loc[sig_ids]
    sel = (si.pert_idose == label).values & si.cell_iname.isin(P.lines).values & si.cmap_name.isin(P.comps).values
    li = {l: i for i, l in enumerate(P.lines)}; ci = {c: i for i, c in enumerate(P.comps)}
    ii = si.cell_iname.map(li).values; jj = si.cmap_name.map(ci).values
    ns = pd.to_numeric(si.nsample, errors="coerce").values

    def mean_over(m):
        acc = np.zeros(P.R.shape); cnt = np.zeros(P.R.shape[:2])
        for r in np.flatnonzero(m):
            acc[ii[r], jj[r]] += Z[r]; cnt[ii[r], jj[r]] += 1
        with np.errstate(invalid="ignore"):
            out = acc / cnt[:, :, None]
        out[cnt == 0] = np.nan
        return out, cnt
    allR, call = mean_over(sel)
    both = np.isfinite(allR[:, :, 0]) & np.isfinite(P.R[:, :, 0])
    same_avail = bool((np.isfinite(allR[:, :, 0]) == np.isfinite(P.R[:, :, 0])).all())
    maxdiff = float(np.nanmax(np.abs(allR[both] - P.R[both]))) if both.any() else np.inf
    R11, c11 = mean_over(sel & (ns >= min_nsample))
    info = dict(persig_hash_status=st, n_signatures=int(sel.sum()), n_signatures_nsample_ge2=int((sel & (ns >= min_nsample)).sum()),
                rebuild_equals_D_tensor=bool(same_avail and maxdiff < 1e-5), rebuild_max_abs_diff=maxdiff,
                n_cells=int(np.isfinite(P.R[:, :, 0]).sum()), n_cells_S11=int(np.isfinite(R11[:, :, 0]).sum()))
    assert info["rebuild_equals_D_tensor"], f"S11 builder does not reproduce task D's P1 tensor: {info}"
    Q = copy.copy(P)
    Q.R = R11
    Q.avail = np.isfinite(R11[:, :, 0])
    Q.dmso = None
    return Q, info
