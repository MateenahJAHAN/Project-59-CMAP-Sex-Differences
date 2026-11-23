#!/usr/bin/env python3
"""Step 2 (A1, A2): map identifiers, set cut-offs, call states, tabulate mismatches."""
import json
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde, spearmanr
from sklearn.metrics import cohen_kappa_score
from sklearn.mixture import GaussianMixture

sys.path.insert(0, os.path.dirname(__file__))
from common import (CONF, DERIVED, MARKERS, RAW, RES, SEED, STATES, XIST, Y_GENES, log, norm,  # noqa: E402
                    timed)


# ---------------------------------------------------------------- mapping
def unique_match(cands):
    cands = sorted(set(c for c in cands if c))
    if len(cands) == 1:
        return cands[0], ""
    if len(cands) > 1:
        return None, "multiple:" + "|".join(cands)
    return None, ""


def build_mapping(cell):
    ccle_cols = pd.read_csv(os.path.join(DERIVED, "ccle_columns.tsv"), sep="\t")["ccle_name"].tolist()
    ccle_set = set(ccle_cols)
    ccle_by_norm = {}
    for c in ccle_cols:
        ccle_by_norm.setdefault(norm(c.split("_", 1)[0]), []).append(c)
    ml = pd.read_csv(os.path.join(RAW, "model_list_20240110.csv"), low_memory=False)
    assert ml["model_id"].is_unique
    by_rrid, by_cclei, by_name, by_syn = {}, {}, {}, {}
    for r in ml.itertuples():
        if isinstance(r.RRID, str):
            by_rrid.setdefault(r.RRID, []).append(r.model_id)
        if isinstance(r.CCLE_ID, str):
            by_cclei.setdefault(r.CCLE_ID, []).append(r.model_id)
        by_name.setdefault(norm(r.model_name), []).append(r.model_id)
        if isinstance(r.synonyms, str):
            for s in re.split(r"[;,]", r.synonyms):
                if norm(s):
                    by_syn.setdefault(norm(s), []).append(r.model_id)
    ml_i = ml.set_index("model_id")
    cs = pd.read_csv(os.path.join(DERIVED, "cellosaurus_compact.tsv"), sep="\t", keep_default_na=False)
    cs_ac = cs.set_index("AC")
    cs_by_norm = {}
    for r in cs[cs["human"].astype(str) == "True"].itertuples():
        cs_by_norm.setdefault(norm(r.ID), []).append(r.AC)

    out = []
    for r in cell.itertuples():
        notes = []
        cvcl = r.cellosaurus_id if isinstance(r.cellosaurus_id, str) and r.cellosaurus_id else None
        cn = r.ccle_name if isinstance(r.ccle_name, str) and r.ccle_name else None
        nn = norm(r.cell_iname)
        # --- Sanger model
        sm, skey = None, ""
        if cvcl:
            m, why = unique_match(by_rrid.get(cvcl, []))
            if m:
                sm, skey = m, "cellosaurus_id=RRID"
            elif why:
                notes.append(f"sanger RRID {why}")
        if sm is None and cn:
            m, why = unique_match(by_cclei.get(cn, []))
            if m:
                sm, skey = m, "ccle_name=CCLE_ID"
            elif why:
                notes.append(f"sanger CCLE_ID {why}")
        if sm is None and nn:
            m, why = unique_match(by_name.get(nn, []))
            if m:
                sm, skey = m, "name=model_name"
            elif why:
                notes.append(f"sanger name {why}")
            else:
                m, why = unique_match(by_syn.get(nn, []))
                if m:
                    sm, skey = m, "name=synonym"
                elif why:
                    notes.append(f"sanger synonym {why}")
        # --- CCLE column
        cc, ckey = None, ""
        if cn and cn in ccle_set:
            cc, ckey = cn, "ccle_name"
        if cc is None and cvcl:
            cands = [ml_i.at[m, "CCLE_ID"] for m in by_rrid.get(cvcl, [])
                     if isinstance(ml_i.at[m, "CCLE_ID"], str) and ml_i.at[m, "CCLE_ID"] in ccle_set]
            m, why = unique_match(cands)
            if m:
                cc, ckey = m, "cellosaurus_id->RRID->CCLE_ID"
            elif why:
                notes.append(f"ccle via RRID {why}")
        if cc is None and nn:
            m, why = unique_match(ccle_by_norm.get(nn, []))
            if m:
                cc, ckey = m, "name"
            elif why:
                notes.append(f"ccle name {why}")
        # --- Cellosaurus
        ac, akey = None, ""
        if cvcl and cvcl in cs_ac.index:
            ac, akey = cvcl, "cellosaurus_id"
        elif cvcl:
            notes.append(f"cellosaurus_id {cvcl} not in cellosaurus.txt")
        if ac is None and sm is not None and isinstance(ml_i.at[sm, "RRID"], str) and ml_i.at[sm, "RRID"] in cs_ac.index:
            ac, akey = ml_i.at[sm, "RRID"], "sanger RRID"
        if ac is None and nn:
            m, why = unique_match(cs_by_norm.get(nn, []))
            if m:
                ac, akey = m, "name=ID"
            elif why:
                notes.append(f"cellosaurus name {why}")
        out.append(dict(
            cell_iname=r.cell_iname, ccle_name=cc, ccle_map_key=ckey, sanger_model_id=sm, sanger_map_key=skey,
            sanger_model_name=ml_i.at[sm, "model_name"] if sm else None,
            sanger_gender=ml_i.at[sm, "gender"] if sm else None,
            cellosaurus_ac=ac, cellosaurus_map_key=akey,
            cellosaurus_name=cs_ac.at[ac, "ID"] if ac else None,
            cellosaurus_sex=(cs_ac.at[ac, "SX"] or None) if ac else None,
            map_notes="; ".join(notes)))
    mp = pd.DataFrame(out)
    assert len(mp) == len(cell) and mp["cell_iname"].is_unique
    return mp


# ---------------------------------------------------------------- cut-offs
def bimodal_cutoff(x, label):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    grid = np.linspace(x.min(), x.max(), 1024)
    dens = gaussian_kde(x)(grid)  # Scott bandwidth (scipy default)
    mx = argrelextrema(dens, np.greater)[0]
    info = dict(axis=label, n=len(x), method="kde_density_minimum", n_modes=len(mx))
    if len(mx) >= 2:
        top2 = sorted(mx[np.argsort(dens[mx])[-2:]])
        seg = np.arange(top2[0], top2[1] + 1)
        imin = seg[np.argmin(dens[seg])]
        cut = grid[imin]
        info.update(mode_low=grid[top2[0]], mode_high=grid[top2[1]], cutoff=cut,
                    density_at_cutoff=dens[imin], density_mode_low=dens[top2[0]], density_mode_high=dens[top2[1]],
                    dip_ratio=dens[imin] / min(dens[top2[0]], dens[top2[1]]))
    else:
        gm = GaussianMixture(2, random_state=SEED).fit(x.reshape(-1, 1))
        lo, hi = sorted(gm.means_.ravel())
        g2 = np.linspace(lo, hi, 2000)
        post = gm.predict_proba(g2.reshape(-1, 1))
        hi_comp = int(np.argmax(gm.means_.ravel()))
        cut = g2[np.argmin(np.abs(post[:, hi_comp] - 0.5))]
        info.update(method="gmm_fallback", mode_low=lo, mode_high=hi, cutoff=cut)
    info["n_high"] = int((x >= info["cutoff"]).sum())
    info["n_low"] = int((x < info["cutoff"]).sum())
    return info, grid, dens


def call_state(ann, yhi, xhi):
    if pd.isna(yhi) or pd.isna(xhi):
        return "no_data"
    if ann == "M":
        return "XY" if yhi else ("discordant" if xhi else "Y-lost")
    if ann == "F":
        return "discordant" if yhi else ("XX-like" if xhi else "XIST-lost")
    return "XY" if yhi else ("XX-like" if xhi else "ambiguous")


def source_scores(df):
    ygs = [g for g in Y_GENES if g in df.columns]
    df = df.copy()
    df["y_score"] = df[ygs].mean(axis=1)
    df["n_y_genes"] = df[ygs].notna().sum(axis=1)
    return df, ygs


def main():
    cell = pd.read_csv(os.path.join(RAW, "cellinfo_beta.txt"), sep="\t")
    assert cell["cell_iname"].is_unique and len(cell) == 240
    assert set(cell["donor_sex"]) <= {"M", "F", "Unknown"}
    with timed("states:mapping"):
        mp = build_mapping(cell)

    # ---------------- expression
    ccle = pd.read_csv(os.path.join(DERIVED, "ccle_markers.tsv"), sep="\t")
    sang = pd.read_csv(os.path.join(DERIVED, "sanger_markers.tsv"), sep="\t")
    brd = pd.read_csv(os.path.join(DERIVED, "cmp_broad_markers.tsv"), sep="\t")
    assert ccle["ccle_name"].is_unique and sang["model_id"].is_unique and brd["model_id"].is_unique
    ccle, ccle_y = source_scores(ccle)
    sang, sang_y = source_scores(sang)
    brd, _ = source_scores(brd)
    log("states", f"CCLE Y genes {ccle_y}; Sanger Y genes {sang_y}")

    with timed("states:cutoffs"):
        cuts, dens_rows = [], []
        for src, df in (("ccle", ccle), ("sanger", sang), ("cmpbroad", brd)):
            for ax, col in (("XIST", XIST), ("y_score", "y_score")):
                info, grid, dens = bimodal_cutoff(df[col], ax)
                info["source"] = src
                cuts.append(info)
                dens_rows.append(pd.DataFrame({"source": src, "axis": ax, "x": grid, "density": dens}))
        cuts = pd.DataFrame(cuts)[["source", "axis", "n", "method", "n_modes", "mode_low", "mode_high", "cutoff",
                                   "density_at_cutoff", "density_mode_low", "density_mode_high", "dip_ratio",
                                   "n_high", "n_low"]]
        cuts.to_csv(os.path.join(RES, "A_cutoffs.tsv"), sep="\t", index=False)
        pd.concat(dens_rows).to_csv(os.path.join(DERIVED, "kde_curves.tsv"), sep="\t", index=False)
        log("states", "cut-offs:\n" + cuts[["source", "axis", "method", "n_modes", "cutoff", "dip_ratio"]].to_string())
    C = {(r.source, r.axis): r.cutoff for r in cuts.itertuples()}

    # ---------------- assemble per-line table
    cnt = pd.read_csv(os.path.join(DERIVED, "trt_cp_counts.tsv"), sep="\t")
    n0 = len(cell)
    t = cell[["cell_iname", "cell_lineage", "cell_type", "donor_sex", "cellosaurus_id"]].merge(mp, on="cell_iname", how="left")
    t = t.merge(cnt, on="cell_iname", how="left")
    t["n_trt_cp_signatures"] = t["n_trt_cp_signatures"].fillna(0).astype(int)
    assert len(t) == n0 and t["cell_iname"].is_unique
    cc = ccle.rename(columns={g: f"{g}_ccle" for g in MARKERS}).rename(
        columns={"y_score": "y_score_ccle", "n_y_genes": "n_y_genes_ccle"})
    cc = cc[["ccle_name"] + [f"{g}_ccle" for g in MARKERS] + ["y_score_ccle", "n_y_genes_ccle"]]
    t = t.merge(cc, on="ccle_name", how="left")
    assert len(t) == n0
    ss = sang.rename(columns={g: f"{g}_sanger" for g in MARKERS}).rename(
        columns={"y_score": "y_score_sanger", "n_y_genes": "n_y_genes_sanger", "model_id": "sanger_model_id"})
    ss = ss[["sanger_model_id"] + [f"{g}_sanger" for g in MARKERS if f"{g}_sanger" in ss.columns]
            + ["y_score_sanger", "n_y_genes_sanger", "n_sanger_datasets"]]
    t = t.merge(ss, on="sanger_model_id", how="left")
    assert len(t) == n0 and t["cell_iname"].is_unique
    bb = brd.rename(columns={"model_id": "sanger_model_id", XIST: "XIST_cmpbroad", "y_score": "y_score_cmpbroad"})
    t = t.merge(bb[["sanger_model_id", "XIST_cmpbroad", "y_score_cmpbroad"]], on="sanger_model_id", how="left")
    assert len(t) == n0 and t["cell_iname"].is_unique
    # Deviation 1: the Broad-sourced rows of the Sanger file are a fallback only, for lines with no CCLE legacy
    # and no Sanger-sequenced data.
    has_primary = t["XIST_ccle"].notna() | t["XIST_sanger"].notna()
    t.loc[has_primary, ["XIST_cmpbroad", "y_score_cmpbroad"]] = np.nan

    # ---------------- states, with the 3x3 sensitivity grid
    shifts = [-1.0, 0.0, 1.0]
    grid_rows = []
    for src in ("ccle", "sanger", "cmpbroad"):
        x = t[f"XIST_{src}"]
        y = t[f"y_score_{src}"]
        for dy in shifts:
            for dx in shifts:
                yhi = (y >= C[(src, "y_score")] + dy).where(y.notna())
                xhi = (x >= C[(src, "XIST")] + dx).where(x.notna())
                st = [call_state(a, b, c) for a, b, c in zip(t["donor_sex"], yhi, xhi)]
                col = f"state_{src}_dy{dy:+.0f}_dx{dx:+.0f}"
                t[col] = st
                if dy == 0 and dx == 0:
                    t[f"state_{src}"] = st
                    t[f"Yhi_{src}"] = yhi
                    t[f"Xhi_{src}"] = xhi
    # combined
    states, confs, notes = [], [], []
    for _, r in t.iterrows():
        avail = [s for s in ("ccle", "sanger", "cmpbroad") if r[f"state_{s}"] != "no_data"]
        nt = [r["map_notes"]] if isinstance(r["map_notes"], str) and r["map_notes"] else []
        if not avail:
            st, cf = "no_data", "none"
        else:
            calls = {s: r[f"state_{s}"] for s in avail}
            if len(set(calls.values())) > 1:
                st, cf = "ambiguous", "low"
                nt.append("sources disagree: " + ", ".join(f"{k}={v}" for k, v in calls.items()))
            else:
                st = list(calls.values())[0]
                if st == "ambiguous":
                    cf = "low"
                else:
                    stable = []
                    for s in avail:
                        g = {r[f"state_{s}_dy{dy:+.0f}_dx{dx:+.0f}"] for dy in shifts for dx in shifts}
                        stable.append(len(g) == 1)
                    cf = "high" if any(stable) else "medium"
                if len(avail) == 1:
                    nt.append(f"single source: {avail[0]}")
        states.append(st)
        confs.append(cf)
        notes.append("; ".join(nt))
    t["state"] = states
    t["confidence"] = confs
    t["notes"] = notes
    assert set(t["state"]) <= set(STATES), set(t["state"])
    assert set(t["confidence"]) <= set(CONF)
    t["y_score"] = t["y_score_ccle"].where(t["y_score_ccle"].notna(), t["y_score_sanger"])
    t["y_score"] = t["y_score"].where(t["y_score"].notna(), t["y_score_cmpbroad"])
    t["y_score_source"] = np.where(t["y_score_ccle"].notna(), "ccle",
                                   np.where(t["y_score_sanger"].notna(), "sanger",
                                            np.where(t["y_score_cmpbroad"].notna(), "cmpbroad", "")))
    # Y high / XIST high in the combined call (for relabelling): sources must agree
    def comb(r, k):
        v = [r[f"{k}_{s}"] for s in ("ccle", "sanger", "cmpbroad") if pd.notna(r[f"{k}_{s}"])]
        if not v:
            return np.nan
        return bool(v[0]) if len(set(v)) == 1 else np.nan
    t["Yhi"] = t.apply(lambda r: comb(r, "Yhi"), axis=1)
    t["Xhi"] = t.apply(lambda r: comb(r, "Xhi"), axis=1)

    cols = (["cell_iname", "cell_lineage", "cell_type", "n_trt_cp_signatures", "donor_sex", "cellosaurus_sex",
             "sanger_gender", "ccle_name", "sanger_model_id", "sanger_model_name", "cellosaurus_id", "cellosaurus_ac",
             "cellosaurus_name", "ccle_map_key", "sanger_map_key", "cellosaurus_map_key", "n_sanger_datasets"]
            + [f"XIST_{s}" for s in ("ccle", "sanger")]
            + [f"{g}_{s}" for g in Y_GENES for s in ("ccle", "sanger")]
            + ["y_score_ccle", "y_score_sanger", "XIST_cmpbroad", "y_score_cmpbroad", "y_score", "y_score_source",
               "Yhi_ccle", "Xhi_ccle", "Yhi_sanger", "Xhi_sanger", "Yhi_cmpbroad", "Xhi_cmpbroad", "Yhi", "Xhi",
               "state_ccle", "state_sanger", "state_cmpbroad", "state", "confidence", "notes"])
    out = t[cols].rename(columns={"cell_lineage": "lineage", "donor_sex": "cmap_donor_sex"})
    out.to_csv(os.path.join(RES, "A_cellline_sex_states.tsv"), sep="\t", index=False, float_format="%.4f")
    grid_cols = ["cell_iname"] + [c for c in t.columns if re.match(r"state_(ccle|sanger|cmpbroad)_dy", c)]
    t[grid_cols].to_csv(os.path.join(RES, "A_state_sensitivity_per_line.tsv"), sep="\t", index=False)

    # ---------------- summaries
    with timed("states:summaries"):
        summaries(t, out, C, ccle, sang)


def crosstab_save(a, b, name):
    ct = pd.crosstab(a.fillna("NA"), b.fillna("NA"))
    ct.to_csv(os.path.join(RES, name), sep="\t")
    return ct


def summaries(t, out, C, ccle, sang):
    R = {}
    # mapping summary
    rows = []
    for k in ("ccle_map_key", "sanger_map_key", "cellosaurus_map_key"):
        vc = t[k].fillna("").replace("", "unmapped").value_counts()
        for key, n in vc.items():
            rows.append(dict(target=k.replace("_map_key", ""), key=key, n_lines=int(n)))
    rows.append(dict(target="any_rnaseq", key="CCLE or Sanger RNA-seq available",
                     n_lines=int((t["state_ccle"].ne("no_data") | t["state_sanger"].ne("no_data")
                                  | t["state_cmpbroad"].ne("no_data")).sum())))
    rows.append(dict(target="ccle_rnaseq", key="CCLE marker data", n_lines=int(t["XIST_ccle"].notna().sum())))
    rows.append(dict(target="sanger_rnaseq", key="Sanger (data_source=Sanger) marker data",
                     n_lines=int(t["XIST_sanger"].notna().sum())))
    rows.append(dict(target="cmpbroad_rnaseq", key="fallback: Sanger file, data_source=Broad (only lines without CCLE/Sanger)",
                     n_lines=int(t["XIST_cmpbroad"].notna().sum())))
    rows.append(dict(target="both_rnaseq", key="both sources",
                     n_lines=int((t["XIST_ccle"].notna() & t["XIST_sanger"].notna()).sum())))
    pd.DataFrame(rows).to_csv(os.path.join(RES, "A_mapping_summary.tsv"), sep="\t", index=False)
    um = t[(t["state"] == "no_data")][["cell_iname", "cell_type", "donor_sex", "cell_lineage", "n_trt_cp_signatures",
                                       "ccle_name", "sanger_model_id", "cellosaurus_ac", "map_notes"]]
    um.to_csv(os.path.join(RES, "A_lines_without_rnaseq.tsv"), sep="\t", index=False)
    unm = t[t["ccle_name"].isna() & t["sanger_model_id"].isna()][
        ["cell_iname", "cell_type", "donor_sex", "cell_lineage", "n_trt_cp_signatures", "cellosaurus_id", "map_notes"]]
    unm.to_csv(os.path.join(RES, "A_unmapped_lines.tsv"), sep="\t", index=False)
    mapped = t[t["ccle_name"].notna() | t["sanger_model_id"].notna()]
    spot = mapped.sample(n=min(10, len(mapped)), random_state=SEED)[
        ["cell_iname", "ccle_name", "ccle_map_key", "sanger_model_name", "sanger_map_key", "cellosaurus_name",
         "cellosaurus_map_key"]]
    spot.to_csv(os.path.join(RES, "A_mapping_spotcheck.tsv"), sep="\t", index=False)

    # annotation agreement
    crosstab_save(t["donor_sex"], t["cellosaurus_sex"], "A_annot_cmap_vs_cellosaurus.tsv")
    crosstab_save(t["donor_sex"], t["sanger_gender"], "A_annot_cmap_vs_sanger.tsv")

    # state counts overall and per source
    sc = pd.crosstab(t["state"], t["donor_sex"], margins=True)
    sc.to_csv(os.path.join(RES, "A_state_counts.tsv"), sep="\t")
    for s in ("ccle", "sanger", "cmpbroad"):
        crosstab_save(t[f"state_{s}"], t["donor_sex"], f"A_state_counts_{s}.tsv")
    pd.crosstab(t["state"], t["confidence"], margins=True).to_csv(os.path.join(RES, "A_state_by_confidence.tsv"), sep="\t")

    # sensitivity grid summary
    shifts = [-1, 0, 1]
    rows = []
    for s in ("ccle", "sanger", "cmpbroad"):
        for dy in shifts:
            for dx in shifts:
                col = f"state_{s}_dy{dy:+d}_dx{dx:+d}"
                vc = t[col].value_counts()
                base = t[f"state_{s}"]
                has = base != "no_data"
                rows.append(dict(source=s, y_shift=dy, xist_shift=dx, y_cutoff=C[(s, "y_score")] + dy,
                                 xist_cutoff=C[(s, "XIST")] + dx,
                                 n_lines_with_data=int(has.sum()),
                                 n_changed_vs_base=int((t[col] != base)[has].sum()),
                                 **{st: int(vc.get(st, 0)) for st in STATES if st != "no_data"}))
    pd.DataFrame(rows).to_csv(os.path.join(RES, "A_state_sensitivity_grid.tsv"), sep="\t", index=False)

    # CCLE vs Sanger agreement (CMap lines with both)
    both = t[t["XIST_ccle"].notna() & t["XIST_sanger"].notna()]
    ag = []
    for k in ("Yhi", "Xhi", "state"):
        a = both[f"{k}_ccle"].astype(str)
        b = both[f"{k}_sanger"].astype(str)
        ag.append(dict(set="CMap lines with both sources", measure=k, n=len(both), n_agree=int((a == b).sum()),
                       pct_agree=100 * (a == b).mean() if len(both) else np.nan,
                       cohen_kappa=cohen_kappa_score(a, b) if len(both) > 1 else np.nan))
    for k, c in (("XIST", "XIST"), ("y_score", "y_score")):
        rho = spearmanr(both[f"{c}_ccle"], both[f"{c}_sanger"]).statistic if len(both) > 2 else np.nan
        ag.append(dict(set="CMap lines with both sources", measure=f"spearman_{k}", n=len(both), spearman_rho=rho))
    # secondary: all models in both files linked by Sanger CCLE_ID
    ml = pd.read_csv(os.path.join(RAW, "model_list_20240110.csv"), low_memory=False)[["model_id", "CCLE_ID"]]
    s2 = sang.merge(ml, on="model_id", how="left")
    assert len(s2) == len(sang)
    s2 = s2.dropna(subset=["CCLE_ID"]).merge(ccle, left_on="CCLE_ID", right_on="ccle_name", suffixes=("_s", "_c"))
    assert s2["model_id"].is_unique and s2["ccle_name"].is_unique
    for k, col in (("Yhi", "y_score"), ("Xhi", "XIST")):
        a = s2[f"{col}_c"] >= C[("ccle", "y_score" if col == "y_score" else "XIST")]
        b = s2[f"{col}_s"] >= C[("sanger", "y_score" if col == "y_score" else "XIST")]
        ag.append(dict(set="all models in both files (Sanger CCLE_ID)", measure=k, n=len(s2),
                       n_agree=int((a == b).sum()), pct_agree=100 * (a == b).mean(),
                       cohen_kappa=cohen_kappa_score(a, b),
                       spearman_rho=spearmanr(s2[f"{col}_c"], s2[f"{col}_s"]).statistic))
    pd.DataFrame(ag).to_csv(os.path.join(RES, "A_ccle_sanger_agreement.tsv"), sep="\t", index=False)
    s2[["model_id", "ccle_name", "XIST_c", "XIST_s", "y_score_c", "y_score_s"]].to_csv(
        os.path.join(DERIVED, "ccle_sanger_overlap_all.tsv"), sep="\t", index=False)
    if len(both):
        crosstab_save(both["state_ccle"], both["state_sanger"], "A_state_ccle_vs_sanger.tsv")

    # ---------------- A2 mismatch
    strict = {"discordant"}
    broad = {"discordant", "Y-lost", "XIST-lost"}
    ann = t[t["donor_sex"].isin(["M", "F"])].copy()
    ann["called"] = ann["state"] != "no_data"
    ann["is_amb"] = ann["state"] == "ambiguous"
    ann["den"] = ann["called"] & ~ann["is_amb"]
    ann["strict"] = ann["state"].isin(strict)
    ann["broad"] = ann["state"].isin(broad)

    def summarise(g):
        d = g[g["den"]]
        sig_d = d["n_trt_cp_signatures"].sum()
        return pd.Series(dict(
            n_lines_annotated=len(g), n_with_state=int(g["called"].sum()), n_ambiguous=int(g["is_amb"].sum()),
            n_denominator=len(d),
            n_XY=int((d.state == "XY").sum()), n_XX_like=int((d.state == "XX-like").sum()),
            n_Y_lost=int((d.state == "Y-lost").sum()), n_XIST_lost=int((d.state == "XIST-lost").sum()),
            n_discordant=int((d.state == "discordant").sum()),
            n_strict_mismatch=int(d["strict"].sum()), n_broad_mismatch=int(d["broad"].sum()),
            pct_strict_mismatch=100 * d["strict"].mean() if len(d) else np.nan,
            pct_broad_mismatch=100 * d["broad"].mean() if len(d) else np.nan,
            trt_cp_sigs_denominator=int(sig_d),
            trt_cp_sigs_strict_mismatch=int(d.loc[d["strict"], "n_trt_cp_signatures"].sum()),
            trt_cp_sigs_broad_mismatch=int(d.loc[d["broad"], "n_trt_cp_signatures"].sum()),
            pct_sigs_strict_mismatch=100 * d.loc[d["strict"], "n_trt_cp_signatures"].sum() / sig_d if sig_d else np.nan,
            pct_sigs_broad_mismatch=100 * d.loc[d["broad"], "n_trt_cp_signatures"].sum() / sig_d if sig_d else np.nan))

    rows = []
    allrow = summarise(ann)
    rows.append(pd.concat([pd.Series(dict(breakdown="all", group="all", annotated_sex="M+F")), allrow]))
    for sx in ("M", "F"):
        rows.append(pd.concat([pd.Series(dict(breakdown="all", group="all", annotated_sex=sx)),
                               summarise(ann[ann.donor_sex == sx])]))
    for lin, g in ann.groupby("cell_lineage"):
        rows.append(pd.concat([pd.Series(dict(breakdown="cell_lineage", group=lin, annotated_sex="M+F")), summarise(g)]))
    bins = [-1, 0, 99, 999, 9999, np.inf]
    labels = ["0", "1-99", "100-999", "1000-9999", ">=10000"]
    ann["sig_bin"] = pd.cut(ann["n_trt_cp_signatures"], bins=bins, labels=labels)
    for b in labels:
        g = ann[ann["sig_bin"] == b]
        rows.append(pd.concat([pd.Series(dict(breakdown="n_trt_cp_signatures", group=b, annotated_sex="M+F")),
                               summarise(g)]))
    mm = pd.DataFrame(rows)
    mm.to_csv(os.path.join(RES, "A_mismatch_by_lineage.tsv"), sep="\t", index=False, float_format="%.2f")
    # Unknown-annotation lines
    unk = t[t["donor_sex"] == "Unknown"]
    pd.crosstab(unk["state"], unk["cell_type"], margins=True).to_csv(
        os.path.join(RES, "A_unknown_annotation_states.tsv"), sep="\t")
    # states by lineage (all lines) for figure
    pd.crosstab(t["cell_lineage"], t["state"]).to_csv(os.path.join(RES, "A_state_counts_by_lineage.tsv"), sep="\t")
    log("states", "state counts:\n" + sc.to_string())
    log("states", "mismatch (all):\n" + mm.iloc[:3].T.to_string())


if __name__ == "__main__":
    with timed("states:all"):
        main()
