#!/usr/bin/env python3
"""Build cell line labels: annotated sex (CMP + Cellosaurus), expressed sex-chromosome state,
lineage, sex-specific lineage flag and Cellosaurus doubling time.

Outputs:
  results/C/C_cellline_sex_labels.tsv
  results/C/C_sex_specific_lineages.tsv
  results/C/C_expression_cutoffs.tsv
  results/C/C_annotation_discordance.tsv
"""
import os
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (DATA, DROP_CANCER_TYPES, RES, SEX_SPECIFIC_CANCER_TYPES, SEX_SPECIFIC_TISSUES,
                    Y_GENES, Step, safe_merge)


# ---------------------------------------------------------------- Cellosaurus
def parse_cellosaurus(path):
    recs = []
    cur = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            code = line[:2]
            if code == "ID":
                cur = {"cvcl_id_name": line[5:].strip(), "AC": None, "AS": [], "SX": None, "DT": []}
            elif cur is None:
                continue
            elif code == "AC":
                cur["AC"] = line[5:].strip()
            elif code == "AS":
                cur["AS"] = [a.strip() for a in line[5:].strip().rstrip(";").split(";") if a.strip()]
            elif code == "SX":
                cur["SX"] = line[5:].strip()
            elif line.startswith("CC   Doubling time:"):
                cur["DT"].append(line[len("CC   Doubling time:"):].strip())
            elif code == "//":
                recs.append(cur)
                cur = None
    return recs


DT_RE = re.compile(r"(\d+(?:\.\d+)?)(?:\s*-\s*(\d+(?:\.\d+)?))?(?:\s*\+-\s*\d+(?:\.\d+)?)?\s*(hours?|days?)\b",
                   re.IGNORECASE)


def parse_doubling(texts):
    """Return (median hours, n values) from Cellosaurus doubling time comments."""
    vals = []
    for t in texts:
        t = re.sub(r"\([^()]*\)", " ", t)  # drop source parentheses
        for m in DT_RE.finditer(t):
            a = float(m.group(1))
            b = float(m.group(2)) if m.group(2) else None
            v = (a + b) / 2 if b is not None else a
            if m.group(3).lower().startswith("day"):
                v *= 24
            vals.append(v)
    if not vals:
        return np.nan, 0
    return float(np.median(vals)), len(vals)


# ---------------------------------------------------------------- cut-offs
def em_gmm2(x, seed=1, iters=500):
    rng = np.random.default_rng(seed)
    mu = np.sort(rng.choice(x, 2, replace=False)).astype(float)
    sd = np.array([x.std(), x.std()])
    w = np.array([0.5, 0.5])
    for _ in range(iters):
        p = np.stack([w[k] * np.exp(-0.5 * ((x - mu[k]) / sd[k]) ** 2) / sd[k] for k in range(2)])
        r = p / p.sum(0)
        w = r.mean(1)
        mu = (r * x).sum(1) / r.sum(1)
        sd = np.sqrt((r * (x - mu[:, None]) ** 2).sum(1) / r.sum(1)) + 1e-6
    grid = np.linspace(mu.min(), mu.max(), 2000)
    p = np.stack([w[k] * np.exp(-0.5 * ((grid - mu[k]) / sd[k]) ** 2) / sd[k] for k in range(2)])
    j = np.argmin(np.abs(p[0] - p[1]))
    return grid[j]


def bimodal_cutoff(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    kde = gaussian_kde(x)
    grid = np.linspace(x.min(), x.max(), 2000)
    d = kde(grid)
    # local maxima on the grid; an end point counts if it is above its single neighbour
    peaks = [i for i in range(len(d))
             if (i == 0 or d[i] > d[i - 1]) and (i == len(d) - 1 or d[i] >= d[i + 1])]
    if len(peaks) >= 2:
        i0, i1 = sorted(sorted(peaks, key=lambda i: d[i])[-2:])
        j = i0 + int(np.argmin(d[i0:i1 + 1]))
        return float(grid[j]), "kde_antimode", float(grid[i0]), float(grid[i1]), len(peaks)
    return float(em_gmm2(x)), "gmm_equal_posterior", np.nan, np.nan, len(peaks)


def expression_scores():
    """Per-model XIST and Y score from CMP (source 1) and CCLE legacy (source 2)."""
    s = pd.read_csv(os.path.join(RES, "C_expr_sanger_cmp_long.tsv.gz"), sep="\t")
    # plan 3b: drop duplicate == 't' rows where a non-duplicate row exists for (model, gene)
    has_f = s[s.duplicate == "f"].groupby(["model_id", "gene_symbol"]).size().rename("n_f").reset_index()
    s = s.merge(has_f, on=["model_id", "gene_symbol"], how="left")
    s = s[(s.duplicate == "f") | s.n_f.isna()]
    s["lx"] = np.log2(s.rsem_tpm + 1)
    w1 = s.groupby(["model_id", "gene_symbol"]).lx.mean().unstack()
    assert w1.index.is_unique
    y1 = [g for g in Y_GENES if g in w1.columns]
    e1 = pd.DataFrame({"xist": w1["XIST"], "yscore": w1[y1].mean(axis=1), "n_y_genes": w1[y1].notna().sum(1)})
    e1["expr_source"] = "CMP_rnaseq_latest_log2TPM1"

    c = pd.read_csv(os.path.join(RES, "C_expr_ccle_legacy_long.tsv.gz"), sep="\t")
    c["lx"] = np.log2(c.rpkm + 1)
    w2 = c.groupby(["ccle_name", "gene_symbol"]).lx.mean().unstack()
    y2 = [g for g in Y_GENES if g in w2.columns]
    e2 = pd.DataFrame({"xist": w2["XIST"], "yscore": w2[y2].mean(axis=1), "n_y_genes": w2[y2].notna().sum(1)})
    e2["expr_source"] = "CCLE_legacy_log2RPKM1"
    return e1, e2, y1, y2


def call_state(yhigh, xhigh, annot):
    if pd.isna(yhigh) or pd.isna(xhigh):
        return "not_measured"
    if yhigh and not xhigh:
        return "XY"
    if xhigh and not yhigh:
        return "XX-like"
    if yhigh and xhigh:
        return "ambiguous"
    if annot == "Male":
        return "Y-lost"
    if annot == "Female":
        return "XIST-lost"
    return "ambiguous"


def main():
    with Step("labels: load CMP model list"):
        m = pd.read_csv(os.path.join(DATA, "model_list_20240110.csv"), low_memory=False)
        assert m.model_id.is_unique
        m = m[["model_id", "model_name", "tissue", "cancer_type", "gender", "BROAD_ID", "CCLE_ID",
               "COSMIC_ID", "RRID", "model_type", "growth_properties"]].copy()
        assert set(m.gender.dropna()) <= {"Male", "Female", "Unknown"}
        m["gender"] = m.gender.fillna("Unknown")

    with Step("labels: parse Cellosaurus"):
        recs = parse_cellosaurus(os.path.join(DATA, "cellosaurus.txt"))
        cel = pd.DataFrame(recs)
        assert cel.AC.is_unique
        dtp = cel.DT.apply(parse_doubling)
        cel["doubling_time_h"] = [v for v, _ in dtp]
        cel["n_doubling_values"] = [n for _, n in dtp]
        cel["has_doubling_comment"] = cel.DT.apply(len) > 0
        by_ac = cel.set_index("AC")[["SX", "doubling_time_h", "n_doubling_values", "has_doubling_comment"]]
        sec = cel[["AC", "AS"]].explode("AS").dropna()
        sec = sec[~sec.AS.isin(by_ac.index)]
        sec = sec.drop_duplicates("AS", keep=False)  # ambiguous secondary accessions are not used

    with Step("labels: annotated sex"):
        mm = m[m.RRID.notna()].copy()
        mm = safe_merge(mm, by_ac.reset_index().rename(columns={"AC": "RRID"}), on="RRID", how="left",
                        validate="m:1", name="CMP RRID -> Cellosaurus AC")
        mm["cello_match"] = np.where(mm.SX.notna() | mm.has_doubling_comment.notna(), "AC", None)
        # secondary accessions for unmatched
        um = mm.has_doubling_comment.isna()
        if um.any():
            alt = sec.set_index("AS").AC
            mm.loc[um, "RRID_primary"] = mm.loc[um, "RRID"].map(alt)
            got = mm.RRID_primary.notna() & um
            for col in ["SX", "doubling_time_h", "n_doubling_values", "has_doubling_comment"]:
                mm.loc[got, col] = mm.loc[got, "RRID_primary"].map(by_ac[col])
            mm.loc[got, "cello_match"] = "AS"
        m = safe_merge(m, mm[["model_id", "SX", "doubling_time_h", "n_doubling_values", "cello_match"]],
                       on="model_id", how="left", validate="1:1", name="CMP -> Cellosaurus fields")
        m["cellosaurus_sex"] = m.SX.where(m.SX.isin(["Male", "Female"]), "Unknown")
        m["cellosaurus_sx_raw"] = np.where(m.SX.notna(), m.SX,
                                           np.where(m.cello_match.notna(), "no_SX_line", "no_entry"))

        def consensus(a, b):
            if a in ("Male", "Female") and b in ("Male", "Female"):
                return a if a == b else "conflict"
            if a in ("Male", "Female"):
                return a
            if b in ("Male", "Female"):
                return b
            return "Unknown"
        m["annotated_sex"] = [consensus(a, b) for a, b in zip(m.gender, m.cellosaurus_sex)]
        assert set(m.annotated_sex) <= {"Male", "Female", "conflict", "Unknown"}

    with Step("labels: expression states"):
        e1, e2, y1, y2 = expression_scores()
        cut_rows = []
        cuts = {}
        for src, e in [("CMP", e1), ("CCLE", e2)]:
            for var in ["xist", "yscore"]:
                cval, method, p0, p1, npk = bimodal_cutoff(e[var].values)
                cuts[(src, var)] = cval
                cut_rows.append({"expr_source": src, "variable": var, "cutoff_log2": cval, "method": method,
                                 "lower_peak": p0, "upper_peak": p1, "n_local_maxima": npk,
                                 "n_models": int(e[var].notna().sum()),
                                 "n_above": int((e[var] > cval).sum()),
                                 "y_genes_used": ",".join(y1 if src == "CMP" else y2) if var == "yscore" else ""})
        pd.DataFrame(cut_rows).to_csv(os.path.join(RES, "C_expression_cutoffs.tsv"), sep="\t", index=False)
        e1 = e1.add_prefix("cmp_")
        e1.index.name = "model_id"
        m = safe_merge(m, e1.reset_index(), on="model_id", how="left", validate="1:1", name="CMP expr")
        e2 = e2.add_prefix("ccle_")
        e2.index.name = "CCLE_ID"
        mc = m[m.CCLE_ID.notna()][["model_id", "CCLE_ID"]]
        assert mc.CCLE_ID.is_unique
        mc = safe_merge(mc, e2.reset_index(), on="CCLE_ID", how="left", validate="1:1", name="CCLE expr")
        m = safe_merge(m, mc.drop(columns="CCLE_ID"), on="model_id", how="left", validate="1:1",
                       name="CCLE expr -> models")
        for src, pre in [("CMP", "cmp_"), ("CCLE", "ccle_")]:
            yh = (m[pre + "yscore"] > cuts[(src, "yscore")]).where(m[pre + "yscore"].notna())
            xh = (m[pre + "xist"] > cuts[(src, "xist")]).where(m[pre + "xist"].notna())
            m[pre + "state"] = [call_state(a, b, s) for a, b, s in zip(yh, xh, m.annotated_sex)]
        m["expressed_state"] = np.where(m.cmp_state != "not_measured", m.cmp_state, m.ccle_state)
        m["expressed_state_source"] = np.where(m.cmp_state != "not_measured", "CMP",
                                               np.where(m.ccle_state != "not_measured", "CCLE", "none"))
        allowed = {"XY", "XX-like", "Y-lost", "XIST-lost", "ambiguous", "not_measured"}
        assert set(m.expressed_state) <= allowed

    with Step("labels: lineage and screen membership"):
        m["sex_specific_lineage"] = m.tissue.isin(SEX_SPECIFIC_TISSUES) | m.cancer_type.isin(SEX_SPECIFIC_CANCER_TYPES)
        m["lineage"] = m.cancer_type.where(~m.cancer_type.isin(DROP_CANCER_TYPES) & m.cancer_type.notna())
        g = pd.read_excel(os.path.join(DATA, "GDSC2_fitted_dose_response_27Oct23.xlsx"),
                          usecols=["SANGER_MODEL_ID"]).SANGER_MODEL_ID.dropna().unique()
        p = pd.read_csv(os.path.join(DATA, "prism_secondary_dose_response_curve_parameters.csv"),
                        usecols=["depmap_id"]).depmap_id.dropna().unique()
        m["in_gdsc2"] = m.model_id.isin(g)
        m["in_prism"] = m.BROAD_ID.isin(p) & m.BROAD_ID.notna()
        assert m.in_gdsc2.sum() == len(set(g)), "every GDSC2 model must be in CMP"
        m["log2_doubling_time_h"] = np.log2(m.doubling_time_h)

    with Step("labels: write outputs"):
        out = m[m.in_gdsc2 | m.in_prism].copy()
        cols = ["model_id", "model_name", "BROAD_ID", "CCLE_ID", "COSMIC_ID", "RRID", "in_gdsc2", "in_prism",
                "tissue", "cancer_type", "lineage", "sex_specific_lineage", "gender", "cellosaurus_sx_raw",
                "cellosaurus_sex", "cello_match", "annotated_sex",
                "cmp_xist", "cmp_yscore", "cmp_n_y_genes", "cmp_state", "ccle_xist", "ccle_yscore",
                "ccle_n_y_genes", "ccle_state", "expressed_state", "expressed_state_source",
                "doubling_time_h", "n_doubling_values", "log2_doubling_time_h", "growth_properties"]
        out[cols].to_csv(os.path.join(RES, "C_cellline_sex_labels.tsv"), sep="\t", index=False)
        ss = (m[m.sex_specific_lineage & (m.in_gdsc2 | m.in_prism)]
              .groupby(["tissue", "cancer_type"]).agg(n_models=("model_id", "size"),
                                                     n_in_gdsc2=("in_gdsc2", "sum"),
                                                     n_in_prism=("in_prism", "sum")).reset_index())
        ss["rule"] = np.where(ss.tissue.isin(SEX_SPECIFIC_TISSUES), "tissue", "cancer_type")
        ss.to_csv(os.path.join(RES, "C_sex_specific_lineages.tsv"), sep="\t", index=False)
        disc = out[(out.gender != out.cellosaurus_sex) & out.gender.isin(["Male", "Female"])
                   & out.cellosaurus_sex.isin(["Male", "Female"])]
        disc[["model_id", "model_name", "RRID", "gender", "cellosaurus_sx_raw", "expressed_state",
              "in_gdsc2", "in_prism"]].to_csv(os.path.join(RES, "C_annotation_discordance.tsv"), sep="\t",
                                              index=False)
        # summary counts used in the report
        rows = []
        def add(table, sub, a, b):
            for (i, j), n in pd.crosstab(sub[a], sub[b], margins=True, margins_name="total").stack().items():
                rows.append({"table": table, "row": i, "col": j, "n": int(n)})
        prim = out[out.annotated_sex.isin(["Male", "Female"]) & out.lineage.notna() & ~out.sex_specific_lineage]
        add("screened_lines: CMP gender x Cellosaurus SX", out, "gender", "cellosaurus_sx_raw")
        add("screened_lines: annotated_sex x expressed_state", out, "annotated_sex", "expressed_state")
        add("primary_eligible_lines: annotated_sex x expressed_state", prim, "annotated_sex", "expressed_state")
        for scr in ["in_gdsc2", "in_prism"]:
            add(f"primary_eligible_lines {scr}: annotated_sex x expressed_state", prim[prim[scr]],
                "annotated_sex", "expressed_state")
        both = out[(out.cmp_state != "not_measured") & (out.ccle_state != "not_measured")]
        add("lines with both expression sources: CMP state x CCLE state", both, "cmp_state", "ccle_state")
        rows.append({"table": "lines with both expression sources: CMP state equals CCLE state", "row": "agree",
                     "col": "n", "n": int((both.cmp_state == both.ccle_state).sum())})
        add("screened_lines: expressed_state_source", out.assign(all="all"), "expressed_state_source", "all")
        add("primary_eligible_lines: has doubling time x screen", prim.assign(
            has_dt=prim.doubling_time_h.notna(), screen=np.where(prim.in_gdsc2 & prim.in_prism, "both",
                                                                np.where(prim.in_gdsc2, "GDSC2_only", "PRISM_only"))),
            "has_dt", "screen")
        pd.DataFrame(rows).to_csv(os.path.join(RES, "C_label_summary.tsv"), sep="\t", index=False)
        print(out.annotated_sex.value_counts().to_dict())
        print(pd.crosstab(out.annotated_sex, out.expressed_state))


if __name__ == "__main__":
    main()
