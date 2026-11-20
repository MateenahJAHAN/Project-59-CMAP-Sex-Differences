"""Exposure designs and line sets for task G.

A design is everything one run of the engine needs for one tested term: Xo (lines x columns; covariates first,
then the tested terms), the tested column k, the tested columns of the model (tcols), the analysis-set mask,
the permutation blocks and the 7i rule.

DECIDE choices applied here:
  7f  which lines enter the factorial contrast;
  7g  TP53 coding;
  7i  minimum group size and distinct relabellings (per compound in gcore.compound_design, per contrast here);
  7c  whether the strict term exists.
Sensitivity variants (S1 to S10) are expressed as design settings. Only this module and the merge steps of
ganalyses.py read state or label columns.
"""
import numpy as np, pandas as pd
import gcore as C

TASKC_YX = {"XY": (1.0, 0.0), "XX-like": (0.0, 1.0), "Y-lost": (0.0, 0.0), "XIST-lost": (0.0, 0.0)}
LT_COLS = ["cell_iname", "lineage", "taskA_state", "taskC_state", "source_conflict", "elig_factorial", "fact_Y", "fact_X",
           "elig_strict", "strict_XY", "elig_annotated", "annot_M", "cmap_donor_sex"]


def rule_for(cfg):
    if cfg.decide("7i") == "A":
        return dict(min_group=cfg.const("min_group_7iA"), min_relab=cfg.const("min_relabellings_7iA"))
    return dict(min_group=cfg.const("min_group_7iB"), min_relab=0)


def tp53_code(cfg):
    return "7gA" if cfg.decide("7g") == "A" else "7gB"


def tp53_columns(cov, code):
    """TP53 covariate columns from the state-free covariate table (7g). 7gA: CMP call with three levels (mutant and
    unknown indicators); 7gB: task D's TP53_cov (CMP, else CCLE; NaN drops the line); ccle: CCLE MAF call (S10)."""
    if code == "none":
        return np.zeros((len(cov), 0)), []
    if code == "7gA":
        v = cov.TP53_cmp.astype(str).values
        assert set(v) <= {"mutant", "wild_type", "unknown"}
        return np.c_[(v == "mutant").astype(float), (v == "unknown").astype(float)], ["TP53_cmp_mutant", "TP53_cmp_unknown"]
    col = {"7gB": "TP53_cov", "ccle": "TP53_ccle", "cmp": "TP53_cmp"}[code]
    v = cov[col].astype(str).values
    x = np.where(v == "mutant", 1.0, np.where(v == "wild_type", 0.0, np.nan))
    return x[:, None], [col]


def tp53_strings(cov, code):
    """TP53 as strings for task D2's precision fit (mutant / wild_type / unknown)."""
    col = {"7gA": "TP53_cmp", "7gB": "TP53_cov", "ccle": "TP53_ccle"}[code]
    return cov[col].astype(str).values


def primary_settings(cfg):
    return dict(f=cfg.decide("7f"), tp53=tp53_code(cfg), rule=rule_for(cfg), block="lineage", y_only=False,
                extra_cov=(), line_variant=None, qc_drop=False, same_lines=None, label="primary")


class Aligned:
    """The line table and covariates aligned to the lines of panel P."""

    def __init__(self, P, lt):
        missing = [c for c in LT_COLS if c not in lt.columns]
        assert not missing, f"line table lacks {missing}"
        a = pd.DataFrame({"cell_iname": P.lines}).merge(lt[LT_COLS], on="cell_iname", how="left", validate="1:1")
        assert len(a) == len(P.lines) and (a.lineage.values == P.cov.lineage.values).all()
        for c in ["fact_Y", "fact_X", "strict_XY", "annot_M"]:
            assert set(a[c].dropna().unique()) <= {0.0, 1.0}, c
        for c in ["elig_factorial", "elig_strict", "elig_annotated", "source_conflict"]:
            assert set(a[c].dropna().unique()) <= {True, False}, c
        self.P, self.a, self.cov = P, a, P.cov

    # ---------- covariate columns ----------
    def tp53(self, code):
        return tp53_columns(self.cov, code)

    def extra(self, names):
        cols = []
        for nm in names:
            if nm == "log2_doubling_time":
                v = self.cov.doubling_time_h.to_numpy(dtype=float)
                cols.append(np.where(v > 0, np.log2(np.where(v > 0, v, 1.0)), np.nan))
            elif nm == "ploidy":
                cols.append(self.cov.ploidy.to_numpy(dtype=float))
            else:
                raise ValueError(nm)
        return (np.column_stack(cols) if cols else np.zeros((len(self.cov), 0))), list(names)

    def block(self, kind):
        if kind == "lineage":
            return self.P.lin
        if kind == "none":
            return np.array(["all"] * len(self.P.lines), dtype="U64")
        if kind == "cmp_cancer_type":
            return np.array([str(x) if pd.notna(x) else "__no_cancer_type__" for x in self.cov.cmp_cancer_type], dtype="U64")
        raise ValueError(kind)

    # ---------- exposures ----------
    def factorial_values(self, s):
        a = self.a
        Y, X = a.fact_Y.to_numpy(dtype=float).copy(), a.fact_X.to_numpy(dtype=float).copy()
        elig = a.elig_factorial.fillna(False).values.astype(bool)
        disc = (a.taskA_state == "discordant").values
        conf = a.source_conflict.fillna(False).values.astype(bool)
        tc = a.taskC_state.astype(str).values
        f = s["f"] if s["line_variant"] is None else {"D44": "B", "F41": "A", "taskC_conflict": "A"}[s["line_variant"]]
        if f == "B":
            base = elig.copy()
        elif f == "A":
            base = elig & ~disc & ~conf
        else:  # C: drop source-conflict lines whose Task C call changes Y or X
            chg = np.array([cf and (t in TASKC_YX) and (TASKC_YX[t] != (y, x)) for cf, t, y, x in zip(conf, tc, Y, X)])
            base = elig & ~chg
        if s["line_variant"] == "taskC_conflict":
            add = conf & np.isin(tc, list(TASKC_YX))
            for i in np.flatnonzero(add):
                Y[i], X[i] = TASKC_YX[tc[i]]
            base = base | add
        return Y, X, base

    def design(self, term, s):
        """term in factY, factX, annotated, strict. s: settings dict."""
        T, tnames = self.tp53(s["tp53"])
        Ex, enames = self.extra(s["extra_cov"])
        cov = np.c_[T, Ex] if Ex.shape[1] else T
        cnames = tnames + enames
        if term in ("factY", "factX"):
            Y, X, base = self.factorial_values(s)
            if s["y_only"]:
                assert term == "factY"
                Xo = np.c_[cov, Y]; k = cov.shape[1]; tcols = [k]; names = cnames + ["Y"]
            else:
                Xo = np.c_[cov, Y, X]; iy = cov.shape[1]
                k = iy if term == "factY" else iy + 1
                tcols = [iy, iy + 1]; names = cnames + ["Y", "X"]
        elif term == "annotated":
            base = self.a.elig_annotated.fillna(False).values.astype(bool)
            Xo = np.c_[cov, self.a.annot_M.to_numpy(dtype=float)]; k = cov.shape[1]; tcols = [k]; names = cnames + ["annot_M"]
        elif term == "strict":
            base = self.a.elig_strict.fillna(False).values.astype(bool)
            Xo = np.c_[cov, self.a.strict_XY.to_numpy(dtype=float)]; k = cov.shape[1]; tcols = [k]; names = cnames + ["XY"]
        else:
            raise ValueError(term)
        if s["qc_drop"]:
            q = self.cov[f"qc_pass_frac_{self.P.pn}"].to_numpy(dtype=float)
            base = base & (q >= s["qc_threshold"])
        lin = self.block(s["block"])
        if s["block"] == "cmp_cancer_type":
            base = base & (lin != "__no_cancer_type__")
        if s.get("same_lines") is not None:
            mask = s["same_lines"] & np.all(np.isfinite(Xo), axis=1)
            if s["block"] == "cmp_cancer_type":
                mask = C.analysis_set(lin, Xo, tcols, mask & base)
        else:
            mask = C.analysis_set(lin, Xo, tcols, base)
        d = dict(term=term, settings=s.get("label", ""), Xo=Xo, k=k, tcols=tcols, mask=mask, lin=lin, rule=s["rule"],
                 base=base, col_names=names)
        d.update(self.contrast_check(d))
        return d

    def contrast_check(self, d):
        """K3 at contrast level: groups and distinct relabellings over the whole analysis set (7i A), else task D's
        rule (each group at least 2 lines)."""
        m = d["mask"]
        t = d["Xo"][m, d["k"]]
        n1, n0 = int(np.sum(t == 1)), int(np.sum(t == 0))
        mg, mr = d["rule"]["min_group"], d["rule"]["min_relab"]
        nrel = C.n_relabellings(d["lin"][m], t) if m.sum() else 0
        ok, why = True, ""
        if min(n1, n0) < mg:
            ok, why = False, f"a group has fewer than {mg} lines"
        elif mr and nrel < mr:
            ok, why = False, f"fewer than {mr} distinct within-lineage relabellings"
        return dict(n_lines=int(m.sum()), n_g1=n1, n_g0=n0, n_lineages=len(set(d["lin"][m])), n_relabellings=nrel,
                    contrast_ok=ok, contrast_reason=why)

    # ---------- line sets and counts for the label-free gates ----------
    def set_and_counts(self, family, s):
        """Membership (no group column) and per-lineage group counts of the analysis set of a contrast family."""
        term = {"factorial": "factY", "annotated": "annotated", "strict": "strict"}[family]
        d = self.design(term, s)
        m = d["mask"]
        lines = [self.P.lines[i] for i in np.flatnonzero(m)]
        lin = self.P.lin
        rows = []
        for l in sorted(set(lin[m])):
            idx = m & (lin == l)
            if family == "factorial":
                Y = d["Xo"][idx, d["tcols"][0]]; X = d["Xo"][idx, d["tcols"][1]]
                g = [int(np.sum((Y == 1) & (X == 0))), int(np.sum((Y == 0) & (X == 1))), int(np.sum((Y == 0) & (X == 0))),
                     int(np.sum((Y == 1) & (X == 1)))]
            else:
                t = d["Xo"][idx, d["k"]]
                g = [int(np.sum(t == 1)), int(np.sum(t == 0)), 0, 0]
            rows.append(dict(lineage=l, n_lines=int(idx.sum()), n_g1=g[0], n_g2=g[1], n_g3=g[2], n_g4=g[3]))
        return lines, pd.DataFrame(rows), d
