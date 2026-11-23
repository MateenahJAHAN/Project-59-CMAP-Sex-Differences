"""Analysis 4 (../PLAN.md): sensitivity of E3 to line misclassification (post hoc; joins the real states to responses,
as the frozen run did).

Flagged lines: revision_v2/analyses/partA_dna_transfer/partB_lines_transfer.tsv (column flag = yes) when it exists,
else the line table (task A confidence not high, source conflict or discordant). Variants, for the factorial Y and X
terms on level 3 (primary scale) and level 5: (a) all flagged lines excluded; (b) each flagged line excluded alone;
(c) each flagged line reassigned (a Y+ or X+ line to neither; a neither line to its alternative group, if any);
(d) PLAN.md deviation 4: A549 as Y+ and X+. The analysis set is rebuilt with the frozen rule (gcore.analysis_set), and
E3 is scored with ganalyses.e3_row (rule 7i A). For variant (a) the panel-level test is also run (B = 10,000).
Writes results/misclass_e3.tsv and results/misclass_flags.tsv.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as K  # noqa: E402
import time

import numpy as np
import pandas as pd

C, GS, GA = K.C, K.GS, K.GA
TRANSFER = os.environ.get("P59_TRANSFER_TABLE", "")


def flagged(lt, f41):
    if "P59_LINCS_WORK" in os.environ:
        measured = os.path.join(os.environ["P59_LINCS_WORK"], "results", "transfer", "sensitivity_lines.tsv")
        source = "results/transfer/sensitivity_lines.tsv" if os.path.isfile(measured) else "sensitivity_lines.tsv"
        table = measured if os.path.isfile(measured) else os.path.join(K.HERE.parent, source)
        flagged_lines = pd.read_csv(table, sep="\t")
        if len(f41) != 41 or not set(flagged_lines.cell_iname).issubset(f41):
            raise ValueError("Sensitivity lines do not match the primary panel")
        return flagged_lines.assign(flag_source=source, flag_source_sha256=K.sha256_file(table))
    if os.path.exists(TRANSFER):
        t = pd.read_csv(TRANSFER, sep="\t")
        assert set(t.line) == f41, "transfer table lines differ from the 41-line set"
        f = t[t.flag == "yes"][["line", "group", "flag_reasons", "alt_groups"]].rename(columns={"line": "cell_iname"})
        return f.assign(flag_source=TRANSFER, flag_source_sha256=K.sha256_file(TRANSFER))
    s = lt[lt.cell_iname.isin(f41)]
    m = (s.taskA_confidence != "high") | s.source_conflict.fillna(False).astype(bool) | (s.taskA_state == "discordant")
    f = s[m][["cell_iname"]].assign(group="", flag_reasons="line table: confidence, conflict or discordant", alt_groups=np.nan)
    return f.assign(flag_source="taskD/D_line_table.tsv", flag_source_sha256="")


def main():
    t_all = time.time()
    D = K.Data(with_level3=True)
    P, cfg = D.P, D.cfg
    lt = pd.read_csv(D.I.real_line_table_path(), sep="\t")
    al = GS.Aligned(P, lt)
    prim = GS.primary_settings(cfg)
    sl = pd.read_csv(os.path.join(K.GRES, "G_sets_lines.tsv"), sep="\t")
    f41 = set(sl[(sl.panel == "P1") & (sl.family == "factorial")].cell_iname)
    fl = flagged(lt, f41)
    fl.to_csv(os.path.join(K.RES, "misclass_flags.tsv"), sep="\t", index=False)
    s5, s3 = D.s0()
    S = {"level3": np.nanmean(s3, axis=1), "level5": np.nanmean(s5, axis=1)}
    pos = {l: i for i, l in enumerate(P.lines)}
    flines = list(fl.cell_iname)
    rows = []
    for term in K.TERMS:
        d = al.design(term, prim)
        iy = d["tcols"][0]
        variants = [("primary (41 lines)", d["Xo"], d["base"], "")]
        excl = d["base"].copy()
        for l in flines:
            excl[pos[l]] = False
        variants.append((f"all flagged lines excluded ({len(flines)})", d["Xo"], excl, ";".join(flines)))
        for l in flines:
            b = d["base"].copy(); b[pos[l]] = False
            variants.append((f"{l} excluded", d["Xo"], b, l))
        for l in flines:
            Xo = d["Xo"].copy(); i = pos[l]
            y, x = Xo[i, iy], Xo[i, iy + 1]
            if y == 1 or x == 1:
                Xo[i, iy], Xo[i, iy + 1] = 0.0, 0.0
                variants.append((f"{l} reassigned to neither", Xo, d["base"], l))
        if "A549" in flines:
            Xo = d["Xo"].copy(); i = pos["A549"]
            Xo[i, iy], Xo[i, iy + 1] = 1.0, 1.0
            variants.append(("A549 reassigned to Y+ and X+ (deviation 4)", Xo, d["base"], "A549"))
        for name, Xo, base, changed in variants:
            mask = C.analysis_set(d["lin"], Xo, d["tcols"], base)
            dd = dict(d, Xo=Xo, mask=mask, base=base)
            cc = al.contrast_check(dd)
            for scale, tn in (("level3", "R3_P1_d111"), ("level5", "P1_resp")):
                e = GA.e3_row(D.ctx, "P1", term, dd, S[scale], scale, tn, "misclassification sensitivity")
                rows.append(dict(term=term, variant=name, lines_changed=changed, scale=scale, n_lines=e.get("n_lines"),
                                 n_g1=cc["n_g1"], n_g0=cc["n_g0"], n_lineages=cc["n_lineages"], df=e.get("df"), mean_S=e.get("mean_S"),
                                 estimate=e.get("estimate"), lo90=e.get("lo90"), hi90=e.get("hi90"), lo95=e.get("lo95"), hi95=e.get("hi95"),
                                 difference_declared=e.get("difference_declared"), contrast_ok=cc["contrast_ok"]))
            if name.startswith("all flagged"):
                df, permT, recs = C.run_term(P, Xo, d["k"], mask, cfg.run("B_full"), rule=d["rule"])
                dfF, obs, pp, _ = C.finish_table(df, permT, compounds_only=P.comps)
                t = dfF[dfF.in_ranking]
                rows.append(dict(term=term, variant=name, lines_changed=changed, scale="level5 (panel-level test)",
                                 n_lines=int(mask.sum()), n_g1=cc["n_g1"], n_g0=cc["n_g0"], n_lineages=cc["n_lineages"],
                                 panel_stat=obs, panel_p=pp, n_tested=len(t), n_p_lt_005=int((t.p < 0.05).sum()),
                                 n_q_lt_01=int((t.q < 0.1).sum()), B=cfg.run("B_full"), contrast_ok=cc["contrast_ok"]))
    out = pd.DataFrame(rows)
    out["flag_source"] = fl.flag_source.iloc[0] if len(fl) else ""
    out.to_csv(os.path.join(K.RES, "misclass_e3.tsv"), sep="\t", index=False)
    print(fl.to_string(index=False))
    print(out[["term", "variant", "scale", "n_lines", "n_g1", "n_g0", "estimate", "lo95", "hi95", "panel_p"]].to_string(index=False))
    print(f"done in {time.time() - t_all:.0f} s")


if __name__ == "__main__":
    main()
