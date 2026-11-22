"""Analysis 7 (../PLAN.md, deviation 7): E3 on the real labels with the response size S built without the weak-template
compounds, and from the compounds above the panel median activity only (post hoc; joins the real states to responses, as
the frozen run did).

Weak templates: the bottom quartile of the P1 compounds by template norm (|mean response over the available P1 lines|)
on each scale, the label-free classes of analysis 2 (common.template_classes). Activity: the frozen label-free
classification of task D (taskD/D_activity_scores.tsv, P1, above_panel_median), used by the frozen run for its activity
subset of H1. E3 is scored with ganalyses.e3_row on the frozen designs (rule 7i A). The all-compound rows must
reproduce results/G/G_H3.tsv. Writes results/e3_weak_templates.tsv.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as K  # noqa: E402

import numpy as np
import pandas as pd

GS, GA = K.GS, K.GA


def main():
    D = K.Data(with_level3=True)
    P, cfg = D.P, D.cfg
    lt = pd.read_csv(D.I.real_line_table_path(), sep="\t")
    al = GS.Aligned(P, lt)
    prim = GS.primary_settings(cfg)
    s5, s3 = D.s0()
    _, weak5, _ = K.template_classes(P.R)
    _, weak3, _ = K.template_classes(D.R3)
    act = D.I.activity()
    act = act[act.panel == "P1"].set_index("cmap_name").above_panel_median.astype(bool)
    above = np.array([j for j, c in enumerate(P.comps) if bool(act.get(c, False))])
    allc = np.arange(len(P.comps))
    h3 = pd.read_csv(os.path.join(K.GRES, "G_H3.tsv"), sep="\t", float_precision="round_trip")
    rows = []
    for term in ("factY", "factX", "annotated"):
        d = al.design(term, prim)
        for scale, s, weak, tn in (("level3", s3, weak3, "R3_P1_d111"), ("level5", s5, weak5, "P1_resp")):
            keep_w = np.setdiff1d(allc, weak)
            for name, cols in (("all compounds (primary)", allc), (f"weak-template compounds left out ({len(weak)} of {len(allc)})", keep_w),
                               (f"compounds above the panel median activity only ({len(above)} of {len(allc)})", above)):
                e = GA.e3_row(D.ctx, "P1", term, d, K.mean_over(s, cols), scale, tn, "weak-template sensitivity")
                r = dict(term=term, scale=scale, variant=name, n_compounds=len(cols), n_lines=e["n_lines"], df=e["df"], mean_S=e["mean_S"],
                         estimate=e["estimate"], lo90=e["lo90"], hi90=e["hi90"], lo95=e["lo95"], hi95=e["hi95"],
                         difference_declared=e["difference_declared"])
                if name.startswith("all"):
                    ref = h3[(h3.panel == "P1") & (h3.term == term) & (h3.scale == scale) & (h3.role.isin(["primary", "S12"]))].iloc[0]
                    r["reproduces_G_H3"] = bool(all(abs(e[c] - ref[c]) <= 1e-9 * max(1.0, abs(ref[c])) for c in ("estimate", "lo95", "hi95", "mean_S")))
                    assert r["reproduces_G_H3"], (term, scale)
                rows.append(r)
    out = pd.DataFrame(rows)
    out["weak_template_definition"] = "bottom quartile of the P1 compounds by |mean response over the available P1 lines|, on the same scale"
    out["activity_source"] = "taskD/D_activity_scores.tsv, panel P1, column above_panel_median"
    out.to_csv(os.path.join(K.RES, "e3_weak_templates.tsv"), sep="\t", index=False)
    print(out[["term", "scale", "variant", "n_compounds", "n_lines", "estimate", "lo95", "hi95", "difference_declared"]].to_string(index=False))


if __name__ == "__main__":
    main()
