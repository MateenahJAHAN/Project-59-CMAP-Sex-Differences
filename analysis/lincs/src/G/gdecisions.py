"""Decision table of task G (draft section 6): each pre-registered item, its rule and its outcome, with the
wording fixed before unblinding. Output: results/G/G_decisions.tsv."""
import numpy as np, pandas as pd

NOWORD = "Wording in every case: 'lines that express Y genes', 'lines that express XIST', XY, XX-like, Y-lost, XIST-lost; never 'female cells', 'male cells' or 'sex does not matter'; no causal or patient claim."


def build(cfg, gates, h1, h2, h3, e45, rep_panel, leads, sec, fails_stop):
    a1 = cfg.const("H1_alpha_per_term")
    rows = []
    add = lambda **k: rows.append(k)
    pc1 = gates[gates.gate == "PC1"].iloc[0] if (gates.gate == "PC1").any() else None
    pc1_pass = bool(pc1 is not None and pc1.passes)
    add(section="5a", item="PC1 (K2)", rule="AMG-232 p < 0.05 and rank in the top 5% (lineage + TP53, CMP), level 5",
        outcome="pass" if pc1_pass else "fail", consequence="" if pc1_pass else "stop and debug; if it still fails, P1 is uninformative and no LINCS state claim is made")
    k4 = gates[gates.gate == "K4"]
    k4_pass = bool(len(k4) and k4.passes.iloc[0])
    add(section="5a, 6 K4", item="K4 (level 3)", rule="PC1 passes on the level-3 tensor of 7d and CDKN1A keeps its level-5 sign",
        outcome="pass" if k4_pass else ("fail" if len(k4) else "not run"),
        consequence="E3 on level 3" if k4_pass else "E3 on level 5, labelled biased toward zero, no equivalence claim")
    for _, g in gates[gates.gate == "calibration"].iterrows():
        add(section="5b, 6 K6", item=f"calibration {g.panel} {g.term}", rule="mean share of compounds with p < 0.05 in [0.03, 0.07]",
            outcome="pass" if g.passes else "fail", consequence="" if g.passes else "stop and fix the code; do not unblind")
    dm = gates[gates.gate == "DMSO"]
    if len(dm):
        add(section="5c, 7b", item="DMSO gate (7b C)", rule="P1 DMSO pseudo-compound p >= 0.05 under the PC1 model",
            outcome="pass" if dm.passes.iloc[0] else "fail", consequence="" if dm.passes.iloc[0] else "Part B does not run")
    opt7b = cfg.decide("7b")
    for _, r in h1[(h1.panel == cfg.const("primary_panel")) & h1.term.isin(["factY", "factX"]) & (h1.analysis == "primary")].iterrows():
        if not r.contrast_ok or r.n_tested == 0:
            add(section="6 K3", item=f"H1 {r.term}", rule="testability (K3)", outcome="not analysed", consequence=str(r.contrast_reason)); continue
        rej = bool(r.panel_p < a1)
        words = []
        if rej:
            words.append("panel-level structure in P1")
            if opt7b == "B" and np.isfinite(r.dmso_p) and r.dmso_p < 0.05:
                words.append("described as a shift in average response that also shows in vehicle z-scores; only the H2 specificity rule can support a compound claim (7b B)")
            e = e45[(e45.panel == r.panel) & (e45.term == r.term) & (e45.comparison == "ploidy within-lineage median split")]
            if len(e) and np.isfinite(e.E5_comparison.iloc[0]) and e.E5_comparison.iloc[0] >= r.E5_excess:
                words.append("not separable from ploidy")
            rp = rep_panel[rep_panel.term == r.term]
            if len(rp):
                words.append(rp.label.iloc[0])
        else:
            if not pc1_pass:
                words.append("PC1 failed: P1 uninformative")
            else:
                t = h2[(h2.panel == r.panel) & (h2.term == r.term) & (h2.analysis == "primary") & h2.in_ranking]
                e1 = float(np.nanmedian(t.E1_hw90_hc3_sd_median)) if len(t) else np.nan
                e2 = float(np.nanmedian((t.E2_hi90 - t.E2_lo90) / 2)) if len(t) else np.nan
                e3 = h3[(h3.panel == r.panel) & (h3.term == r.term) & (h3.role == "primary")]
                e3s = (f"E3 {e3.estimate.iloc[0]:.1f} (90% CI {e3.lo90.iloc[0]:.1f} to {e3.hi90.iloc[0]:.1f}; {e3.unit.iloc[0]})"
                       if len(e3) and "estimate" in e3 and np.isfinite(e3.estimate.iloc[0]) else "E3 not estimable")
                who = "lines that express Y genes" if r.term == "factY" else "lines that express XIST"
                words.append(f"null wording: within lineage, at 1.0 uM and 24 h, {who} respond to the same compounds as lines that do not, "
                             f"within the stated bounds: median per-compound E1 90% half-width {e1:.2f} SD per gene, E2 90% half-width "
                             f"{e2:.1f} ({t.E2_unit.iloc[0] if len(t) else ''}); {e3s}. A bound, not proof of no effect, in cancer cell lines only.")
        add(section="4, 6", item=f"H1 {r.term} (P1)", rule=f"panel statistic permutation p < {a1}", outcome="rejected" if rej else "not rejected",
            consequence="; ".join(words))
    for term in ("factY", "factX"):
        t = h2[(h2.panel == cfg.const("primary_panel")) & (h2.term == term) & (h2.analysis == "primary")]
        if not len(t):
            continue
        nf = int(t.finding.sum()); nold = int((t.meets_rules_1_to_4 & t.old_list).sum())
        e4 = bool(t.rule3_E4.fillna(False).astype(bool).any())
        notes = [f"{nold} old-list compounds meet rules 1 to 4 (not findings from P1 alone)" if nold else "",
                 "findings are candidates for follow-up, not claims about patients" if nf else "",
                 "" if e4 else "the discovery count is inside the joint null's 95th percentile (E4): compatible with chance"]
        add(section="4, 6", item=f"H2 {term} (P1)", rule="finding only if rules 1 to 4 hold; old-list compounds need replication",
            outcome=f"{nf} findings", consequence="; ".join(x for x in notes if x))
    for _, r in h3[(h3.panel == cfg.const("primary_panel")) & (h3.role == "primary")].iterrows():
        va = h3[(h3.panel == r.panel) & (h3.term == r.term) & (h3.role == "viability-adjusted")]
        gone = bool(len(va) and "difference_disappears_after_viability_adjustment" in va
                    and va.difference_disappears_after_viability_adjustment.iloc[0] == True)
        add(section="4, 6 (H3)" if r.term != "annotated" else "4 (H4)", item=f"E3 {r.term} ({r.scale})",
            rule="difference if the 95% CI excludes 0; equivalence only if the label-free probability is at least 0.8",
            outcome=("difference declared" if r.get("difference_declared", False) == True else "no difference declared"),
            consequence=str(r.get("equivalence", "")) + ("; the difference disappears after the PRISM adjustment: a difference in "
                                                         "sensitivity, not a transcriptional program" if gone else ""))
    add(section="4 (H4)", item="annotated sex", rule="reported beside H1 to H3; no claim of its own", outcome="reported", consequence="see G_H4.tsv")
    for _, r in leads.iterrows():
        if r.panel != cfg.const("primary_panel"):
            continue
        add(section="1, 4 (leads)", item=f"{r.lead} {r.compound} ({r.term}, {r.alternative})", rule="alpha 0.05; reported with and without Holm",
            outcome=("p < 0.05" if r.get("p", np.nan) < 0.05 else "p >= 0.05"), consequence=f"Holm-adjusted p {r.get('p_holm', np.nan)}")
    add(section="1 (secondary)", item="secondary family", rule="BH within each readout; labelled secondary", outcome="reported", consequence="see G_secondary_*.tsv")
    add(section="6", item="P1b and P2", rule="reported without claims of their own, except the P2 label of 7a", outcome="reported", consequence="")
    add(section="6", item="wording", rule=NOWORD, outcome="applies", consequence="")
    if fails_stop:
        add(section="5", item="run status", rule="gates must pass before any result is read", outcome="stopped", consequence="; ".join(fails_stop))
    return pd.DataFrame(rows)
