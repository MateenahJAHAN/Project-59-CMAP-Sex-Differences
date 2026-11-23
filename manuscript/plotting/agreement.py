"""Join independent expression and STR measurements for Figure 3."""
import numpy as np
import pandas as pd
import layout as V
from fig3_agreement import CONFIDENT, PATTERNS, STR_LEVELS, STR_SOURCES, Y_LEVELS


def pattern(y, x, yc, xc):
    Y, X = y >= yc, x >= xc
    return {(True, False): "Y only", (False, True): "XIST only", (True, True): "both", (False, False): "neither"}[(Y, X)]


def from_sources(sources):
    def read(key):
        frame = V.read_tsv(sources[key])
        frame["ln"] = frame["_line"]
        return frame

    t = read("tahoe")
    tc = read("tahoe_cut")
    t_yc = float(V.one_row(tc, level="lineplate", score="y_score")["cutoff"])
    t_xc = float(V.one_row(tc, level="lineplate", score="XIST")["cutoff"])
    t["pattern"] = [pattern(y, x, t_yc, t_xc) for y, x in zip(t.y_score, t.XIST)]
    t["ycall"] = np.where(t.y_score >= t_yc, Y_LEVELS[0], Y_LEVELS[1])
    c = read("c_labels")
    cc = read("c_cut")
    ccut = {(r.expr_source, r.variable): (float(r.cutoff_log2), int(r.ln)) for r in cc.itertuples()}

    def c_pattern(r):
        s = r.expressed_state_source
        if s not in ("CMP", "CCLE"):
            return None
        x, y = (r.cmp_xist, r.cmp_yscore) if s == "CMP" else (r.ccle_xist, r.ccle_yscore)
        return pattern(y, x, ccut[(s, "yscore")][0], ccut[(s, "xist")][0])
    c["pattern"] = [c_pattern(r) for r in c.itertuples()]
    assert c.RRID.is_unique
    f1c = read("f1_cut")
    e_yc = float(V.one_row(f1c, dataset="EmeraldBay", score="y_score")["cutoff"])
    e_xc = float(V.one_row(f1c, dataset="EmeraldBay", score="XIST")["cutoff"])
    s_yc = float(V.one_row(f1c, dataset="sci-Plex (all samples)", score="y_score")["cutoff"])
    s_xc = float(V.one_row(f1c, dataset="sci-Plex (all samples)", score="XIST")["cutoff"])
    e = read("emerald")
    e["pattern"] = [pattern(y, x, e_yc, e_xc) for y, x in zip(e.y_score, e.XIST)]
    sp = read("sciplex")
    sp["pattern"] = [pattern(y, x, s_yc, s_xc) for y, x in zip(sp.y_score, sp.XIST)]
    st_all = V.read_tsv(sources["str"], dtype=str, keep_default_na=False)
    st_all["ln"] = st_all["_line"]
    st = st_all[st_all.in_tahoe == "True"]
    by_src = read("str_by_source").set_index("source")
    str_counts = {}
    for key, _ in STR_SOURCES:
        q = st_all[st_all[f"{key}_Y"].isin(["Y expressed", "Y below cut-off"]) &
                   st_all.amelogenin_class.isin(STR_LEVELS)]
        xy, yexp = q.amelogenin_class == "X,Y", q[f"{key}_Y"] == "Y expressed"
        recount = dict(n_compared=len(q), n_agree=int((xy == yexp).sum()), n_type_A=int((xy & ~yexp).sum()),
                       n_type_B=int((~xy & yexp).sum()))
        row = by_src.loc[key]
        assert all(int(row[k]) == v for k, v in recount.items()), (key, recount)
        assert (q[f"{key}_vs_STR"].map({"agree": "n_agree", "type A": "n_type_A", "type B": "n_type_B"})
                .value_counts().to_dict() == {k: v for k, v in recount.items() if k != "n_compared" and v})
        str_counts[key] = dict(recount, line=int(row.ln))
    # the atlases' Y calls in the DNA table are those of the scores and cut-offs read here
    emap = st_all.set_index("accession")
    for r in e.itertuples():
        call = emap.loc[r.cellosaurus_id, "EmeraldBay_Y"]
        assert call == ("Y expressed" if r.y_score >= e_yc else "Y below cut-off"), r.cell_name

    rows = []        # one row per line (or line and experiment) and panel
    # A: Tahoe-100M against external RNA-seq
    cmap = c.set_index("RRID")
    for r in t.itertuples():
        if r.cellosaurus_id not in cmap.index or cmap.loc[r.cellosaurus_id, "pattern"] is None:
            continue
        q = cmap.loc[r.cellosaurus_id]
        rows.append(dict(panel="A", name=r.cell_name, accession=r.cellosaurus_id, row_call=r.pattern,
                         col_call=q.pattern, row_source=f"{sources['tahoe'].name}:{r.ln}",
                         col_source=f"{sources['c_labels'].name}:{q.ln} ({q.expressed_state_source})", note=""))
    # B: EmeraldBay against Tahoe-100M
    tmap = t.set_index("cellosaurus_id")
    for r in e.itertuples():
        if r.cellosaurus_id not in tmap.index:
            continue
        q = tmap.loc[r.cellosaurus_id]
        assert q.cell_name == r.cell_name, (q.cell_name, r.cell_name)
        conf = "confident in both" if (r.confidence in CONFIDENT and q.confidence in CONFIDENT) else \
            "low confidence in one atlas"
        rows.append(dict(panel="B", name=r.cell_name, accession=r.cellosaurus_id, row_call=r.pattern,
                         col_call=q.pattern, row_source=f"{sources['emerald'].name}:{r.ln}",
                         col_source=f"{sources['tahoe'].name}:{q.ln}", note=conf))
    # C: sci-Plex (each line in each experiment) against external RNA-seq
    for r in sp.itertuples():
        q = cmap.loc[r.cellosaurus_id]
        rows.append(dict(panel="C", name=f"{r.cell_name} in {r.sample}", accession=r.cellosaurus_id,
                         row_call=r.pattern, col_call=q.pattern, row_source=f"{sources['sciplex'].name}:{r.ln}",
                         col_source=f"{sources['c_labels'].name}:{q.ln} ({q.expressed_state_source})", note=""))
    # D: Tahoe-100M Y score against STR Amelogenin
    smap = st.set_index("accession")
    assert smap.index.is_unique
    for r in t.itertuples():
        q = smap.loc[r.cellosaurus_id]
        cls = q.amelogenin_class
        if cls not in STR_LEVELS:
            rows.append(dict(panel="D", name=r.cell_name, accession=r.cellosaurus_id, row_call=r.ycall,
                             col_call=cls, row_source=f"{sources['tahoe'].name}:{r.ln}",
                             col_source=f"{sources['str'].name}:{q.ln}", note="not compared"))
            continue
        assert (q["Tahoe-100M_Y"] == "Y expressed") == (r.ycall == Y_LEVELS[0]), r.cell_name
        rows.append(dict(panel="D", name=r.cell_name, accession=r.cellosaurus_id, row_call=r.ycall, col_call=cls,
                         row_source=f"{sources['tahoe'].name}:{r.ln}", col_source=f"{sources['str'].name}:{q.ln}",
                         note=""))
    lines = pd.DataFrame(rows)
    lines = lines[lines.note != "not compared"].copy()
    diag = {"A": dict(zip(PATTERNS, PATTERNS)), "B": dict(zip(PATTERNS, PATTERNS)), "C": dict(zip(PATTERNS, PATTERNS)),
            "D": dict(zip(Y_LEVELS, STR_LEVELS))}
    lines["same"] = [diag[p].get(a) == b for p, a, b in zip(lines.panel, lines.row_call, lines.col_call)]
    cuts = dict(tahoe=(t_yc, t_xc), emerald=(e_yc, e_xc), sciplex=(s_yc, s_xc),
                ext={k: v[0] for k, v in ccut.items()})
    return dict(lines=lines, cuts=cuts, str_counts=str_counts)
