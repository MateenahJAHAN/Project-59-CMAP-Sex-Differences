"""Step 4: states, cut-offs, ambient, mosaic tests, identity check and agreement (plan sections 3-10).
Uses vehicle cells only, except the identity check, which uses label-free deep-cell profile sums."""
import sys, os, glob, json, zlib, time
import numpy as np, pandas as pd
from scipy.stats import poisson
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "verification/source_replay/part_a/scripts"))
from common import (MARKERS, Y_GENES, YP_GENES, YQ_GENES, X_GAMETOLOGS, SEED, NBINS, NGROUP, MIN_UNIT_CELLS,
                    CCLE_CUT, STATES, ARM_STATES, pseudobulk_scores, kde_cutoff, call_state, call_arm,
                    em_mixture, mh_or, l2cpm)

ROOT = sys.argv[1]
RES = os.path.join(ROOT, "results/F1"); os.makedirs(RES, exist_ok=True)
SCORES = ["XIST", "y_score", "yp_score", "yq_score"]
PILOT_A549, PILOT_SW48, SHARE_TOL = 0.427, 0.437, 0.15
consort = []


def C(dataset, step, n_cells, n_lines, note=""):
    consort.append(dict(dataset=dataset, step=step, n_cells=int(n_cells), n_lines=int(n_lines), note=note))


def seed_for(*parts):
    return [SEED, zlib.crc32("|".join(map(str, parts)).encode())]


def natkey(s):
    import re
    m = re.findall(r"\d+", str(s))
    return (int(m[-1]) if m else 0, str(s))


# ---------------- reference tables ----------------
A = pd.read_csv(os.path.join(ROOT, "inputs/A_cellline_sex_states.tsv"), sep="\t")
Cc = pd.read_csv(os.path.join(ROOT, "inputs/C_cellline_sex_labels.tsv"), sep="\t")
T = pd.read_csv(os.path.join(ROOT, "inputs/tahoe_pilot_states.tsv"), sep="\t")
assert Cc.RRID.dropna().is_unique and T.cell_line_id.is_unique
A_DUP = sorted(A.cellosaurus_id[A.cellosaurus_id.duplicated(keep=False) & A.cellosaurus_id.notna()].unique())
A_map = A.dropna(subset=["cellosaurus_id"]).drop_duplicates("cellosaurus_id", keep=False).set_index("cellosaurus_id")
assert A_map.index.is_unique
c_sex = Cc.dropna(subset=["RRID"]).set_index("RRID").cellosaurus_sex
t_sex = T.set_index("cell_line_id").pie_sex


def annotation(cvcl):
    if cvcl in c_sex.index and c_sex[cvcl] in ("Male", "Female"):
        return c_sex[cvcl], "C_cellosaurus_sex"
    if cvcl in t_sex.index and pd.notna(t_sex[cvcl]):
        v = t_sex[cvcl]
        return ("Unspecified" if v == "Sex unspecified" else v), "pilot_pie_sex"
    if cvcl in c_sex.index:
        return c_sex[cvcl], "C_cellosaurus_sex"
    return "Unknown", "none"


# ---------------- generic per-dataset analysis ----------------
def unit_and_line_calls(veh, ann, cuts, dataset, has_batches=True):
    """veh: vehicle cells with columns line, batch, depth, markers. ann: line -> (annotation, source).
    Returns line table and unit table."""
    units = pseudobulk_scores(veh, ["line", "batch"])
    units["qualifies"] = (units.n_cells >= MIN_UNIT_CELLS) & has_batches
    lines = pseudobulk_scores(veh, ["line"])
    lines["annotation"] = lines.line.map(lambda l: ann[l][0]); lines["annotation_source"] = lines.line.map(lambda l: ann[l][1])
    rows = []
    for _, r in lines.iterrows():
        ann, ann_src = r["annotation"], r["annotation_source"]
        st = call_state(ann, r.XIST, r.y_score, cuts["XIST"], cuts["y_score"])
        grid = {call_state(ann, r.XIST, r.y_score, cuts["XIST"] + dx, cuts["y_score"] + dy)
                for dx in (-1, 0, 1) for dy in (-1, 0, 1)}
        u = units[(units.line == r.line) & units.qualifies]
        ucalls = [call_state(ann, a, b, cuts["XIST"], cuts["y_score"]) for a, b in zip(u.XIST, u.y_score)]
        if has_batches:
            n_pl = len(u); agree = sum(c == st for c in ucalls)
        else:
            n_pl = 1; agree = 1
        gs = len(grid) == 1
        if n_pl == 0:
            conf = "low"
        elif agree == n_pl and gs:
            conf = "high"
        elif agree / n_pl >= 0.8 and has_batches:
            conf = "medium"
        else:
            conf = "low"
        arm = call_arm(r.yp_score, r.yq_score, cuts["yp_score"], cuts["yq_score"])
        y_neg = (r.y_score < cuts["y_score"]) and all(u.y_score < cuts["y_score"])
        x_neg = (r.XIST < cuts["XIST"]) and all(u.XIST < cuts["XIST"])
        rows.append(dict(line=r.line, state=st, arm_state=arm, confidence=conf, n_plates=n_pl, plates_agree=agree,
                         grid_stable=gs, grid_states=";".join(sorted(grid)), y_negative_line=y_neg, x_negative_line=x_neg,
                         unit_states=";".join(ucalls)))
    calls = pd.DataFrame(rows)
    assert set(calls.state) <= STATES and set(calls.arm_state) <= ARM_STATES
    out = lines.merge(calls, on="line", how="left", validate="1:1")
    assert len(out) == len(lines)
    units["unit_state"] = [call_state(lines.set_index("line").annotation[l], a, b, cuts["XIST"], cuts["y_score"])
                           for l, a, b in zip(units.line, units.XIST, units.y_score)]
    units["dataset"] = dataset
    return out, units


def backgrounds(veh, lines_tab, marker_cols, neg_flag, require_arm_lost=False):
    """Plan 6: marker UMI rate in vehicle cells of marker-negative lines, per bg_batch.
    require_arm_lost (sensitivity, see report Deviations): Y-negative lines must also have arm state 'Y lost'."""
    sel = lines_tab[neg_flag].astype(bool)
    if require_arm_lost:
        sel &= lines_tab.arm_state == "Y lost"
    neg_lines = set(lines_tab.loc[sel, "line"])
    v = veh[veh.line.isin(neg_lines)]
    s = v.assign(m=v[marker_cols].sum(axis=1)).groupby("bg_batch").agg(umi=("m", "sum"), depth=("depth", "sum"), n_cells=("m", "size"))
    s["rate"] = np.where(s.umi > 0, s.umi / s.depth, 0.5 / s.depth)
    s["zero_floor_used"] = s.umi == 0
    return s, neg_lines


def mosaic_for_line(vl, marker, bgmap, batches_odd_even, dataset, line, rng_parts):
    """(a), (b), (c) for one line's vehicle cells vl. marker: 'Y' or 'XIST'."""
    cols = Y_GENES if marker == "Y" else ["XIST"]
    y = vl[cols].sum(axis=1).values.astype(float); d = vl.depth.values.astype(float)
    b = vl.bg_batch.map(bgmap).values.astype(float)
    assert np.isfinite(b).all() and (b > 0).all()
    fit = em_mixture(y, d, b)
    res = dict(dataset=dataset, line=line, marker=marker, n_vehicle_cells=len(vl), mean_depth=float(d.mean()),
               marker_umi_per_1k=float(1000 * y.sum() / d.sum()), background_per_1k=float(1000 * np.average(b, weights=d)),
               pi_neg=fit["pi"], lam_per_1k=1000 * fit["lam"], em_iters=fit["iters"],
               loglik_gain=fit["loglik"] - fit["loglik_one"], zero_frac=float(np.mean(y == 0)))
    odd, even = batches_odd_even
    for nm, bs in (("odd", odd), ("even", even)):
        m = vl.batch.isin(bs).values
        res[f"n_{nm}"] = int(m.sum()); res[f"{nm}_batches"] = ";".join(map(str, bs))
        res[f"pi_{nm}"] = em_mixture(y[m], d[m], b[m])["pi"] if m.sum() >= 20 else np.nan
    res["stability_diff"] = abs(res["pi_odd"] - res["pi_even"]) if np.isfinite(res["pi_odd"]) and np.isfinite(res["pi_even"]) else np.nan
    res["b_holds"] = bool(res["stability_diff"] <= 0.10) if np.isfinite(res["stability_diff"]) else None
    if marker == "Y":
        strata = pd.qcut(d, 10, labels=False, duplicates="drop")
        yq = vl[YQ_GENES].sum(axis=1).values; yp = vl[YP_GENES].sum(axis=1).values
        res["or_y"] = mh_or(yq == 0, yp == 0, strata)
        X = vl[X_GAMETOLOGS].sum(axis=1).values.astype(np.int64)
        rng = np.random.default_rng(rng_parts)
        f = min(1.0, y.sum() / max(X.sum(), 1)); q = yq.sum() / max(y.sum(), 1)
        xt = rng.binomial(X, f); xa = rng.binomial(xt, q); xb = xt - xa
        res["or_x_matched"] = mh_or(xa == 0, xb == 0, strata)
        xa2 = rng.binomial(X, 0.5); xb2 = X - xa2
        res["or_x_unthinned_5050"] = mh_or(xa2 == 0, xb2 == 0, strata)
        res["x_thin_factor"] = f; res["yq_share"] = q
        res["frac_yq0"] = float(np.mean(yq == 0)); res["frac_yp0"] = float(np.mean(yp == 0))
        if np.isfinite(res["or_y"]) and np.isfinite(res["or_x_matched"]) and res["or_x_matched"] > 0:
            res["c_holds"] = bool(res["or_y"] > 2 * res["or_x_matched"])
        else:
            res["c_holds"] = None
    return res


def identity(acc, cnt, lines, line, marker, keep_genes, others=None):
    """Plan 7d. acc rows = line_index * NGROUP + group. others: optional dict name -> profile (cross-experiment)."""
    li = lines.index(line); base = li * NGROUP
    if marker == "Y":
        neg_g, pos_g = 0, list(range(1, 1 + NBINS))
    else:
        neg_g, pos_g = 1 + NBINS, list(range(2 + NBINS, 2 + 2 * NBINS))
    n_neg = int(cnt[base + neg_g]); n_pos_bins = cnt[[base + g for g in pos_g]]; n_pos = int(n_pos_bins.sum())
    res = dict(n_deep_neg=n_neg, n_deep_pos=n_pos)
    if n_neg < 20 or n_pos < 40:
        res.update(d_holds=None, d_note="fewer than 20 deep negative cells" if n_neg < 20 else "too few positive cells")
        return res

    def prof(rows):
        v = acc[rows].sum(0).astype(float)[keep_genes]
        return np.log2(1e6 * v / v.sum() + 1)

    k = int(np.clip(round(NBINS * n_neg / n_pos), 1, NBINS // 2))
    neg = prof([base + neg_g])
    s1 = prof([base + g for g in pos_g[:k]]); s2 = prof([base + g for g in pos_g[k:2 * k]])
    posall = prof([base + g for g in pos_g])
    r_ns2 = np.corrcoef(neg, s2)[0, 1]; r_split = np.corrcoef(s1, s2)[0, 1]; r_nall = np.corrcoef(neg, posall)[0, 1]
    other = {}
    for j, l2 in enumerate(lines):
        if l2 == line:
            continue
        other[l2] = np.corrcoef(neg, prof([j * NGROUP + g for g in range(0, 1 + NBINS)]))[0, 1]
    if others:
        for nm, p in others.items():
            other[nm] = np.corrcoef(neg, p)[0, 1]
    best = max(other, key=other.get) if other else None
    res.update(k_bins=k, n_s1=int(n_pos_bins[:k].sum()), n_s2=int(n_pos_bins[k:2 * k].sum()), r_neg_s2=r_ns2,
               r_split_s1_s2=r_split, r_neg_ownpos=r_nall, best_other_line=best,
               r_neg_best_other=other.get(best, np.nan), n_other_lines=len(other))
    res["d1_holds"] = bool(r_ns2 >= r_split - 0.01)
    res["d2_holds"] = bool(r_nall - res["r_neg_best_other"] >= 0.10) if best else None
    res["d_holds"] = (res["d1_holds"] and res["d2_holds"]) if res["d2_holds"] is not None else None
    return res


def flag(r):
    crit = [r.get("a_holds"), r.get("b_holds"), r.get("d_holds")] + ([r.get("c_holds")] if r["marker"] == "Y" else [])
    if any(c is False for c in crit):
        return "not mosaic"
    if any(c is None for c in crit):
        return "not evaluable"
    return "mosaic"


def foreign_shares(veh, lines_tab, ybg, xbg):
    out = []
    for _, r in lines_tab.iterrows():
        vl = veh[veh.line == r.line]
        ysum = vl[Y_GENES].sum(axis=1).values; xs = vl["XIST"].values; d = vl.depth.values
        fy = fx = ey = ex_ = np.nan
        flags = np.zeros(len(vl), bool)
        if r.y_negative_line:
            fy = float(np.mean(ysum >= 3)); flags |= ysum >= 3
            ey = float(np.mean(poisson.sf(2, vl.bg_batch.map(ybg).values * d)))
        if r.x_negative_line:
            fx = float(np.mean(xs >= 3)); flags |= xs >= 3
            ex_ = float(np.mean(poisson.sf(2, vl.bg_batch.map(xbg).values * d))) if xbg is not None else np.nan
        fc = float(flags.mean()) if (r.y_negative_line or r.x_negative_line) else np.nan
        out.append(dict(line=r.line, foreign_share_y=fy, foreign_share_y_expected_ambient=ey, foreign_share_xist=fx,
                        foreign_share_xist_expected_ambient=ex_, foreign_cell_share=fc))
    return pd.DataFrame(out)


def analyze(dataset, veh, ann_map, cuts, acc, cnt, acc_lines, keep_genes, has_batches, odd_even_fn,
            bg_override=None, others=None, force=(), min_cells=1000, arm_bg=False):
    veh = veh.copy()
    ann = {l: ann_map[l] for l in veh.line.unique()}
    lt, units = unit_and_line_calls(veh, ann, cuts, dataset, has_batches)
    # ambient
    ybg, yneg = backgrounds(veh, lt, Y_GENES, "y_negative_line", require_arm_lost=arm_bg)
    xbg, xneg = backgrounds(veh, lt, ["XIST"], "x_negative_line")
    amb_rows = []
    for nm, bg, negl, cols in (("Y", ybg, yneg, Y_GENES), ("XIST", xbg, xneg, ["XIST"])):
        for bb, r in bg.iterrows():
            amb_rows.append(dict(dataset=dataset, marker=nm, level="batch_background", bg_batch=bb, line="(all negative lines)",
                                 n_cells=int(r.n_cells), marker_umis=float(r.umi), depth=float(r.depth),
                                 rate_per_1k=1000 * r.rate, zero_floor_used=bool(r.zero_floor_used)))
        v = veh[veh.line.isin(negl)]
        g = v.assign(m=v[cols].sum(axis=1)).groupby(["line", "bg_batch"]).agg(umi=("m", "sum"), depth=("depth", "sum"), n=("m", "size"))
        for (l, bb), r in g.iterrows():
            amb_rows.append(dict(dataset=dataset, marker=nm, level="line_x_batch", bg_batch=bb, line=l, n_cells=int(r.n),
                                 marker_umis=float(r.umi), depth=float(r.depth), rate_per_1k=1000 * r.umi / r.depth,
                                 zero_floor_used=False))
    ybgmap = ybg.rate.to_dict(); xbgmap = xbg.rate.to_dict()
    src_y = src_x = "own"
    if bg_override is not None:
        if len(ybgmap) == 0:
            ybgmap = {b: bg_override["Y"] for b in veh.bg_batch.unique()}; src_y = "override"
        if len(xbgmap) == 0:
            xbgmap = {b: bg_override["XIST"] for b in veh.bg_batch.unique()}; src_x = "override"
    # batches without a negative line: fall back to the dataset pooled rate
    for mp_, bg in ((ybgmap, ybg), (xbgmap, xbg)):
        pooled = (bg.umi.sum() / bg.depth.sum()) if len(bg) and bg.umi.sum() > 0 else None
        for b in veh.bg_batch.unique():
            if b not in mp_ and pooled is not None:
                mp_[b] = pooled
    fs = foreign_shares(veh, lt, ybgmap, xbgmap if len(xbgmap) else None)
    lt = lt.merge(fs, on="line", how="left", validate="1:1")
    # mosaic tests
    mos = []
    nveh = veh.line.value_counts()
    for _, r in lt.iterrows():
        l = r.line
        y_hi = r.y_score >= cuts["y_score"]
        tests = []
        if nveh[l] >= min_cells or l in force:
            if y_hi and r.state in ("XY", "discordant"):
                tests.append(("Y", "test"))
            elif nveh[l] >= min_cells or ("Y", l) in force:
                tests.append(("Y", "control"))
            if r.state == "XX-like":
                tests.append(("XIST", "test"))
            elif nveh[l] >= min_cells or ("XIST", l) in force:
                tests.append(("XIST", "control"))
        for marker, role in tests:
            if (marker, l) in force and nveh[l] < min_cells:
                role = "forced (below 1,000 cells)"
            bgm = ybgmap if marker == "Y" else xbgmap
            if not bgm:
                continue
            vl = veh[veh.line == l]
            oe = odd_even_fn(vl)
            res = mosaic_for_line(vl, marker, bgm, oe, dataset, l, seed_for(dataset, l, marker))
            res["role"] = role; res["state"] = r.state; res["confidence"] = r.confidence
            res["background_source"] = (src_y if marker == "Y" else src_x) + ("; Y-negative lines restricted to arm state Y lost (sensitivity)" if arm_bg and marker == "Y" else "")
            res["a_holds"] = bool(0.15 <= res["pi_neg"] <= 0.85)
            res.update(identity(acc, cnt, acc_lines, l, marker, keep_genes, others))
            res["mosaic_flag"] = flag(res) if role == "test" else ("not assigned (" + role + ")")
            mos.append(res)
    mos = pd.DataFrame(mos)
    return lt, units, pd.DataFrame(amb_rows), mos, ybgmap, xbgmap


def keep_gene_mask(acc, cnt, nlines, exclude_cols, extra_mask=None):
    tot = np.zeros(acc.shape[1])
    for li in range(nlines):
        tot += acc[[li * NGROUP + g for g in range(0, 1 + NBINS)]].sum(0)
    cpm = 1e6 * tot / tot.sum()
    m = cpm >= 10
    m[list(exclude_cols)] = False
    if extra_mask is not None:
        m &= extra_mask
    return m


def main():
    t0 = time.time()
    runlog = []
    # ============ EmeraldBay ============
    meta = pd.read_parquet(os.path.join(ROOT, "data/F1/emerald_cell_meta.parquet"))
    C("EmeraldBay", "all cells in 116 shards (metadata scan)", len(meta), meta.cell_line.nunique())
    tf_n = int((meta.drug == "DMSO_TF").sum()); t0_n = int((meta.drug == "DMSO_T0").sum())
    files = sorted(glob.glob(os.path.join(ROOT, "data/F1/emerald_proc/cells_*.parquet")))
    cells = pd.concat([pd.read_parquet(f) for f in files], ignore_index=True)
    assert cells.BARCODE_SUB_LIB_ID.is_unique
    logs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(ROOT, "data/F1/emerald_proc/shardlog_*.json")))]
    assert len(logs) == 116 and sum(l["rows"] for l in logs) == len(meta)
    special = sum(l["special_entries"] for l in logs); negs = sum(l["negative_entries"] for l in logs)
    nonint = sum(l["nonint_entries"] for l in logs)
    # join check against the metadata scan
    chk = cells[["BARCODE_SUB_LIB_ID", "cell_line", "drug", "sample"]].merge(
        meta[["BARCODE_SUB_LIB_ID", "cell_line", "drug", "sample"]], on="BARCODE_SUB_LIB_ID", how="left", validate="1:1",
        suffixes=("", "_m"))
    assert len(chk) == len(cells) and (chk.cell_line == chk.cell_line_m).all() and (chk.drug == chk.drug_m).all()
    assert (cells.drug == "DMSO_TF").sum() == tf_n and (cells.drug == "DMSO_T0").sum() == t0_n
    veh_all = cells.rename(columns={"cell_line": "line"})
    veh_all["batch"] = veh_all["sample"]; veh_all["bg_batch"] = veh_all["sample"]
    veh = veh_all[veh_all.drug == "DMSO_TF"].copy()
    C("EmeraldBay", "vehicle cells (drug == DMSO_TF)", len(veh), veh.line.nunique(), "DMSO_T0 excluded from primary analyses")
    C("EmeraldBay", "quality filter", len(veh), veh.line.nunique(), "none provided by the dataset; no cells removed")
    assert (veh.depth > 0).all()
    cl = pd.read_parquet(os.path.join(ROOT, "data/F1/meta/cell_line_metadata.parquet")).drop_duplicates("Cell_ID_Cellosaur")
    assert cl.Cell_ID_Cellosaur.is_unique and len(cl) == 52
    ann_map = {l: annotation(l) for l in cl.Cell_ID_Cellosaur}
    # cut-offs from line x sample units
    units0 = pseudobulk_scores(veh, ["line", "sample"])
    q = units0[units0.n_cells >= MIN_UNIT_CELLS]
    C("EmeraldBay", "line x sample units with >= 10 vehicle cells (enter KDE)", q.n_cells.sum(), q.line.nunique(),
      f"{len(q)} of {len(units0)} units")
    cut_rows = []; cuts = {}
    for s in SCORES:
        c, info = kde_cutoff(q[s].values)
        if not np.isfinite(c):
            c = CCLE_CUT[s]
        cuts[s] = c
        cut_rows.append(dict(dataset="EmeraldBay", score=s, cutoff=c, source=info.get("method"), **{k: v for k, v in info.items() if k != "method"}))
    # identity accumulators
    accs = sorted(glob.glob(os.path.join(ROOT, "data/F1/emerald_proc/acc_w*.npz")))
    acc = None; cnt = None; done = []
    for f in accs:
        z = np.load(f)
        acc = z["acc"].astype(np.float64) if acc is None else acc + z["acc"]
        cnt = z["cnt"] if cnt is None else cnt + z["cnt"]; done += list(z["done"])
    assert len(done) == 116 and len(set(done)) == 116, len(done)
    acc_lines = sorted(cl.Cell_ID_Cellosaur)
    g = pd.read_parquet(os.path.join(ROOT, "data/F1/meta/gene_metadata.parquet"))
    excl = g[g.gene_symbol.isin(["XIST"] + Y_GENES)].token_id.values
    valid = np.zeros(acc.shape[1], bool); valid[g.token_id.values] = True
    keep = keep_gene_mask(acc, cnt, len(acc_lines), excl, valid)
    n_deep = int(cnt[[li * NGROUP + gg for li in range(len(acc_lines)) for gg in range(0, 1 + NBINS)]].sum())
    C("EmeraldBay", "identity check: deep cells (>= 2,500 UMIs), all conditions", n_deep, len(acc_lines),
      f"{int(keep.sum())} genes with CPM >= 10 used")

    def oe_em(vl):
        bs = sorted(vl.batch.unique(), key=natkey)
        return [bs[i] for i in range(0, len(bs), 2)], [bs[i] for i in range(1, len(bs), 2)]

    force = {("Y", "CVCL_0023"), ("XIST", "CVCL_1724"), "CVCL_0023", "CVCL_1724"}
    lt, units, amb, mos, ybgmap, xbgmap = analyze("EmeraldBay", veh, ann_map, cuts, acc, cnt, acc_lines, keep, True, oe_em, force=force)
    _, _, amb_sens, mos_sens, ybg_sens, _ = analyze("EmeraldBay", veh, ann_map, cuts, acc, cnt, acc_lines, keep, True, oe_em,
                                                     force=force, arm_bg=True)
    mos_sens = mos_sens[mos_sens.marker == "Y"].copy()
    # exploratory: DMSO_T0
    t0v = veh_all[veh_all.drug == "DMSO_T0"].copy()
    expl = []
    for l, mk in (("CVCL_0023", "Y"), ("CVCL_1724", "XIST")):
        vl = t0v[t0v.line == l]
        cols = Y_GENES if mk == "Y" else ["XIST"]
        bg = ybgmap if mk == "Y" else xbgmap
        # T0 sample has its own batch; use the pooled vehicle background of the TF samples (exploratory)
        b = np.full(len(vl), np.mean(list(bg.values())))
        if len(vl) >= 20:
            fit = em_mixture(vl[cols].sum(axis=1).values, vl.depth.values, b)
            expl.append(dict(dataset="EmeraldBay DMSO_T0 (exploratory)", line=l, marker=mk, n_cells=len(vl),
                             samples=";".join(sorted(vl["sample"].unique())), pi_neg=fit["pi"], lam_per_1k=1000 * fit["lam"],
                             background_per_1k=1000 * b[0], note="background = mean of DMSO_TF sample backgrounds"))
    # states table
    clm = cl.set_index("Cell_ID_Cellosaur")
    st = lt.rename(columns={"line": "cellosaurus_id"})
    st.insert(0, "cell_name", st.cellosaurus_id.map(clm.cell_name))
    st.insert(2, "organ", st.cellosaurus_id.map(clm.Organ))
    st = st.rename(columns={"n_cells": "n_vehicle_cells"})
    cols = ["cell_name", "cellosaurus_id", "organ", "annotation", "annotation_source", "n_vehicle_cells", "n_plates",
            "XIST", "y_score", "yp_score", "yq_score", "state", "arm_state", "confidence", "plates_agree", "grid_stable",
            "foreign_cell_share", "foreign_share_y", "foreign_share_y_expected_ambient", "foreign_share_xist",
            "foreign_share_xist_expected_ambient", "y_negative_line", "x_negative_line", "grid_states", "unit_states",
            "depth"] + [s + "_l2cpm" for s in MARKERS]
    st = st[cols].sort_values(["state", "cell_name"])
    assert len(st) == 52 and st.cellosaurus_id.is_unique
    st.to_csv(os.path.join(RES, "F1_emerald_states.tsv"), sep="\t", index=False)
    units.to_csv(os.path.join(RES, "F1_emerald_units.tsv"), sep="\t", index=False)
    mos.insert(1, "cell_name", mos.line.map(clm.cell_name))
    mos.to_csv(os.path.join(RES, "F1_emerald_mosaic.tsv"), sep="\t", index=False)
    mos_sens.insert(1, "cell_name", mos_sens.line.map(clm.cell_name))
    mos_sens.to_csv(os.path.join(RES, "F1_emerald_mosaic_sensitivity_background.tsv"), sep="\t", index=False)
    amb_sens[(amb_sens.marker == "Y")].assign(analysis="sensitivity: Y-negative lines restricted to arm state Y lost").to_csv(
        os.path.join(RES, "F1_ambient_sensitivity_background.tsv"), sep="\t", index=False)
    C("EmeraldBay", "lines with >= 1,000 vehicle cells (mosaic test)", veh[veh.line.isin(veh.line.value_counts()[lambda s: s >= 1000].index)].shape[0],
      int((veh.line.value_counts() >= 1000).sum()))
    em_lt = st
    # per-cell vehicle marker table (for figures and re-use)
    out_cells = veh_all[["BARCODE_SUB_LIB_ID", "line", "drug", "sample", "depth", "n_genes"] + MARKERS]
    out_cells.to_csv(os.path.join(RES, "F1_emerald_vehicle_cells.tsv.gz"), sep="\t", index=False, compression="gzip")

    if os.environ.get("F1_EMERALD_ONLY"):  # debugging switch; run_all.sh never sets it
        print(mos.to_string()); return
    # ============ sci-Plex ============
    sp_states, sp_mos, sp_units, sp_amb, sp_mos_sens = [], [], [], [], []
    name_to_rrid = {}
    for nm, cname in (("A549", "A549"), ("MCF7", "MCF7"), ("K562", "K-562")):
        hit = Cc[Cc.model_name == cname]
        assert len(hit) == 1, cname
        name_to_rrid[nm] = hit.RRID.iloc[0]
    sp_ann = {nm: annotation(r) for nm, r in name_to_rrid.items()}
    em_cut_tab = cuts.copy()
    sp_cuts = CCLE_CUT.copy()
    for s in SCORES:
        cut_rows.append(dict(dataset="sci-Plex (all samples)", score=s, cutoff=sp_cuts[s], source="Task A CCLE (fallback; too few lines for KDE)"))
    bg3 = None; bg4 = None; prof3 = {}
    for tag in ("sciPlex4", "sciPlex3", "sciPlex2"):
        sm = json.load(open(os.path.join(ROOT, f"data/F1/sciplex_proc/{tag}_labels_summary.json")))
        v = pd.read_parquet(os.path.join(ROOT, f"data/F1/sciplex_proc/{tag}_vehicle_cells.parquet"))
        z = np.load(os.path.join(ROOT, f"data/F1/sciplex_proc/{tag}_acc.npz"), allow_pickle=True)  # files written by process_sciplex.py in this run
        acc_s, cnt_s, lines_s = z["acc"], z["cnt"], list(z["lines"])
        human = z["human"]; sym = z["symbol"]; gid = z["gene_id"]
        C(tag, "cells in matrix / pData", sm["n_cells_file"], 0, f"pData UMI floor (min n.umi) = {sm['min_n_umi_pdata']}")
        C(tag, "cells with a line label (hash filter)", sm["n_assigned"], len(sm["assigned_by_line"]))
        C(tag, "vehicle cells", len(v), v.line.nunique())
        assert len(v) == sm["n_vehicle"] and v.cell.is_unique
        assert set(v.line) <= set(name_to_rrid)
        excl = np.where(np.isin(sym, ["XIST"] + Y_GENES) & human)[0]
        keep_s = keep_gene_mask(acc_s, cnt_s, len(lines_s), excl, human)
        nd = int(cnt_s[[li * NGROUP + gg for li in range(len(lines_s)) for gg in range(0, 1 + NBINS)]].sum())
        C(tag, "identity check: deep cells (>= 2,500 UMIs), all conditions", nd, len(lines_s), f"{int(keep_s.sum())} genes used")
        if tag == "sciPlex2":
            def oe(vl):
                wells = sorted(vl.well.unique()); r = np.random.default_rng(seed_for("sciPlex2", "wells")); r.shuffle(wells)
                h = len(wells) // 2
                return wells[:h], wells[h:]
            v = v.assign(batch=v.well)  # batch = hash well, used only for the random-half stability test
            hb = False
            # cross-experiment comparison profiles: sciPlex3 MCF7 and K562, aligned on gene id
            others = {}
            idx3 = {gname: i for i, gname in enumerate(prof3["gene_id"])}
            take = np.array([idx3.get(x, -1) for x in gid])
            for l2 in ("MCF7", "K562"):
                p3 = prof3[l2]
                vec = np.where(take >= 0, p3[np.clip(take, 0, None)], 0.0)
                vv = vec[keep_s]; others["sciPlex3 " + l2] = np.log2(1e6 * vv / vv.sum() + 1)
            override = bg3
        else:
            def oe(vl):
                bs = sorted(vl.batch.unique(), key=natkey)
                return [bs[i] for i in range(0, len(bs), 2)], [bs[i] for i in range(1, len(bs), 2)]
            hb = True; others = None
            # sciPlex3 has no XIST-negative line, so its XIST background is the pooled sciPlex4 rate (deviation, see report)
            override = {"Y": np.nan, "XIST": bg4["XIST"]} if tag == "sciPlex3" else None
        lt_s, un_s, amb_s, mos_s, yb, xb = analyze(tag, v, sp_ann, sp_cuts, acc_s, cnt_s, lines_s, keep_s, hb, oe,
                                                     bg_override=override, others=others, force={("Y", "A549"), "A549"})
        _, _, _, mos_s2, _, _ = analyze(tag, v, sp_ann, sp_cuts, acc_s, cnt_s, lines_s, keep_s, hb, oe,
                                        bg_override=override, others=others, force={("Y", "A549"), "A549"}, arm_bg=True)
        sp_mos_sens.append(mos_s2[mos_s2.marker == "Y"])
        if tag == "sciPlex2":
            lt_s["n_plates"] = 1; lt_s["note"] = "no plate field: one pseudobulk per line; stability uses random halves of hash wells"
        # secondary: EmeraldBay cut-offs
        lt_s["state_emerald_cutoffs"] = [call_state(sp_ann[l][0], a, b, em_cut_tab["XIST"], em_cut_tab["y_score"]) for l, a, b in zip(lt_s.line, lt_s.XIST, lt_s.y_score)]
        lt_s["arm_state_emerald_cutoffs"] = [call_arm(a, b, em_cut_tab["yp_score"], em_cut_tab["yq_score"]) for a, b in zip(lt_s.yp_score, lt_s.yq_score)]
        lt_s.insert(0, "sample", tag)
        sp_states.append(lt_s); sp_mos.append(mos_s); un_s["dataset"] = tag; sp_units.append(un_s); sp_amb.append(amb_s)
        if tag == "sciPlex4":
            ab4 = amb_s[amb_s.level == "batch_background"]
            bg4 = {mk: ab4[ab4.marker == mk].marker_umis.sum() / ab4[ab4.marker == mk].depth.sum() for mk in ("Y", "XIST")}
        if tag == "sciPlex3":
            # pooled backgrounds for sciPlex2 / sciPlex4 fallback
            ab = amb_s[amb_s.level == "batch_background"]
            bg3 = {mk: ab[ab.marker == mk].marker_umis.sum() / ab[ab.marker == mk].depth.sum() for mk in ("Y", "XIST")}
            prof3["gene_id"] = gid
            for l2 in lines_s:
                li = lines_s.index(l2)
                prof3[l2] = acc_s[[li * NGROUP + gg for gg in range(0, 1 + NBINS)]].sum(0).astype(float)
    sps = pd.concat(sp_states, ignore_index=True)
    sps = sps.rename(columns={"line": "cell_name", "n_cells": "n_vehicle_cells"})
    sps.insert(2, "cellosaurus_id", sps.cell_name.map(name_to_rrid))
    sps.insert(3, "organ", sps.cellosaurus_id.map(Cc.set_index("RRID").tissue))
    keepc = ["sample", "cell_name", "cellosaurus_id", "organ", "annotation", "annotation_source", "n_vehicle_cells", "n_plates",
             "XIST", "y_score", "yp_score", "yq_score", "state", "arm_state", "confidence", "plates_agree", "grid_stable",
             "foreign_cell_share", "foreign_share_y", "foreign_share_y_expected_ambient", "foreign_share_xist",
             "foreign_share_xist_expected_ambient", "y_negative_line", "x_negative_line", "state_emerald_cutoffs",
             "arm_state_emerald_cutoffs", "grid_states", "unit_states", "depth"] + [s + "_l2cpm" for s in MARKERS]
    sps = sps[keepc]
    sps.to_csv(os.path.join(RES, "F1_sciplex_states.tsv"), sep="\t", index=False)
    spm = pd.concat(sp_mos, ignore_index=True)
    spm.to_csv(os.path.join(RES, "F1_sciplex_mosaic.tsv"), sep="\t", index=False)
    spm_sens = pd.concat(sp_mos_sens, ignore_index=True)
    spm_sens.to_csv(os.path.join(RES, "F1_sciplex_mosaic_sensitivity_background.tsv"), sep="\t", index=False)
    pd.concat(sp_units, ignore_index=True).to_csv(os.path.join(RES, "F1_sciplex_units.tsv"), sep="\t", index=False)
    amb_all = pd.concat([amb] + sp_amb, ignore_index=True)
    amb_all.to_csv(os.path.join(RES, "F1_ambient.tsv"), sep="\t", index=False)
    pd.DataFrame(cut_rows).to_csv(os.path.join(RES, "F1_cutoffs.tsv"), sep="\t", index=False)
    pd.DataFrame(expl).to_csv(os.path.join(RES, "F1_exploratory_T0.tsv"), sep="\t", index=False)

    # ============ agreement ============
    ag = []

    def comp(name, df, a, b, rule, thr=None, subset=""):
        df = df.dropna(subset=[a, b])
        n = len(df); k = int((df[a] == df[b]).sum())
        pct = k / n if n else np.nan
        res = "descriptive"
        if thr is not None:
            res = "pass" if (n > 0 and pct >= thr) else "fail"
        dis = df[df[a] != df[b]]
        ag.append(dict(comparison=name, subset=subset, n_compared=n, n_agree=k, pct_agree=pct, value=np.nan, reference=np.nan,
                       abs_diff=np.nan, rule=rule, result=res,
                       disagreements=";".join(f"{r.get('cell_name', '')}:{r[a]} vs {r[b]}" for _, r in dis.iterrows())))

    em = em_lt.copy()
    tp = T.set_index("cell_line_id")
    em["pilot_state"] = em.cellosaurus_id.map(tp.tahoe_state); em["pilot_conf"] = em.cellosaurus_id.map(tp.confidence)
    sh = em[em.pilot_state.notna()]
    hm = sh[sh.confidence.isin(["high", "medium"]) & sh.pilot_conf.isin(["high", "medium"])]
    comp("EmeraldBay vs Tahoe pilot state", hm, "state", "pilot_state", ">= 90% agreement (pre-registered)", 0.90,
         "lines with high or medium confidence in both")
    comp("EmeraldBay vs Tahoe pilot state", sh, "state", "pilot_state", "descriptive", None, "all shared lines")
    em["C_state"] = em.cellosaurus_id.map(Cc.dropna(subset=["RRID"]).set_index("RRID").expressed_state)
    comp("EmeraldBay vs Task C expressed_state", em[em.C_state.notna() & (em.C_state != "not_measured")], "state", "C_state",
         "descriptive (no threshold in task)", None, "lines in Task C with a measured state")
    em["A_state"] = em.cellosaurus_id.map(A_map.state)
    comp("EmeraldBay vs Task A state", em[em.A_state.notna() & (em.A_state != "no_data")], "state", "A_state",
         "descriptive", None, "lines in Task A with data")

    def share(name, mosdf, line_col, line, marker, ref, subset):
        r = mosdf[(mosdf[line_col] == line) & (mosdf.marker == marker)]
        if len(r) == 0:
            ag.append(dict(comparison=name, subset=subset, n_compared=0, rule=f"|share - {ref}| <= {SHARE_TOL}", result="not run",
                           reference=ref)); return
        r = r.iloc[0]
        dlt = abs(r.pi_neg - ref)
        ag.append(dict(comparison=name, subset=subset, n_compared=int(r.n_vehicle_cells), n_agree=np.nan, pct_agree=np.nan,
                       value=r.pi_neg, reference=ref, abs_diff=dlt, rule=f"|share - {ref}| <= {SHARE_TOL} (pre-registered)",
                       result="pass" if dlt <= SHARE_TOL else "fail",
                       disagreements=f"mosaic_flag={r.mosaic_flag}; role={r.role}"))

    share("EmeraldBay A549 Y-negative share vs Tahoe pilot", mos, "line", "CVCL_0023", "Y", PILOT_A549, "DMSO_TF vehicle cells")
    share("EmeraldBay A549 Y-negative share vs Tahoe pilot [sensitivity background]", mos_sens, "line", "CVCL_0023", "Y", PILOT_A549,
          "DMSO_TF vehicle cells; Y background from lines with arm state Y lost (deviation, sensitivity)")
    share("EmeraldBay SW48 XIST-negative share vs Tahoe pilot", mos, "line", "CVCL_1724", "XIST", PILOT_SW48, "DMSO_TF vehicle cells")
    sps2 = sps.copy()
    sps2["A_state"] = sps2.cellosaurus_id.map(A_map.state)
    sps2["C_state"] = sps2.cellosaurus_id.map(Cc.dropna(subset=["RRID"]).set_index("RRID").expressed_state)
    for tag in ("sciPlex2", "sciPlex3", "sciPlex4"):
        s_ = sps2[sps2["sample"] == tag]
        comp(f"{tag} vs Task A state", s_, "state", "A_state", "descriptive", None, ",".join(s_.cell_name))
        comp(f"{tag} vs Task C expressed_state", s_, "state", "C_state", "descriptive", None, ",".join(s_.cell_name))
        comp(f"{tag} CCLE cut-offs vs EmeraldBay cut-offs", s_, "state", "state_emerald_cutoffs", "descriptive (sensitivity)", None, ",".join(s_.cell_name))
        share(f"{tag} A549 Y-negative share vs Tahoe pilot", spm[spm.dataset == tag], "line", "A549", "Y", PILOT_A549, "vehicle cells")
        share(f"{tag} A549 Y-negative share vs Tahoe pilot [sensitivity background]", spm_sens[spm_sens.dataset == tag], "line", "A549", "Y",
              PILOT_A549, "vehicle cells; Y background from lines with arm state Y lost (deviation, sensitivity)")
    agd = pd.DataFrame(ag)
    agd.to_csv(os.path.join(RES, "F1_agreement.tsv"), sep="\t", index=False)
    pd.DataFrame(consort).to_csv(os.path.join(RES, "F1_consort.tsv"), sep="\t", index=False)
    json.dump({t: json.load(open(os.path.join(ROOT, f"data/F1/sciplex_proc/{t}_labels_summary.json"))) for t in ("sciPlex2", "sciPlex3", "sciPlex4")},
              open(os.path.join(RES, "F1_sciplex_labels_summary.json"), "w"), indent=1)
    json.dump(dict(special_token_entries=special, negative_value_entries=negs, nonint_entries=nonint,
                   n_vehicle_TF=tf_n, n_vehicle_T0=t0_n, emerald_cutoffs=cuts, ccle_cutoffs=CCLE_CUT,
                   emerald_y_background_per_1k={k: 1000 * v for k, v in ybgmap.items()},
                   emerald_y_background_per_1k_sensitivity={k: 1000 * v for k, v in ybg_sens.items()},
                   emerald_xist_background_per_1k={k: 1000 * v for k, v in xbgmap.items()},
                   sciplex3_pooled_background_per_1k={k: 1000 * v for k, v in bg3.items()},
                   sciplex4_pooled_background_per_1k={k: 1000 * v for k, v in bg4.items()},
                   n_genes_identity_emerald=int(keep.sum()), seconds=round(time.time() - t0, 1)),
              open(os.path.join(RES, "F1_run_summary.json"), "w"), indent=1, default=float)
    print(agd[["comparison", "subset", "n_compared", "n_agree", "pct_agree", "value", "abs_diff", "result"]].to_string())


if __name__ == "__main__":
    main()
