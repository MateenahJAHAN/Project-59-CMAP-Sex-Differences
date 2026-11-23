"""A3: depth-matched shares of marker-negative cells across cultures (vehicle cells only).

Rule, written and committed before this script ran:
  depth_matching.md (historical commits 10f74ec and 7ac5e72).
Estimators: the local Tahoe caller's functions (partA_local/code/mosaic_models.py), unchanged.
Pairs whose input table is missing are skipped and reported as pending.
"""
import sys, json, hashlib, time
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tahoe"))
from common import Y8, SEED, load_cells                      # noqa: E402
from mosaic_models import fit_nb_mix, fit_pois_mix, deep_zero  # noqa: E402

T0 = time.time()
ROOT = Path(sys.argv[1]).resolve()
F1 = ROOT / "results/F1"
LOC = ROOT / "results/tahoe/latest"
SCI = ROOT / "results/F1/sciplex_vehicle_cells.tsv.gz"
OUT = ROOT / "results/depth_matching"
OUT.mkdir(parents=True, exist_ok=True)
RULE = Path(__file__).with_name("depth_matching.md")
TOL, NREP, MIN_DEEP = 0.15, 10, 100
PAIRS = [  # (pair index, line, Cellosaurus ID, marker, dataset, role) as in the rule's table
    (0, "A549", "CVCL_0023", "Y", "EmeraldBay", "test"),
    (1, "A549", "CVCL_0023", "Y", "sciPlex2", "test"),
    (2, "A549", "CVCL_0023", "Y", "sciPlex3", "test"),
    (3, "A549", "CVCL_0023", "Y", "sciPlex4", "test"),
    (4, "SW48", "CVCL_1724", "XIST", "EmeraldBay", "test"),
    (5, "NCI-H460", "CVCL_0459", "Y", "EmeraldBay", "context"),
    (6, "LoVo", "CVCL_0399", "Y", "EmeraldBay", "context"),
]
ESTS = ["a3000", "a8000", "b2000", "c"]
LOG = []


def log(s):
    s = time.strftime("%H:%M:%S") + " " + s
    print(s, flush=True)
    LOG.append(s)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def as_bool(s):
    m = s.map({True: True, False: False, "True": True, "False": False})
    assert m.notna().all(), s[m.isna()]
    return m.astype(bool)


def as_int(a):
    a = np.asarray(a, float)
    assert np.isfinite(a).all() and (a == np.round(a)).all() and (a >= 0).all()
    return a.astype(np.int64)


# ---------------- Tahoe (local caller's cells, filter and per-plate background)
S = pd.read_csv(LOC / "partA_local_states.tsv", sep="\t").set_index("cellosaurus_id")
bgt = pd.read_csv(LOC / "partA_local_plate_background.tsv", sep="\t")
bg_y_t = bgt[bgt.line_set.str.startswith("Y-negative lines with both")].set_index("unit")["per_umi"]
bg_x_t = bgt[bgt.marker == "XIST"].set_index("unit")["per_umi"]
assert len(bg_y_t) == len(bg_x_t) == 14
cells = load_cells(["drug", "cell_line_id", "cell_name", "plate", "total_umi", "XIST", "y8"])


def tahoe_side(cid, marker):
    r = S.loc[cid]
    c = cells[cells.cell_line_id == cid]
    keep = np.ones(len(c), bool)
    if not r.y_high and marker != "Y":
        keep &= (c.y8 < 3).to_numpy()
    if not r.xist_high and marker != "XIST":
        keep &= (c.XIST < 3).to_numpy()
    c = c[keep]
    assert len(c) == r.n_vehicle_cells, (cid, len(c), r.n_vehicle_cells)
    bg = c.plate.map(bg_y_t if marker == "Y" else bg_x_t).to_numpy(float)
    assert np.isfinite(bg).all() and (bg > 0).all()
    return dict(d=as_int(c.total_umi), y=as_int(c.y8 if marker == "Y" else c.XIST), bg=bg)


# ---------------- other datasets: same filter rule and background definition (rule, "Data")
def background(c, st, unitcol):
    """per-UMI background per unit: marker UMIs / UMIs over vehicle cells with < 3 marker UMIs of marker-negative lines."""
    ylines = st.index[as_bool(st.y_negative_line) & (st.arm_state == "Y lost")]
    xlines = st.index[as_bool(st.x_negative_line)]
    cy = c[c.cid.isin(ylines) & (c.y8 < 3)]
    cx = c[c.cid.isin(xlines) & (c.XIST < 3)]
    by = cy.groupby(unitcol).y8.sum() / cy.groupby(unitcol).depth.sum()
    bx = cx.groupby(unitcol).XIST.sum() / cx.groupby(unitcol).depth.sum()
    return by, bx, sorted(ylines), sorted(xlines)


OTHER, BG, BGINFO = {}, {}, {}
emc = pd.read_csv(F1 / "F1_emerald_vehicle_cells.tsv.gz", sep="\t")
emc = emc[emc.drug == "DMSO_TF"].copy()
emc["cid"] = emc.line
emc["y8"] = emc[Y8].sum(axis=1)
ems = pd.read_csv(F1 / "F1_emerald_states.tsv", sep="\t").set_index("cellosaurus_id")
assert ems.index.is_unique
OTHER["EmeraldBay"] = (emc, ems, "sample")
by, bx, yl, xl = background(emc, ems, "sample")
BG["EmeraldBay"] = (by, bx)
BGINFO["EmeraldBay"] = dict(y_lines=len(yl), xist_lines=len(xl), y_bg_per_1000=(1000 * by).round(5).to_dict(),
                            xist_bg_per_1000=(1000 * bx).round(5).to_dict())

if not SCI.is_file():
    raise FileNotFoundError(SCI)
SCI_OK = True
if SCI_OK:
    spc = pd.read_csv(SCI, sep="\t")
    sps = pd.read_csv(F1 / "F1_sciplex_states.tsv", sep="\t")
    name2cid = dict(zip(sps.cell_name, sps.cellosaurus_id))
    spc["cid"] = spc.line.map(lambda v: v if str(v).startswith("CVCL_") else name2cid[v])
    spc["y8"] = spc[Y8].sum(axis=1)
    for ds in ["sciPlex2", "sciPlex3", "sciPlex4"]:
        c = spc[spc.dataset == ds].copy()
        st = sps[sps["sample"] == ds].set_index("cellosaurus_id")
        assert st.index.is_unique and len(c) > 0
        OTHER[ds] = (c, st, "bg_batch")
        if ds != "sciPlex2":
            by, bx, yl, xl = background(c, st, "bg_batch")
            BG[ds] = (by, bx)
            BGINFO[ds] = dict(y_lines=yl, y_bg_per_1000=(1000 * by).round(5).to_dict())
    # sciPlex2 has no Y-negative line: pooled sciPlex3 rate (addendum 1)
    c3, st3, _ = OTHER["sciPlex3"]
    yl3 = st3.index[as_bool(st3.y_negative_line) & (st3.arm_state == "Y lost")]
    cy = c3[c3.cid.isin(yl3) & (c3.y8 < 3)]
    pooled = cy.y8.sum() / cy.depth.sum()
    c2 = OTHER["sciPlex2"][0]
    BG["sciPlex2"] = (pd.Series(pooled, index=sorted(c2.bg_batch.unique())), pd.Series(dtype=float))
    BGINFO["sciPlex2"] = dict(y_bg_per_1000_pooled_sciPlex3=round(1000 * pooled, 5))
else:
    log("sci-Plex per-cell table not present yet: sci-Plex pairs pending")


def other_side(ds, cid, marker):
    c, st, unitcol = OTHER[ds]
    r = st.loc[cid]
    c = c[c.cid == cid]
    keep = np.ones(len(c), bool)
    if marker == "Y" and as_bool(pd.Series([r.x_negative_line])).iloc[0]:
        keep &= (c.XIST < 3).to_numpy()
    if marker == "XIST" and as_bool(pd.Series([r.y_negative_line])).iloc[0]:
        keep &= (c.y8 < 3).to_numpy()
    n_before = len(c)
    c = c[keep]
    by, bx = BG[ds]
    bg = c[unitcol].map(by if marker == "Y" else bx).to_numpy(float)
    assert np.isfinite(bg).all() and (bg > 0).all(), (ds, cid, marker)
    return dict(d=as_int(c.depth), y=as_int(c.y8 if marker == "Y" else c.XIST), bg=bg, n_before_filter=n_before)


# ---------------- estimators and thinning
def estimates(s, rng):
    d, y, bg = s["d"].astype(float), s["y"].astype(float), s["bg"]
    pos = y >= 1
    r_pos = y[pos].sum() / d[pos].sum()
    fb = fit_nb_mix(y, d, bg)
    nb = (fb["r_high"], fb["k_high"])
    a3, a8 = deep_zero(y, d, 3000, r_pos, nb, rng), deep_zero(y, d, 8000, r_pos, nb, rng)
    fc = fit_pois_mix(y, d, bg)
    return dict(n_cells=len(d), mean_depth=float(d.mean()), median_depth=float(np.median(d)),
                a3000=a3["zero_frac"], n_ge3000=a3["n_deep"], a8000=a8["zero_frac"], n_ge8000=a8["n_deep"],
                b2000=fb["pi_2000"], b_converged=fb["converged"], b_at_bound=fb["at_bound"],
                c=fc["pi"], c_converged=fc["converged"])


def thin(s, target, rng):
    d, y = s["d"], s["y"]
    n = len(d)
    order = np.argsort(d, kind="mergesort")
    q = (np.arange(1, n + 1) - 0.5) / n
    t = np.empty(n)
    t[order] = np.quantile(np.asarray(target, float), q, method="linear")
    p = np.minimum(1.0, t / d)
    y2 = rng.binomial(y, p)
    o2 = rng.binomial(d - y, p)
    d2 = y2 + o2
    assert (d2 > 0).all()
    return dict(d=d2, y=y2, bg=s["bg"])


rows, reps, unmatched = [], [], []
done_unmatched = set()
for k, name, cid, marker, ds, role in PAIRS:
    if ds not in OTHER:
        rows.append(dict(pair=k, line=name, marker=marker, dataset=ds, role=role, estimator="all", status="pending (no input table)"))
        continue
    tah, oth = tahoe_side(cid, marker), other_side(ds, cid, marker)
    for side_name, s in (("Tahoe", tah), (ds, oth)):
        key = (side_name, name, marker)
        if key not in done_unmatched:
            e = estimates(s, np.random.default_rng(SEED + 100 * k + 99))
            unmatched.append(dict(dataset=side_name, line=name, marker=marker, **e))
            done_unmatched.add(key)
    deeper = "Tahoe" if tah["d"].mean() > oth["d"].mean() else ds
    thick, shallow = (tah, oth) if deeper == "Tahoe" else (oth, tah)
    e_sh = estimates(shallow, np.random.default_rng(SEED + 100 * k + 99))
    e_th = []
    for r in range(NREP):
        rng = np.random.default_rng(SEED + 100 * k + r)
        tt = thin(thick, shallow["d"], rng)
        e = estimates(tt, rng)
        e_th.append(e)
        reps.append(dict(pair=k, line=name, marker=marker, dataset=ds, thinned_side=deeper, replicate=r, **e))
    E = pd.DataFrame(e_th)
    log(f"pair {k} {name} {marker} Tahoe vs {ds}: thinned {deeper}; mean depth {thick['d'].mean():.0f} -> "
        f"{E.mean_depth.median():.0f} (target {shallow['d'].mean():.0f}); cells {len(thick['d'])} vs {len(shallow['d'])}")
    for est in ESTS:
        med = float(E[est].median())
        if est in ("a3000", "a8000"):
            nn = "n_ge" + est[1:]
            used = bool(E[nn].min() >= MIN_DEEP and e_sh[nn] >= MIN_DEEP)
            why = "" if used else f"fewer than {MIN_DEEP} cells above {est[1:]} UMIs (thinned min {int(E[nn].min())}, other {int(e_sh[nn])})"
        elif est == "b2000":
            okt = (E.b_converged & ~E.b_at_bound)
            used = bool(okt.all() and e_sh["b_converged"] and not e_sh["b_at_bound"])
            why = "" if used else f"fit not converged or at a bound ({int(okt.sum())} of {NREP} thinned replicates ok; other side ok={e_sh['b_converged'] and not e_sh['b_at_bound']})"
        else:
            used, why = True, ""
        tah_v = med if deeper == "Tahoe" else e_sh[est]
        oth_v = e_sh[est] if deeper == "Tahoe" else med
        diff = abs(tah_v - oth_v) if used else np.nan
        rows.append(dict(pair=k, line=name, marker=marker, dataset=ds, role=role, estimator=est, status="computed",
                         thinned_side=deeper, share_tahoe=tah_v, share_other=oth_v, abs_diff=diff,
                         thinned_rep_min=float(E[est].min()), thinned_rep_max=float(E[est].max()),
                         used=used, not_used_reason=why, agree=(diff <= TOL) if used else np.nan,
                         n_cells_tahoe=len(tah["d"]), n_cells_other=len(oth["d"]),
                         mean_depth_tahoe_before=float(tah["d"].mean()), mean_depth_other_before=float(oth["d"].mean()),
                         mean_depth_thinned_after=float(E.mean_depth.median()),
                         c_converged_all=bool(E.c_converged.all() and e_sh["c_converged"])))

R = pd.DataFrame(rows)
P = pd.DataFrame(reps)
U = pd.DataFrame(unmatched)
R.to_csv(OUT / "A3_pairs.tsv", sep="\t", index=False, float_format="%.6f")
P.to_csv(OUT / "A3_replicates.tsv", sep="\t", index=False, float_format="%.6f")
U.to_csv(OUT / "A3_unmatched.tsv", sep="\t", index=False, float_format="%.6f")

# ---------------- decisions (rule, "Decision rule")
dec = {}
for k, name, cid, marker, ds, role in PAIRS:
    sub = R[(R.pair == k) & (R.status == "computed") & (R.used == True)]  # noqa: E712
    if (R[(R.pair == k)].status == "computed").any():
        dec[f"{k}:{name}:{ds}"] = dict(role=role, n_estimators_used=int(len(sub)),
                                       agree=(bool(sub.agree.all()) if len(sub) else None),
                                       max_abs_diff=(float(sub.abs_diff.max()) if len(sub) else None))
    else:
        dec[f"{k}:{name}:{ds}"] = dict(role=role, agree=None, status="pending")
a549_em = dec["0:A549:EmeraldBay"]["agree"]
sp = [dec[f"{k}:A549:sciPlex{j}"]["agree"] for k, j in ((1, 2), (2, 3), (3, 4))]
sp_known = [v for v in sp if v is not None]
sciplex_culture = (sum(bool(v) for v in sp) >= 2) if len(sp_known) == 3 else None
summary = dict(rule_file=str(RULE), rule_sha256=sha(RULE), seed_base=SEED, n_replicates=NREP, tolerance=TOL,
               min_deep_cells=MIN_DEEP, inputs_sha256={
                   "cells_markers.parquet": sha(ROOT / "tahoe/dmso/cells_markers.parquet"),
                   "partA_local_states.tsv": sha(LOC / "partA_local_states.tsv"),
                   "partA_local_plate_background.tsv": sha(LOC / "partA_local_plate_background.tsv"),
                   "F1_emerald_vehicle_cells.tsv.gz": sha(F1 / "F1_emerald_vehicle_cells.tsv.gz"),
                   "F1_emerald_states.tsv": sha(F1 / "F1_emerald_states.tsv"),
                   **({"sciplex_vehicle_cells.tsv.gz": sha(SCI), "F1_sciplex_states.tsv": sha(F1 / "F1_sciplex_states.tsv")}
                      if SCI_OK else {})},
               background=BGINFO, pairs=dec,
               a549_emerald_agrees=a549_em, a549_sciplex_experiments_agree=sp,
               a549_sciplex_culture_agrees=sciplex_culture,
               a549_A3_pass=(bool(a549_em and sciplex_culture) if sciplex_culture is not None and a549_em is not None else None),
               a549_k_of_3=(1 + int(bool(a549_em)) + int(bool(sciplex_culture))) if sciplex_culture is not None else None,
               sw48_A3_pass=dec["4:SW48:EmeraldBay"]["agree"],
               runtime_s=round(time.time() - T0, 1))
(OUT / "A3_summary.json").write_text(json.dumps(summary, indent=1, default=str))
show = R[R.status == "computed"][["pair", "line", "dataset", "estimator", "thinned_side", "share_tahoe", "share_other",
                                  "abs_diff", "used", "agree", "not_used_reason"]]
log("\n" + show.to_string(index=False))
log(json.dumps({k: v for k, v in summary.items() if k not in ("inputs_sha256", "background")}, indent=1, default=str))
(OUT / "A3_log.txt").write_text("\n".join(LOG) + "\n")
