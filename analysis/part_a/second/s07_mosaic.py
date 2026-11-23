"""Step 6: mosaic tests (a) to (f) for every eligible line with at least 1,000 vehicle cells.
Y is tested in lines whose final line Y score is above the Y cut-off, XIST in lines whose final XIST is
above the XIST cut-off. Vehicle cells only."""
import json, sys, time
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.special import gammaln, expit, logsumexp
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import runlog, SEED, Y8, YP, YQ, X8, MARKERS
import sparse_io

t0 = time.time()
R, D = "results/F2", "data/F2"
a = pd.read_parquet(f"{D}/cells_analysis.parquet")
st = pd.read_csv(f"{R}/F2_states.tsv", sep="\t").set_index("cell_name")
cut = pd.read_csv(f"{R}/F2_cutoffs.tsv", sep="\t")
cf = cut[cut.cutoff_set == "final_line_x_plate"].set_index("score").cutoff
amb = pd.read_csv(f"{R}/F2_ambient.tsv", sep="\t")
BG = {"Y": dict(zip(amb.plate, amb.y_per_1000_umis / 1000)), "XIST": dict(zip(amb.plate, amb.xist_per_1000_umis / 1000))}
a["X8"] = a[X8].sum(axis=1)
a["plate_n"] = a.plate.str.replace("plate", "").astype(int)
DEPTHS = [1000, 2000, 4000, 8000]
NB_BOOT = 1000


def pois_logpmf(y, lam):
    return y * np.log(lam) - lam - gammaln(y + 1)


def nb_logpmf(y, mean, k):
    return gammaln(y + k) - gammaln(k) - gammaln(y + 1) + k * np.log(k / (k + mean)) + y * np.log(mean / (k + mean))


def fit_poisson_mix(y, d, bg):
    """(c) pi Poisson(bg d) + (1-pi) Poisson(lam d), bg fixed per cell; EM."""
    lam = max(y.sum() / d.sum(), 1e-6) * 1.5
    pi_ = 0.3
    ll_old = -np.inf
    lo = pois_logpmf(y, bg * d)
    for it in range(2000):
        hi = pois_logpmf(y, lam * d)
        l1, l2 = np.log(pi_) + lo, np.log1p(-pi_) + hi
        ll = np.logaddexp(l1, l2)
        r = np.exp(l1 - ll)
        tot = ll.sum()
        pi_ = min(max(r.mean(), 1e-9), 1 - 1e-9)
        lam = ((1 - r) * y).sum() / ((1 - r) * d).sum()
        if abs(tot - ll_old) < 1e-8:
            break
        ll_old = tot
    return dict(pi=pi_, lam=lam, ll=tot, iters=it + 1)


def fit_nb_mix(y, d, bg, start_pi=0.3):
    """(b) low Poisson(bg d), high NB(mean mu d, size k), low share logistic(alpha + beta (ln d - ln 2000))."""
    ld = np.log(d) - np.log(2000)
    lo = pois_logpmf(y, bg * d)
    pos = y > 0
    mu0 = y[pos].sum() / d[pos].sum() if pos.any() else 1e-4

    def nll(th):
        al, be, lmu, lk = th
        p = expit(al + be * ld)
        hi = nb_logpmf(y, np.exp(lmu) * d, np.exp(lk))
        return -np.logaddexp(np.log(np.clip(p, 1e-300, 1)) + lo, np.log(np.clip(1 - p, 1e-300, 1)) + hi).sum()

    best = None
    sp = min(max(start_pi, 0.02), 0.98)
    for th0 in ([np.log(sp / (1 - sp)), 0.0, np.log(mu0), np.log(5.0)], [0.0, 0.0, np.log(mu0), np.log(1.0)], [-2.0, 0.0, np.log(mu0), np.log(20.0)]):
        res = minimize(nll, th0, method="L-BFGS-B", bounds=[(-15, 15), (-10, 10), (np.log(1e-7), np.log(1.0)), (np.log(1e-3), np.log(1e4))])
        if best is None or res.fun < best.fun:
            best = res
    al, be, lmu, lk = best.x
    out = dict(alpha=al, beta=be, mu=np.exp(lmu), k=np.exp(lk), nll=best.fun, converged=bool(best.success))
    for dd in DEPTHS:
        out[f"share_{dd}"] = float(expit(al + be * (np.log(dd) - np.log(2000))))
    return out


def mh_or(A, B, strata):
    num = den = 0.0
    for s in np.unique(strata):
        m = strata == s
        a_ = np.sum(A[m] & B[m]); b_ = np.sum(A[m] & ~B[m]); c_ = np.sum(~A[m] & B[m]); d_ = np.sum(~A[m] & ~B[m])
        if min(a_, b_, c_, d_) == 0:
            a_, b_, c_, d_ = a_ + .5, b_ + .5, c_ + .5, d_ + .5
        n = a_ + b_ + c_ + d_
        num += a_ * d_ / n; den += b_ * c_ / n
    return num / den


# ---------- identity gene set and centroids ----------
PB = np.load(f"{D}/pb_line_plate.npy")
ck = json.load(open(f"{D}/reduce_ckpt.json"))
LINES, PLATES = ck["lines"], ck["plates"]
gm = pd.read_parquet(f"{D}/meta/gene_metadata.parquet")
sym = gm.gene_symbol.astype(str)
excl = sym.isin(Y8 + ["TXLNGY", "PRKY", "XIST", "TSIX"]) | sym.str.startswith("MT-") | sym.str.startswith("TTTY")
tot = PB.sum(axis=0)
cpm_all = 1e6 * tot / tot.sum()
GSET = np.flatnonzero((~excl.to_numpy()) & (cpm_all >= 1))
LPB = PB.reshape(len(LINES), len(PLATES), -1).sum(axis=1)[:, GSET].astype(float)
LPROF = np.log2(1e6 * LPB / LPB.sum(axis=1, keepdims=True) + 1)
var = LPROF.var(axis=0)
TOP = np.argsort(-var)[:2000]
CENT = LPROF[:, TOP]
json.dump(dict(n_gene_set=int(len(GSET)), n_excluded=int(excl.sum()), n_top_var=len(TOP)), open(f"{R}/F2_identity_geneset.json", "w"))


def prof(counts):
    c = np.asarray(counts, float)
    return np.log2(1e6 * c / c.sum() + 1)


def corr(x, y):
    return float(np.corrcoef(x, y)[0, 1])


def identity(line, cl, rng):
    li = LINES.index(line)
    sp = sparse_io.load_line(li)  # barcodes, csr (cells x genes)
    pos_in = pd.Index(sp["barcodes"]).get_indexer(cl.barcode)
    deep = cl.umis.to_numpy() >= 2500
    neg = cl.marker.to_numpy() == 0
    have = pos_in >= 0
    assert have[deep].all(), "deep cell missing from sparse store"
    out = dict(n_deep=int(deep.sum()), n_deep_neg=int((deep & neg).sum()), n_deep_pos=int((deep & ~neg).sum()))
    if out["n_deep_neg"] < 20:
        out["f_status"] = "not_testable"
        return out
    X = sp["csr"][:, GSET]
    dn = np.flatnonzero(deep & neg); dp = np.flatnonzero(deep & ~neg)
    depth = cl.umis.to_numpy()
    edges = np.quantile(depth[deep], np.linspace(0, 1, 21))
    bins = np.clip(np.searchsorted(edges, depth, side="right") - 1, 0, 19)
    P1, P2, repl = [], [], 0
    for b in range(20):
        nb = dn[bins[dn] == b]; pb = dp[bins[dp] == b]
        if len(nb) == 0:
            continue
        if len(pb) >= 2 * len(nb):
            pick = rng.choice(pb, 2 * len(nb), replace=False)
        else:
            pick = rng.choice(pb, 2 * len(nb), replace=True) if len(pb) else np.array([], int); repl += 1
        P1.extend(pick[:len(nb)]); P2.extend(pick[len(nb):])
    rows = lambda idx: np.asarray(X[pos_in[np.asarray(idx)]].sum(axis=0)).ravel()
    pn, p1, p2 = prof(rows(dn)), prof(rows(P1)), prof(rows(P2))
    r_neg, r_split = corr(pn, p1), corr(p1, p2)
    others = {LINES[j]: corr(pn, LPROF[j]) for j in range(len(LINES)) if j != li}
    best_other = max(others, key=others.get)
    out.update(f_status="tested", bins_with_replacement=repl, r_neg_vs_pos=r_neg, r_split_half_pos=r_split,
               r_neg_vs_own_line_pseudobulk=corr(pn, LPROF[li]), best_other_line=best_other, r_neg_vs_best_other=others[best_other],
               f_holds=bool(r_neg >= r_split - 0.01 and r_neg >= others[best_other] + 0.10))
    # single deep cells: nearest line centroid over the 2,000 most variable genes
    def nearest(idx):
        M = X[pos_in[idx]][:, TOP].toarray().astype(float)
        totg = np.asarray(X[pos_in[idx]].sum(axis=1)).ravel()
        Lc = np.log2(1e6 * M / totg[:, None] + 1)
        Lc = (Lc - Lc.mean(1, keepdims=True)) / (Lc.std(1, keepdims=True) + 1e-12)
        C = (CENT - CENT.mean(1, keepdims=True)) / CENT.std(1, keepdims=True)
        sc = Lc @ C.T / Lc.shape[1]
        return np.argmax(sc, axis=1)
    nn = nearest(dn)
    ps = rng.choice(dp, min(len(dn), len(dp)), replace=False)
    nnp = nearest(ps)
    out.update(single_neg_share_own_line=float((nn == li).mean()), single_pos_share_own_line=float((nnp == li).mean()),
               single_neg_most_common_other=LINES[np.bincount(nn[nn != li], minlength=len(LINES)).argmax()] if (nn != li).any() else "")
    return out


ncell = a.groupby("line").size()
JOBS = []
for line in sorted(ncell.index):
    if ncell[line] < 1000 or not st.loc[line, "eligible"]:
        continue
    if st.loc[line, "y_score"] > cf["y_score"]:
        JOBS.append((line, "Y"))
    if st.loc[line, "xist"] > cf["xist"]:
        JOBS.append((line, "XIST"))


def run_job(job):
    line, mk = job
    rng = np.random.default_rng([SEED, LINES.index(line), 0 if mk == "Y" else 1])
    curves = []
    if True:
        g = a[a.line == line].copy().reset_index(drop=True)
        g["marker"] = g.Y8 if mk == "Y" else g.XIST
        y = g.marker.to_numpy().astype(float); d = g.umis.to_numpy().astype(float)
        bg = g.plate.map(BG[mk]).to_numpy()
        r = dict(cell_name=line, marker=mk, state=st.loc[line, "state"], n_cells=len(g))
        # (a)
        pos = y > 0
        rate = y[pos].sum() / d[pos].sum()
        for thr in (3000, 8000):
            m = d >= thr
            n = int(m.sum()); z = int((y[m] == 0).sum())
            r[f"a_n_cells_ge{thr}"] = n
            if n:
                f0 = z / n
                bs = rng.binomial(n, f0, NB_BOOT) / n
                r[f"a_zero_frac_ge{thr}"] = f0
                r[f"a_ci_low_ge{thr}"], r[f"a_ci_high_ge{thr}"] = np.percentile(bs, [2.5, 97.5])
                r[f"a_expected_zero_frac_ge{thr}"] = float(np.exp(-rate * d[m]).mean())
        r["a_marker_rate_positive_cells"] = rate
        # (c)
        c = fit_poisson_mix(y, d, bg)
        r.update(c_share=c["pi"], c_lambda=c["lam"], c_iters=c["iters"])
        # (b)
        b = fit_nb_mix(y, d, bg, start_pi=c["pi"])
        r.update({f"b_{k}": v for k, v in b.items()})
        # (d)
        for par, mpar in [("odd", g.plate_n.to_numpy() % 2 == 1), ("even", g.plate_n.to_numpy() % 2 == 0)]:
            bb = fit_nb_mix(y[mpar], d[mpar], bg[mpar], start_pi=c["pi"])
            cc = fit_poisson_mix(y[mpar], d[mpar], bg[mpar])
            r[f"d_b_share2000_{par}"] = bb["share_2000"]; r[f"d_c_share_{par}"] = cc["pi"]; r[f"d_n_{par}"] = int(mpar.sum())
        r["d_abs_diff_b2000"] = abs(r["d_b_share2000_odd"] - r["d_b_share2000_even"])
        # (e)
        if mk == "Y":
            strata = pd.qcut(d, 10, labels=False, duplicates="drop")
            A = g.Yq.to_numpy() == 0; B = g.Yp.to_numpy() == 0
            h1 = rng.binomial(g.X8.to_numpy(), 0.5); h2 = g.X8.to_numpy() - h1
            r["e_or_yq0_yp0"] = mh_or(A, B, strata)
            r["e_or_xhalves_control"] = mh_or(h1 == 0, h2 == 0, strata)
            r["e_ratio"] = r["e_or_yq0_yp0"] / r["e_or_xhalves_control"]
        # (f)
        r.update({f"f_{k}" if not k.startswith("f_") else k: v for k, v in identity(line, g, rng).items()})
        # rule
        rb = 0.15 <= r["b_share_2000"] <= 0.85
        rd = r["d_abs_diff_b2000"] <= 0.10
        re_ = (r["e_ratio"] > 2) if mk == "Y" else True
        rf = bool(r.get("f_holds", False))
        r.update(rule_b=rb, rule_d=rd, rule_e=re_ if mk == "Y" else np.nan, rule_f=rf, flagged_b_d_e=bool(rb and rd and re_), mosaic=bool(rb and rd and re_ and rf))
        sh = [r.get("a_zero_frac_ge3000"), r.get("a_zero_frac_ge8000")] + [r[f"b_share_{dd}"] for dd in DEPTHS] + [r["c_share"]]
        sh = [s for s in sh if s is not None and not np.isnan(s)]
        r["share_range_min"], r["share_range_max"] = min(sh), max(sh)
        # depth curve
        edges = [0, 1000, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000, 1e9]
        for lo_, hi_ in zip(edges[:-1], edges[1:]):
            m = (d >= lo_) & (d < hi_)
            if m.sum() == 0:
                continue
            dm = np.median(d[m])
            pi_b = float(expit(b["alpha"] + b["beta"] * (np.log(d[m]) - np.log(2000))).mean())
            curves.append(dict(cell_name=line, marker=mk, depth_lo=lo_, depth_hi=hi_, n_cells=int(m.sum()), median_umis=dm,
                               zero_frac=float((y[m] == 0).mean()), expected_zero_frac_homogeneous=float(np.exp(-rate * d[m]).mean()),
                               b_model_low_share=pi_b,
                               b_model_zero_frac=float((pi_b * np.exp(-bg[m] * d[m]) + (1 - pi_b) * np.exp(nb_logpmf(0, b["mu"] * d[m], b["k"]))).mean())))
        print(line, mk, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k in
              ["a_zero_frac_ge3000", "a_zero_frac_ge8000", "b_share_2000", "c_share", "d_abs_diff_b2000", "e_ratio", "f_r_neg_vs_pos", "f_r_split_half_pos", "f_r_neg_vs_best_other", "mosaic"]}, flush=True)
    return r, curves


if __name__ == "__main__":
    # Keep each fit in one process. Forking an initialized numerical runtime can deadlock on macOS.
    res = [run_job(job) for job in JOBS]
    rows = [r for r, _ in res]
    curves = [c for _, cs in res for c in cs]
    M = pd.DataFrame(rows)
    M.to_csv(f"{R}/F2_mosaic.tsv", sep="\t", index=False)
    pd.DataFrame(curves).to_csv(f"{R}/F2_depth_curve.tsv", sep="\t", index=False)
    runlog("s07_mosaic", t0, f"tests {len(M)}")
