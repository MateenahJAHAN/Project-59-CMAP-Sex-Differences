"""Secondary family of task G (draft sections 1 and 4; REC 4d). Every output is labelled secondary.

- TF layer (7n): TF activity = weighted sum of landmark responses over the TF's targets divided by the square root
  of the sum of squared weights; CollecTRI signed regulons (weight = task D's consensus sign; edges without a
  consensus sign have weight 0); TFs with at least 10 landmark targets. Same statistics and permutations as H1
  and H2 (per compound T_c over TFs; panel mean), BH within the readout.
- Focused one-sided test (7e): mean t^2 over the Lopes-Ramos TFs (PANDA-LCL top-30 regulons under 7e A), against
  its permutation null (upper tail) and against random TF sets of the same size drawn from all TFs of the regulon
  table (competitive). MYBL1 has no PANDA-LCL regulon and is reported as not testable.
- X-linked score (X term): per compound, mean t over the X-linked landmarks minus the mean t over all landmarks;
  two-sided permutation p; the panel score is the mean over compounds; competitive p from random landmark sets
  of the same size drawn from all 978 landmarks.
- TF-versus-gene signal-to-null ratio: descriptive (panel statistic / median of its permutation null).
"""
import numpy as np, pandas as pd
import gcore as C

DMSO = "DMSO (pseudo-compound)"


def collectri_weights(reg, genes, min_targets, family):
    n_all = reg.groupby("tf").size()
    elig = n_all[n_all >= min_targets].index
    single = sorted(t for t in elig if "_" not in t)
    cplx = sorted(t for t in elig if "_" in t)
    tfs = {"single": single, "complex": cplx, "all": single + cplx}[family]
    return weights_from_edges(reg[reg.tf.isin(tfs)], genes, tfs, "tf", "target", "sign", n_all)


def weights_from_edges(edges, genes, tfs, tf_col, tg_col, w_col, n_all=None):
    gi = {g: i for i, g in enumerate(genes)}
    W = np.zeros((len(genes), len(tfs)))
    ti = {t: i for i, t in enumerate(tfs)}
    for t, g, w in zip(edges[tf_col], edges[tg_col], edges[w_col]):
        W[gi[g], ti[t]] = float(w)
    nrm = np.sqrt((W ** 2).sum(0))
    keep = nrm > 0
    info = pd.DataFrame(dict(tf=tfs, n_landmark_targets=[int(n_all[t]) if n_all is not None else int((W[:, i] != 0).sum()) for i, t in enumerate(tfs)],
                             n_weighted_targets=(W != 0).sum(0), used=keep))
    return W[:, keep] / nrm[keep], [t for t, k in zip(tfs, keep) if k], info


def tf_override(P, Wn):
    return {j: P.R[:, j, :] @ Wn for j in range(len(P.comps))}


def tf_layer(P, d, Wn, B, h1_sub, with_perm=True):
    df, permT, recs = C.run_term(P, d["Xo"], d["k"], d["mask"], B, Y_override=tf_override(P, Wn), want_perm=with_perm,
                                 lin=d["lin"], rule=d["rule"])
    if not df.tested.any():
        return None
    if with_perm:
        df, obs, pp, null = C.finish_table(df, permT, compounds_only=P.comps)
        sub = h1_sub(df)
        o, p, nl, n = C.panel_stat(df, permT, sub)
        return dict(df=df, panel_stat=o, panel_p=p, n=n, null_median=float(np.median(nl)) if nl is not None else np.nan, recs=recs)
    df = df.copy()
    df["q"] = np.nan; df["rank"] = np.nan
    sub = h1_sub(df.assign(in_ranking=df.tested))
    return dict(df=df, panel_stat=float(df.loc[sub & df.tested.values, "T_obs"].mean()), panel_p=np.nan, n=int((sub & df.tested.values).sum()),
                null_median=np.nan, recs=recs)


def focused(P, d, W_all, tf_all, set_tfs, B, n_rand, seed, h1_sub):
    """One-sided focused test over set_tfs; competitive null from random sets of the same size from tf_all."""
    idx = [tf_all.index(t) for t in set_tfs]
    res = tf_layer(P, d, W_all[:, idx], B, h1_sub)
    if res is None:
        return None
    obs_t = tf_layer(P, d, W_all, 1, h1_sub, with_perm=False)
    T2 = np.array([r["_t"] ** 2 for r in obs_t["recs"] if r["tested"]])           # compounds x all TFs
    names = [r["compound"] for r in obs_t["recs"] if r["tested"]]
    sub = h1_sub(obs_t["df"].assign(in_ranking=obs_t["df"].tested))
    keep = np.array([sub[obs_t["df"].compound.tolist().index(c)] for c in names])
    T2 = T2[keep]
    rng = np.random.default_rng(seed)
    k = len(idx)
    rand = np.array([T2[:, rng.choice(len(tf_all), k, replace=False)].mean() for _ in range(n_rand)])
    obs = float(T2[:, idx].mean())
    res["competitive_obs"] = obs
    res["competitive_p"] = float((1 + np.sum(rand >= obs)) / (n_rand + 1))
    res["n_random_sets"] = n_rand
    return res


def xlinked(P, d, is_x, B, n_rand, seed, h1_compounds):
    nx = int(is_x.sum())
    w = is_x / nx - 1.0 / len(is_x)
    recs = C.run_term_signed(P, d["Xo"], d["k"], d["mask"], B, lambda nm, rec: (lambda t, b: t @ w), lin=d["lin"], rule=d["rule"])
    rows = []
    perm_scores, obs_scores = [], []
    rng = np.random.default_rng(seed)
    rsets = [rng.choice(len(is_x), nx, replace=False) for _ in range(n_rand)]
    Tobs = []
    for r in recs:
        if not r["tested"]:
            rows.append(dict(compound=r["compound"], tested=False, reason=r["reason"])); continue
        sc = float(r["obs"]); null = r["perm"]
        p = float((1 + np.sum(np.abs(null) >= abs(sc) * (1 - 1e-9))) / (B + 1))
        t = r["_t"]
        rnd = np.array([t[s].mean() - t.mean() for s in rsets])
        pc = float((1 + np.sum(np.abs(rnd) >= abs(sc))) / (n_rand + 1))
        rows.append(dict(compound=r["compound"], tested=True, reason="", n_lines=r["n_lines"], score=sc, p_perm_two_sided=p,
                         p_competitive=pc, in_panel_score=r["compound"] in h1_compounds))
        if r["compound"] in h1_compounds:
            perm_scores.append(null); obs_scores.append(sc); Tobs.append(t)
    df = pd.DataFrame(rows)
    if df.tested.any():
        df.loc[df.tested, "q_bh"] = C.bh(df.loc[df.tested, "p_perm_two_sided"].values)
    panel = dict(n_compounds=len(obs_scores), n_xlinked_landmarks=nx)
    if obs_scores:
        ps = np.mean(np.vstack(perm_scores), axis=0)
        ob = float(np.mean(obs_scores))
        Tm = np.vstack(Tobs)
        rnd = np.array([(Tm[:, s].mean(1) - Tm.mean(1)).mean() for s in rsets])
        panel.update(panel_score=ob, panel_p_two_sided=float((1 + np.sum(np.abs(ps) >= abs(ob) * (1 - 1e-9))) / (B + 1)),
                     panel_p_competitive=float((1 + np.sum(np.abs(rnd) >= abs(ob))) / (n_rand + 1)), n_random_sets=n_rand)
    return df, panel
