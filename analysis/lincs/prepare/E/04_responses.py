"""Step 3 of task E: vehicle-referenced responses and line x compound tensors on both scales.

Per treated well and landmark: value minus the mean of the same plate's DMSO wells (plates with
< 2 DMSO wells dropped and counted). Line x compound mean over wells (level 3) and over signatures
(level 5). Each tensor is written with its sha256 (results/E/E_hashes.tsv) before any statistic.

Outputs (data/E/tensors/): R3_<T>.npy, Z5_<T>.npy (lines x compounds x 978, NaN = missing cell),
N3_<T>.npy, N5_<T>.npy (counts), <T>_lines.txt, <T>_compounds.txt, for T in P1, P2, P1_d111, P1_d012.
Also results/E/E_plate_dmso.tsv (DMSO count and drop flag per plate).
"""
import sys, os, time, hashlib
import numpy as np, pandas as pd

ROOT = sys.argv[1]
TD = f"{ROOT}/data/E/tensors"
OUT = f"{ROOT}/results/E"
RUNLOG = f"{ROOT}/logs/E_runlog.txt"
t0 = time.time()


def log(msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\t04_responses\t{msg}"
    print(line, flush=True); open(RUNLOG, "a").write(line + "\n")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


HASHES = pd.read_csv(f"{OUT}/E_hashes.tsv", sep="\t", dtype=str)


def load_checked(name):
    row = HASHES[HASHES.object == name].iloc[0]
    path = f"{ROOT}/{row.local_path}"
    assert sha256(path) == row.sha256, f"sha256 mismatch for {name}"
    return np.load(path) if path.endswith(".npy") else pd.read_csv(path, sep="\t", dtype=str)


def record_hash(name, path, shape, dtype, note):
    global HASHES
    row = pd.DataFrame([dict(object=name, local_path=os.path.relpath(path, ROOT), shape="x".join(map(str, shape)),
                             dtype=str(dtype), bytes=str(os.path.getsize(path)), sha256=sha256(path),
                             written_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), note=note)])
    HASHES = pd.concat([HASHES[HASHES.object != name], row], ignore_index=True)
    HASHES.to_csv(f"{OUT}/E_hashes.tsv", sep="\t", index=False)


wells = pd.read_csv(f"{OUT}/E_wells.tsv.gz", sep="\t", dtype={"sig_id": str, "pert_idose": str})
sigs = pd.read_csv(f"{OUT}/E_panel_signatures.tsv.gz", sep="\t", dtype={"sig_id": str, "pert_idose": str})
plate_rows = []


def cell_tensor(values, line_idx, comp_idx, nL, nC):
    """Mean of rows of `values` per (line, compound); NaN where no row."""
    G = values.shape[1]
    S = np.zeros((nL * nC, G)); N = np.zeros(nL * nC, dtype=np.int64)
    k = line_idx * nC + comp_idx
    np.add.at(S, k, values); np.add.at(N, k, 1)
    with np.errstate(invalid="ignore", divide="ignore"):
        M = S / N[:, None]
    M[N == 0] = np.nan
    return M.reshape(nL, nC, G), N.reshape(nL, nC)


for p in ["P1", "P2"]:
    X = load_checked(f"l3_trt_{p}").astype(np.float64)
    xr = load_checked(f"l3_trt_{p}_rows")
    C = load_checked(f"l3_ctl_{p}")
    cr = load_checked(f"l3_ctl_{p}_rows")
    Zs = load_checked(f"l5_{p}").astype(np.float64)
    zr = load_checked(f"l5_{p}_rows")
    wt = wells[(wells.panel == p) & (wells.role == "treated")].set_index("sample_id")
    wd = wells[(wells.panel == p) & (wells.role == "dmso")].set_index("sample_id")
    assert len(xr) == len(wt) and set(xr.sample_id) == set(wt.index)
    assert len(cr) == len(wd) and set(cr.sample_id) == set(wd.index)
    xm = wt.loc[xr.sample_id].reset_index()
    cm = wd.loc[cr.sample_id].reset_index()
    assert (xm.row_index.astype(str).values == xr.row_index.values).all()

    # plate DMSO means
    dcount = cm.groupby("det_plate").size()
    plates = sorted(set(xm.det_plate))
    keep_pl = [pl for pl in plates if dcount.get(pl, 0) >= 2]
    drop_pl = [pl for pl in plates if dcount.get(pl, 0) < 2]
    for pl in plates:
        plate_rows.append(dict(panel=p, det_plate=pl, n_dmso=int(dcount.get(pl, 0)),
                               n_treated=int((xm.det_plate == pl).sum()), dropped=pl in drop_pl))
    pidx = {pl: i for i, pl in enumerate(keep_pl)}
    cmask = cm.det_plate.isin(pidx).values
    ci = cm.det_plate[cmask].map(pidx).values
    Ssum = np.zeros((len(keep_pl), C.shape[1])); np.add.at(Ssum, ci, C[cmask])
    Dmean = Ssum / np.bincount(ci, minlength=len(keep_pl))[:, None]
    xmask = xm.det_plate.isin(pidx).values
    log(f"{p}: plates {len(plates)}, dropped (<2 DMSO) {len(drop_pl)}, treated wells dropped {int((~xmask).sum())}")
    xm = xm[xmask].reset_index(drop=True)
    r = X[xmask] - Dmean[xm.det_plate.map(pidx).values]

    lines = sorted(set(sigs[sigs.panel == p].cell_iname))
    comps = sorted(set(sigs[sigs.panel == p].cmap_name))
    li = {x: i for i, x in enumerate(lines)}; co = {x: i for i, x in enumerate(comps)}
    zm = sigs[sigs.panel == p].set_index("sig_id").loc[zr.sig_id].reset_index()
    assert len(zm) == len(Zs)

    variants = [(p, np.ones(len(xm), bool), np.ones(len(zm), bool))]
    if p == "P1":
        variants += [("P1_d111", (xm.pert_idose == "1.11 uM").values, (zm.pert_idose == "1.11 uM").values),
                     ("P1_d012", (xm.pert_idose == "0.12 uM").values, (zm.pert_idose == "0.12 uM").values)]
    for T, mw, ms in variants:
        R3, N3 = cell_tensor(r[mw], xm.cell_iname[mw].map(li).values, xm.cmap_name[mw].map(co).values, len(lines), len(comps))
        Z5, N5 = cell_tensor(Zs[ms], zm.cell_iname[ms].map(li).values, zm.cmap_name[ms].map(co).values, len(lines), len(comps))
        assert ((N3 > 0) == (N5 > 0)).all(), f"{T}: observed cells differ between scales"
        open(f"{TD}/{T}_lines.txt", "w").write("\n".join(lines) + "\n")
        open(f"{TD}/{T}_compounds.txt", "w").write("\n".join(comps) + "\n")
        for nm, arr, note in [("R3", R3, "level-3 vehicle-referenced response, line x compound x landmark mean"),
                              ("Z5", Z5, "level-5 z, line x compound x landmark mean"),
                              ("N3", N3, "level-3 treated wells per cell"), ("N5", N5, "level-5 signatures per cell")]:
            path = f"{TD}/{nm}_{T}.npy"
            np.save(path, arr.astype(np.float64) if nm in ("R3", "Z5") else arr)
            record_hash(f"{nm}_{T}", path, arr.shape, arr.dtype, note)
        for nm in ("lines", "compounds"):
            record_hash(f"{T}_{nm}", f"{TD}/{T}_{nm}.txt", (len(lines) if nm == "lines" else len(comps),), "txt", "axis labels")
        log(f"{T}: tensor {R3.shape}, observed cells {int((N3 > 0).sum())} of {N3.size}, wells {int(N3.sum())}, signatures {int(N5.sum())}")
    # well-level responses for NC2 (treated) and DMSO wells with plate index, hashed
    np.save(f"{TD}/r3_wells_{p}.npy", r)
    record_hash(f"r3_wells_{p}", f"{TD}/r3_wells_{p}.npy", r.shape, r.dtype, "level-3 treated well responses (kept plates)")
    xm[["sample_id", "det_plate", "cell_iname", "cmap_name", "pert_idose"]].to_csv(f"{TD}/r3_wells_{p}_rows.tsv", sep="\t", index=False)
    record_hash(f"r3_wells_{p}_rows", f"{TD}/r3_wells_{p}_rows.tsv", (len(xm), 5), "tsv", "row labels")

pd.DataFrame(plate_rows).to_csv(f"{OUT}/E_plate_dmso.tsv", sep="\t", index=False)
log(f"done in {time.time()-t0:.0f}s")
