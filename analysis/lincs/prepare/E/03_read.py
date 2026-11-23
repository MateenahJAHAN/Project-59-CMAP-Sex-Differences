"""Step 2 of task E: read the 978 landmark values of every panel well (level 3) and every panel
signature (level 5) by HTTP range requests. Check 5 rows of each file against an h5py read.
Write each matrix to disk and record its sha256 in results/E/E_hashes.tsv before any statistic.

Outputs (data/E/tensors/, not committed):
  l3_trt_<P>.npy (float32), l3_ctl_<P>.npy (float64), l5_<P>.npy (float32), plus row-order tables
  results/E/E_rangeread_check.tsv, results/E/E_hashes.tsv
"""
import sys, os, time, hashlib
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rangeread as rr

ROOT = sys.argv[1]
TD = f"{ROOT}/data/E/tensors"; os.makedirs(TD, exist_ok=True)
OUT = f"{ROOT}/results/E"
RUNLOG = f"{ROOT}/logs/E_runlog.txt"
URLS = {
    "l3_trt": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level3/level3_beta_trt_cp_n1805898x12328.gctx",
    "l3_ctl": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level3/level3_beta_ctl_n188708x12328.gctx",
    "l5_trt": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/level5_beta_trt_cp_n720216x12328.gctx",
}
rng = np.random.default_rng(20261115)


def log(msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\t03_read\t{msg}"
    print(line, flush=True)
    open(RUNLOG, "a").write(line + "\n")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def record_hash(name, path, shape, dtype, note):
    hp = f"{OUT}/E_hashes.tsv"
    row = pd.DataFrame([dict(object=name, local_path=os.path.relpath(path, ROOT), shape="x".join(map(str, shape)),
                             dtype=str(dtype), bytes=os.path.getsize(path), sha256=sha256(path),
                             written_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), note=note)])
    if os.path.exists(hp):
        old = pd.read_csv(hp, sep="\t", dtype=str)
        old = old[old.object != name]
        row = pd.concat([old, row.astype(str)], ignore_index=True)
    row.to_csv(hp, sep="\t", index=False)


# landmark columns
gene_ids = [x.strip() for x in open(f"{ROOT}/data/E/ids/l3_trt_row_ids.txt")]
gi = pd.read_csv(f"{ROOT}/data/E/meta/geneinfo_beta.txt", sep="\t", usecols=["gene_id", "gene_symbol", "feature_space"], dtype=str)
lm = gi[gi.feature_space == "landmark"]
assert len(lm) == 978
pos = {g: i for i, g in enumerate(gene_ids)}
lm = lm.assign(col=lm.gene_id.map(pos)).sort_values("col")
assert lm.col.notna().all()
lm_cols = lm.col.astype(int).values
lm[["gene_id", "gene_symbol", "col"]].to_csv(f"{OUT}/E_landmarks.tsv", sep="\t", index=False)

lays = {k: rr.layout(u) for k, u in URLS.items()}
for k, l in lays.items():
    log(f"layout {k}: offset {l['offset']} shape {l['shape']} dtype {l['dtype']}")
assert lays["l3_trt"]["dtype"] == np.dtype("<f4") and lays["l3_ctl"]["dtype"] == np.dtype("<f8") and lays["l5_trt"]["dtype"] == np.dtype("<f4")

# equality check: range read vs h5py, 5 rows per file
wells = pd.read_csv(f"{OUT}/E_wells.tsv.gz", sep="\t", dtype={"sig_id": str, "pert_idose": str})
sigs = pd.read_csv(f"{OUT}/E_panel_signatures.tsv.gz", sep="\t", dtype={"sig_id": str, "pert_idose": str})
l5ids = pd.Series(np.arange(lays["l5_trt"]["shape"][0]), index=[x.strip() for x in open(f"{ROOT}/data/E/ids/l5_trt_col_ids.txt")])
sigs["row_index"] = sigs.sig_id.map(l5ids)
assert sigs.row_index.notna().all()
sigs["row_index"] = sigs.row_index.astype(int)
checks = []
for key, rows in [("l3_trt", wells[wells.file == "l3_trt"].row_index.values),
                  ("l3_ctl", wells[wells.file == "l3_ctl"].row_index.values),
                  ("l5_trt", sigs.row_index.values)]:
    pick = rng.choice(rows, 5, replace=False)
    a, _ = rr.read_rows(URLS[key], lays[key], pick)
    b = rr.h5py_rows(URLS[key], pick)
    for j, r in enumerate(pick):
        eq = bool(np.array_equal(a[j], b[j], equal_nan=True))
        checks.append(dict(file=key, row_index=int(r), n_values=a.shape[1], identical=eq,
                           max_abs_diff=float(np.nanmax(np.abs(a[j].astype(float) - b[j].astype(float))))))
ck = pd.DataFrame(checks)
ck.to_csv(f"{OUT}/E_rangeread_check.tsv", sep="\t", index=False)
log(f"range vs h5py check: {int(ck.identical.sum())}/{len(ck)} identical")
assert ck.identical.all()

total_bytes = 0
for p in ["P1", "P2"]:
    for key in ["l3_ctl", "l3_trt"]:
        w = wells[(wells.panel == p) & (wells.file == key)].sort_values("row_index").reset_index(drop=True)
        X, nb = rr.read_rows(URLS[key], lays[key], w.row_index.values, cols=lm_cols, log=log, label=f"{p} {key}")
        total_bytes += nb
        assert np.isfinite(X).all(), f"non-finite values in {p} {key}"
        path = f"{TD}/{key}_{p}.npy"
        np.save(path, X)
        w[["sample_id", "row_index"]].to_csv(f"{TD}/{key}_{p}_rows.tsv", sep="\t", index=False)
        record_hash(f"{key}_{p}", path, X.shape, X.dtype, f"level-3 wells x 978 landmarks, rows ordered as {key}_{p}_rows.tsv")
        record_hash(f"{key}_{p}_rows", f"{TD}/{key}_{p}_rows.tsv", (len(w), 2), "tsv", "row order")
    s = sigs[sigs.panel == p].sort_values("row_index").reset_index(drop=True)
    Z, nb = rr.read_rows(URLS["l5_trt"], lays["l5_trt"], s.row_index.values, cols=lm_cols, log=log, label=f"{p} l5")
    total_bytes += nb
    nonfin = int((~np.isfinite(Z)).sum())
    log(f"{p} l5 non-finite values: {nonfin}")
    path = f"{TD}/l5_{p}.npy"
    np.save(path, Z)
    s[["sig_id", "row_index"]].to_csv(f"{TD}/l5_{p}_rows.tsv", sep="\t", index=False)
    record_hash(f"l5_{p}", path, Z.shape, Z.dtype, f"level-5 signatures x 978 landmarks, rows ordered as l5_{p}_rows.tsv")
    record_hash(f"l5_{p}_rows", f"{TD}/l5_{p}_rows.tsv", (len(s), 2), "tsv", "row order")
log(f"total bytes fetched {total_bytes}")
pd.DataFrame([dict(total_bytes_fetched=total_bytes)]).to_csv(f"{OUT}/E_bytes_fetched.tsv", sep="\t", index=False)
