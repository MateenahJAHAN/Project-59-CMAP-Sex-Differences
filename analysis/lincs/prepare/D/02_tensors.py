"""Task D script 2: panel definitions and outcome tensors.

Reads only siginfo, geneinfo and the level-5 gctx files (HTTP range reads).
It reads no line annotation and no sex or state column (cellinfo is not read).

Outputs
  results/D/D_panel_definitions.tsv   one row per panel x candidate compound
  results/D/D_panel_summary.tsv       lines, compounds, signatures, fill per panel
  results/D/D_panel_doses.tsv         numeric pert_dose values behind each label
  results/D/D_panel_lines.tsv         lines per panel with signature counts
  results/D/D_landmarks.tsv           the 978 landmark genes in tensor order
  results/D/D_rangeread_check.tsv     range read vs h5py read for 5 signatures per file
  data/D/tensors/<panel>_resp.npy     float32 (lines, compounds, 978), NaN = missing
  results/D/tensors/<panel>_nsig.tsv  signature counts per line x compound
"""
import os, sys, time, hashlib, json, concurrent.futures as cf
import numpy as np, pandas as pd, requests, fsspec, h5py

ROOT = os.environ["P59_LINCS_WORK"]
RES = os.path.join(ROOT, "results", "D")
TEN = os.path.join(ROOT, "data", "D", "tensors")
os.makedirs(TEN, exist_ok=True); os.makedirs(os.path.join(RES, "tensors"), exist_ok=True)
BASE = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/"
TRT = BASE + "level5_beta_trt_cp_n720216x12328.gctx"
CTL = BASE + "level5_beta_ctl_n58022x12328.gctx"
LOG = os.path.join(ROOT, "logs", "D_runlog.txt")
T0 = time.time()


def log(msg):
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} [02_tensors +{time.time()-T0:.0f}s] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


# ---------------- panels from siginfo ----------------
cols = ["sig_id", "pert_type", "cell_iname", "project_code", "cmap_name", "pert_id",
        "pert_itime", "pert_idose", "pert_dose", "pert_dose_unit", "qc_pass"]
si = pd.read_csv(os.path.join(ROOT, "data/D/cmap/siginfo_beta.txt"), sep="\t", usecols=cols, low_memory=False)
assert si.sig_id.is_unique
log(f"siginfo rows {len(si)}")

PANELS = {
    "P1": dict(project="ASG", time="24 h", gate_labels=["1.11 uM", "0.12 uM"], label="1.11 uM"),
    "P1b": dict(project="ASG", time="24 h", gate_labels=["1.11 uM", "0.12 uM"], label="0.12 uM"),
    "P2": dict(project="CPC", time="6 h", gate_labels=["10 uM"], label="10 uM"),
}
defs, summ, doses, plines, sigsel = [], [], [], [], {}
for pn, P in PANELS.items():
    base = si[(si.project_code == P["project"]) & (si.pert_type == "trt_cp") & (si.pert_itime == P["time"])
              & (si.pert_idose.isin(P["gate_labels"]))]
    lines = sorted(base.cell_iname.unique())
    nl = len(lines)
    cnt = base.groupby("cmap_name").cell_iname.nunique()
    keep = sorted(cnt[cnt >= 0.8 * nl].index)
    gate = base[base.cmap_name.isin(keep)]
    fill_gate = gate.groupby("cmap_name").cell_iname.nunique().sum() / (len(keep) * nl)
    sel = gate[gate.pert_idose == P["label"]].copy()
    fill_lab = sel.groupby("cmap_name").cell_iname.nunique().reindex(keep).fillna(0).sum() / (len(keep) * nl)
    sigsel[pn] = dict(sel=sel, lines=lines, comps=keep)
    pid = base.groupby("cmap_name").pert_id.apply(lambda s: ";".join(sorted(s.unique())))
    nlab = sel.groupby("cmap_name").cell_iname.nunique()
    nsl = sel.groupby("cmap_name").size()
    for c in sorted(cnt.index):
        defs.append([pn, c, pid[c], int(cnt[c]), int(nlab.get(c, 0)), int(nsl.get(c, 0)), c in keep])
    summ.append([pn, P["project"], P["time"], ";".join(P["gate_labels"]), P["label"], nl, len(keep),
                 len(gate), round(fill_gate, 4), len(sel), round(fill_lab, 4)])
    for (lab, d), n in sel.groupby(["pert_idose", "pert_dose"]).size().items():
        doses.append([pn, lab, d, int(n)])
    for l, n in sel.groupby("cell_iname").size().reindex(lines).fillna(0).items():
        q = sel[sel.cell_iname == l].qc_pass
        plines.append([pn, l, int(n)])
    log(f"{pn}: lines {nl} compounds {len(keep)} gate signatures {len(gate)} fill {fill_gate:.4f}; "
        f"label {P['label']} signatures {len(sel)} fill {fill_lab:.4f}")

pd.DataFrame(defs, columns=["panel", "cmap_name", "pert_ids", "n_lines_gate_labels", "n_lines_at_label",
                            "n_signatures_at_label", "included"]).to_csv(os.path.join(RES, "D_panel_definitions.tsv"), sep="\t", index=False)
pd.DataFrame(summ, columns=["panel", "project_code", "pert_itime", "gate_labels", "tensor_label", "n_lines", "n_compounds",
                            "n_signatures_gate", "fill_gate", "n_signatures_label", "fill_label"]).to_csv(os.path.join(RES, "D_panel_summary.tsv"), sep="\t", index=False)
pd.DataFrame(doses, columns=["panel", "pert_idose", "pert_dose", "n_signatures"]).to_csv(os.path.join(RES, "D_panel_doses.tsv"), sep="\t", index=False)
pd.DataFrame(plines, columns=["panel", "cell_iname", "n_signatures_at_label"]).to_csv(os.path.join(RES, "D_panel_lines.tsv"), sep="\t", index=False)

# DMSO pseudo-compound (negative control)
dm = {}
for pn, P in [("P1", PANELS["P1"]), ("P2", PANELS["P2"])]:
    d = si[(si.project_code == P["project"]) & (si.pert_type == "ctl_vehicle") & (si.pert_itime == P["time"])
           & (si.cmap_name == "DMSO") & (si.cell_iname.isin(sigsel[pn]["lines"]))]
    dm[pn] = d
    log(f"{pn} DMSO ctl_vehicle signatures {len(d)} in {d.cell_iname.nunique()} lines")

# ---------------- gctx ids ----------------
fs = fsspec.filesystem("https", block_size=4 * 1024 * 1024)


def read_ids(url):
    with fs.open(url, "rb") as fo, h5py.File(fo, "r") as h:
        col = [x.decode() for x in h["0/META/COL/id"][:]]
        row = [x.decode() for x in h["0/META/ROW/id"][:]]
        ds = h["0/DATA/0/matrix"]
        assert ds.id.get_create_plist().get_layout() == h5py.h5d.CONTIGUOUS
        off = ds.id.get_offset(); shape = ds.shape; dt = ds.dtype
    return col, row, off, shape, dt


ids = {}
for name, url in [("trt", TRT), ("ctl", CTL)]:
    col, row, off, shape, dt = read_ids(url)
    assert len(col) == shape[0] and len(row) == shape[1] and dt == np.float32
    ids[name] = dict(col={c: i for i, c in enumerate(col)}, row=row, off=off, shape=shape, url=url)
    log(f"{name} gctx: shape {shape} offset {off} dtype {dt}")
assert ids["trt"]["row"] == ids["ctl"]["row"]

gi = pd.read_csv(os.path.join(ROOT, "data/D/cmap/geneinfo_beta.txt"), sep="\t", dtype={"gene_id": str})
lm = gi[gi.feature_space == "landmark"]
assert len(lm) == 978, len(lm)
rowpos = {g: i for i, g in enumerate(ids["trt"]["row"])}
lm_pos = np.array([rowpos[g] for g in lm.gene_id])
lm[["gene_id", "gene_symbol"]].assign(gctx_row_index=lm_pos).to_csv(os.path.join(RES, "D_landmarks.tsv"), sep="\t", index=False)
NG = ids["trt"]["shape"][1]
RB = NG * 4


# ---------------- range reads ----------------
def fetch_rows(url, off, idx_sorted):
    """Return dict idx -> landmark vector, reading coalesced runs of rows."""
    runs, s = [], 0
    while s < len(idx_sorted):
        e = s
        while e + 1 < len(idx_sorted) and idx_sorted[e + 1] - idx_sorted[s] < 64 and idx_sorted[e + 1] - idx_sorted[e] <= 4:
            e += 1
        runs.append((idx_sorted[s], idx_sorted[e]))
        s = e + 1
    want = set(idx_sorted)
    sess = {}

    def get(run):
        a, b = run
        import threading
        t = threading.get_ident()
        if t not in sess:
            sess[t] = requests.Session()
        rng = f"bytes={off + a * RB}-{off + (b + 1) * RB - 1}"
        for attempt in range(6):
            try:
                r = sess[t].get(url, headers={"Range": rng}, timeout=120)
                if r.status_code == 206 and len(r.content) == (b - a + 1) * RB:
                    arr = np.frombuffer(r.content, dtype="<f4").reshape(b - a + 1, NG)
                    return {a + i: arr[i, lm_pos].copy() for i in range(b - a + 1) if a + i in want}
                raise IOError(f"status {r.status_code} len {len(r.content)}")
            except Exception as ex:
                time.sleep(2 * (attempt + 1))
        raise IOError(f"failed run {run}")

    out = {}
    with cf.ThreadPoolExecutor(24) as ex:
        for k, d in enumerate(ex.map(get, runs)):
            out.update(d)
            if k % 2000 == 0:
                log(f"  {k}/{len(runs)} runs read")
    return out, len(runs)


def verify(url, off, idxs):
    rows = []
    with fs.open(url, "rb") as fo, h5py.File(fo, "r") as h:
        ds = h["0/DATA/0/matrix"]
        for i in idxs:
            a = ds[i, :][lm_pos]
            r = requests.get(url, headers={"Range": f"bytes={off + i * RB}-{off + (i + 1) * RB - 1}"}, timeout=60)
            b = np.frombuffer(r.content, dtype="<f4")[lm_pos]
            rows.append([url.rsplit("/", 1)[1], int(i), bool(np.array_equal(a, b)), float(np.max(np.abs(a - b)))])
    return rows


need_trt = sorted({ids["trt"]["col"][s] for pn in PANELS for s in sigsel[pn]["sel"].sig_id})
need_ctl = sorted({ids["ctl"]["col"][s] for pn in dm for s in dm[pn].sig_id if s in ids["ctl"]["col"]})
n_missing_ctl = sum(1 for pn in dm for s in dm[pn].sig_id if s not in ids["ctl"]["col"])
log(f"signatures to read: trt {len(need_trt)}, ctl {len(need_ctl)} (DMSO sig_ids not in ctl gctx: {n_missing_ctl})")
rng = np.random.default_rng(20261115)
chk = verify(TRT, ids["trt"]["off"], rng.choice(need_trt, 5, replace=False))
if need_ctl:
    chk += verify(CTL, ids["ctl"]["off"], rng.choice(need_ctl, min(5, len(need_ctl)), replace=False))
pd.DataFrame(chk, columns=["file", "row_index", "identical", "max_abs_diff"]).to_csv(os.path.join(RES, "D_rangeread_check.tsv"), sep="\t", index=False)
assert all(r[2] for r in chk), chk
log("range-read check: all identical to h5py reads")

vals_trt, nr = fetch_rows(TRT, ids["trt"]["off"], need_trt)
log(f"trt rows read {len(vals_trt)} in {nr} requests")
vals_ctl, nr = fetch_rows(CTL, ids["ctl"]["off"], need_ctl) if need_ctl else ({}, 0)
log(f"ctl rows read {len(vals_ctl)} in {nr} requests")


# ---------------- tensors ----------------
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for ch in iter(lambda: f.read(1 << 20), b""):
            h.update(ch)
    return h.hexdigest()


hashes = []
for pn in PANELS:
    sel, lines, comps = sigsel[pn]["sel"], sigsel[pn]["lines"], sigsel[pn]["comps"]
    li = {l: i for i, l in enumerate(lines)}; ci = {c: j for j, c in enumerate(comps)}
    S = np.zeros((len(lines), len(comps), 978), dtype=np.float64)
    N = np.zeros((len(lines), len(comps)), dtype=np.int64)
    for sid, l, c in zip(sel.sig_id, sel.cell_iname, sel.cmap_name):
        v = vals_trt[ids["trt"]["col"][sid]]
        assert np.all(np.isfinite(v))
        S[li[l], ci[c]] += v; N[li[l], ci[c]] += 1
    assert N.sum() == len(sel)
    with np.errstate(invalid="ignore"):
        R = (S / N[:, :, None]).astype(np.float32)
    R[N == 0] = np.nan
    p = os.path.join(TEN, f"{pn}_resp.npy"); np.save(p, R)
    hashes.append([f"data/D/tensors/{pn}_resp.npy", sha(p), os.path.getsize(p), f"{pn} line x compound x landmark mean level-5 z"])
    nsig = pd.DataFrame(N, index=pd.Index(lines, name="cell_iname"), columns=comps)
    q = os.path.join(RES, "tensors", f"{pn}_nsig.tsv"); nsig.to_csv(q, sep="\t")
    hashes.append([f"results/D/tensors/{pn}_nsig.tsv", sha(q), os.path.getsize(q), f"{pn} signature counts (rows lines, columns compounds)"])
    log(f"{pn} tensor {R.shape}, cells filled {(N>0).sum()} of {N.size} ({(N>0).mean():.4f})")
    if pn in dm:
        d = dm[pn][dm[pn].sig_id.isin(ids["ctl"]["col"])]
        D = np.full((len(lines), 978), np.nan, dtype=np.float64); nd = np.zeros(len(lines), dtype=np.int64)
        acc = {}
        for sid, l in zip(d.sig_id, d.cell_iname):
            acc.setdefault(l, []).append(vals_ctl[ids["ctl"]["col"][sid]])
        for l, vs in acc.items():
            D[li[l]] = np.mean(vs, axis=0); nd[li[l]] = len(vs)
        p = os.path.join(TEN, f"{pn}_DMSO_resp.npy"); np.save(p, D.astype(np.float32))
        hashes.append([f"data/D/tensors/{pn}_DMSO_resp.npy", sha(p), os.path.getsize(p), f"{pn} DMSO pseudo-compound line x landmark mean level-5 z"])
        pd.DataFrame({"cell_iname": lines, "n_dmso_signatures": nd}).to_csv(os.path.join(RES, "tensors", f"{pn}_DMSO_nsig.tsv"), sep="\t", index=False)
        log(f"{pn} DMSO pseudo-compound: lines with DMSO {int((nd>0).sum())} of {len(lines)}")
    pd.Series(lines, name="cell_iname").to_csv(os.path.join(RES, "tensors", f"{pn}_lines.tsv"), sep="\t", index=False)
    pd.Series(comps, name="cmap_name").to_csv(os.path.join(RES, "tensors", f"{pn}_compounds.tsv"), sep="\t", index=False)
    for f in [f"{pn}_lines.tsv", f"{pn}_compounds.tsv", f"{pn}_DMSO_nsig.tsv"]:
        q = os.path.join(RES, "tensors", f)
        if os.path.exists(q):
            hashes.append([f"results/D/tensors/{f}", sha(q), os.path.getsize(q), f"{pn} tensor axis / counts"])
pd.DataFrame(hashes, columns=["path", "sha256", "bytes", "description"]).to_csv(os.path.join(RES, "D_tensor_hashes.tsv"), sep="\t", index=False)
log("done")
