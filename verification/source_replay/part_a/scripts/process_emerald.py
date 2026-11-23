"""Step 2: one streaming pass over the 116 EmeraldBay shards.

For every shard: download (or reuse a verified local copy), check SHA-256 against the HF LFS oid,
then per cell compute depth (gene-metadata tokens with non-negative values only) and the UMI counts
of XIST, the 8 Y genes and the 8 X gametologs.

Outputs (under data/F1/emerald_proc/, gitignored):
  cells_<shard>.parquet : vehicle cells only (DMSO_TF and DMSO_T0), one row per cell, marker counts.
  acc_w<k>.npz          : per worker, summed gene counts of deep cells (depth >= 2500) over ALL conditions,
                          by line x group, without any condition label (identity check, plan step 7d).
                          groups: 0 = Y-negative (0 Y UMIs), 1..20 = Y-positive random bins,
                                  21 = XIST-negative, 22..41 = XIST-positive random bins.
  shardlog_<shard>.json : per-shard checks (rows, special-token entries, negative values).
The raw shard is deleted after processing, and the deletion is logged.
"""
import glob
import sys, os, json, time, datetime, fcntl, traceback
import numpy as np, pandas as pd, pyarrow.parquet as pq, scipy.sparse as sp
import multiprocessing as mp
sys.path.insert(0, os.path.dirname(__file__))
from fetch import download, log_download
from common import MARKERS, SEED, DEEP_MIN, NBINS, NGROUP

ROOT = sys.argv[1]
NWORK = int(sys.argv[2]) if len(sys.argv) > 2 else 4
BASE = "https://huggingface.co/datasets/tahoebio/EmeraldBay/resolve/752eb8f0d68765fe554194b6c003dfb782ad661d/"
RAW = os.path.join(ROOT, "data/F1/emerald")
OUT = os.path.join(ROOT, "data/F1/emerald_proc")
DONE = os.path.join(ROOT, "data/F1/emerald_done")
VEH = {"DMSO_TF", "DMSO_T0"}


def log_deletion(path, nbytes, reason):
    p = os.path.join(ROOT, "provenance/F1_deletions.tsv")
    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(p, "a") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        fh.write(f"F1\temerald_pass\t{os.path.relpath(path, ROOT)}\t{nbytes}\t{ts}\t{reason}\n")
        fcntl.flock(fh, fcntl.LOCK_UN)


def setup():
    g = pd.read_parquet(os.path.join(ROOT, "data/F1/meta/gene_metadata.parquet"))
    assert g.token_id.is_unique and g.gene_symbol.notna().all()
    vocab = json.load(open(os.path.join(ROOT, "data/F1/meta/extended_vocab.json")))
    ntok = max(vocab.values()) + 1
    valid = np.zeros(ntok, bool); valid[g.token_id.values] = True
    midx = np.full(ntok, -1, np.int64)
    for i, s in enumerate(MARKERS):
        row = g[g.gene_symbol == s]
        assert len(row) == 1, s
        midx[int(row.token_id.iloc[0])] = i
    cl = pd.read_parquet(os.path.join(ROOT, "data/F1/meta/cell_line_metadata.parquet"))
    lines = sorted(cl.Cell_ID_Cellosaur.unique())
    assert len(lines) == 52
    return valid, midx, ntok, lines


def process_shard(path, shard_i, valid, midx, ntok, line_ix, acc, cnt):
    pf = pq.ParquetFile(path)
    vrows = []
    chk = dict(rows=0, special_entries=0, negative_entries=0, nonint_entries=0, rowgroups=pf.metadata.num_row_groups)
    for rg in range(pf.metadata.num_row_groups):
        t = pf.read_row_group(rg)
        meta = t.select(["drug", "drugname_drugconc", "cell_line", "sample", "BARCODE_SUB_LIB_ID"]).to_pandas()
        n = len(meta); chk["rows"] += n
        ga = t.column("genes").combine_chunks(); ea = t.column("expressions").combine_chunks()
        off = ga.offsets.to_numpy(); assert np.array_equal(off, ea.offsets.to_numpy())
        tok = ga.values.to_numpy(); ex = ea.values.to_numpy()
        if off[0] != 0:  # sliced arrays: rebase
            tok = tok[off[0]:off[-1]]; ex = ex[off[0]:off[-1]]; off = off - off[0]
        assert len(tok) == off[-1] and tok.min() >= 0 and tok.max() < ntok
        cell_of = np.repeat(np.arange(n), np.diff(off))
        ok = valid[tok] & (ex >= 0)
        chk["special_entries"] += int((~valid[tok]).sum()); chk["negative_entries"] += int((ex < 0).sum())
        chk["nonint_entries"] += int((ex[ok] != np.round(ex[ok])).sum())
        depth = np.bincount(cell_of[ok], weights=ex[ok], minlength=n)
        ngenes = np.bincount(cell_of[ok], minlength=n)
        mk = midx[tok]; m = ok & (mk >= 0)
        mc = np.bincount(cell_of[m] * len(MARKERS) + mk[m], weights=ex[m], minlength=n * len(MARKERS)).reshape(n, len(MARKERS))
        li = meta.cell_line.map(line_ix)
        assert li.notna().all(), "unknown cell line"
        li = li.values.astype(np.int64)
        ysum = mc[:, 1:9].sum(1); xist = mc[:, 0]
        # vehicle per-cell table
        isv = meta.drug.isin(VEH).values
        if isv.any():
            d = meta.loc[isv, ["BARCODE_SUB_LIB_ID", "cell_line", "drug", "sample"]].copy()
            d["shard"] = shard_i; d["depth"] = depth[isv]; d["n_genes"] = ngenes[isv]
            for j, s in enumerate(MARKERS):
                d[s] = mc[isv, j].astype(np.int32)
            vrows.append(d)
        # identity accumulators: deep cells, all conditions, no condition label kept
        rng = np.random.default_rng([SEED, shard_i, rg])
        by = rng.integers(0, NBINS, n); bx = rng.integers(0, NBINS, n)
        deep = depth >= DEEP_MIN
        gy = np.where(ysum == 0, 0, 1 + by); gx = np.where(xist == 0, 1 + NBINS, 2 + NBINS + bx)
        di = np.where(deep)[0]
        if len(di):
            X = sp.csr_matrix((np.where(ok, ex, 0.0).astype(np.float32), tok, off), shape=(n, ntok))[di]
            slots = np.concatenate([li[di] * NGROUP + gy[di], li[di] * NGROUP + gx[di]])
            cols = np.concatenate([np.arange(len(di)), np.arange(len(di))])
            M = sp.csr_matrix((np.ones(len(slots), np.float32), (slots, cols)), shape=(acc.shape[0], len(di)))
            P = (M @ X).tocoo()
            acc[P.row, P.col] += P.data
            np.add.at(cnt, slots, 1)
        del t, ga, ea, tok, ex
    v = pd.concat(vrows, ignore_index=True) if vrows else None
    return v, chk


def save_checkpoint(path, acc, cnt, done):
    temp = path + ".tmp.npz"
    np.savez(temp, acc=acc, cnt=cnt, done=np.array(done))
    os.replace(temp, path)
    for name in done:
        open(os.path.join(DONE, name + ".done"), "a").close()


def worker(wid, q, valid, midx, ntok, lines, tree):
    line_ix = {l: i for i, l in enumerate(lines)}
    accp = os.path.join(OUT, f"acc_w{wid}.npz")
    if os.path.exists(accp):
        z = np.load(accp); acc = z["acc"].astype(np.float32); cnt = z["cnt"]; done = list(z["done"])
    else:
        acc = np.zeros((len(lines) * NGROUP, ntok), np.float32); cnt = np.zeros(len(lines) * NGROUP, np.int64); done = []
    since = 0
    while True:
        item = q.get()
        if item is None:
            break
        i, e = item
        name = os.path.basename(e["path"])
        if os.path.exists(os.path.join(DONE, name + ".done")):
            continue
        t0 = time.time()
        local = os.path.join(RAW, name)
        nb, sha, how = download(BASE + e["path"], local, e["lfs"]["oid"], e["lfs"]["size"])
        if how == "downloaded":
            log_download(ROOT, BASE + e["path"], local, nb, sha, "sha256 matches HF LFS oid; deleted after processing")
        metap = os.path.join(OUT, f"meta_{i:05d}.parquet")
        cols = ["drug", "drugname_drugconc", "cell_line", "sample", "BARCODE_SUB_LIB_ID"]
        pq.read_table(local, columns=cols).to_pandas().assign(shard=name).to_parquet(metap, index=False)
        t1 = time.time()
        v, chk = process_shard(local, i, valid, midx, ntok, line_ix, acc, cnt)
        if v is not None:
            v.to_parquet(os.path.join(OUT, f"cells_{i:05d}.parquet"), index=False)
        chk.update(shard=name, sha256=sha, bytes=nb, how=how, n_vehicle=0 if v is None else len(v),
                   t_download_s=round(t1 - t0, 1), t_process_s=round(time.time() - t1, 1), worker=wid)
        json.dump(chk, open(os.path.join(OUT, f"shardlog_{i:05d}.json"), "w"))
        done.append(name); since += 1
        if since >= 6:
            save_checkpoint(accp, acc, cnt, done)
            since = 0
        os.remove(local); log_deletion(local, nb, "raw shard; per-cell vehicle table and deep-cell profile sums saved")
        print(f"w{wid} {name} rows={chk['rows']} veh={chk['n_vehicle']} dl={chk['t_download_s']}s proc={chk['t_process_s']}s", flush=True)
    save_checkpoint(accp, acc, cnt, done)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True); os.makedirs(DONE, exist_ok=True)
    valid, midx, ntok, lines = setup()
    tree = sorted(json.load(open(os.path.join(ROOT, "data/F1/meta/tree_expression.json"))), key=lambda e: e["path"])
    assert len(tree) == 116
    # local (already verified) shards first, then the rest
    order = sorted(range(len(tree)), key=lambda i: (not os.path.exists(os.path.join(RAW, os.path.basename(tree[i]["path"]))), i))
    # A crash can occur after committing the accumulator but before writing its completion markers.
    # Rebuild the queue from the committed arrays, which are the source of truth.
    committed = set()
    for prior in sorted(glob.glob(os.path.join(OUT, "acc_w*.npz"))):
        with np.load(prior, allow_pickle=False) as z:
            names = set(z["done"].tolist())
        if committed & names:
            raise ValueError("Duplicate shards in existing checkpoints")
        committed.update(names)
    for marker in os.listdir(DONE):
        if marker.removesuffix(".done") not in committed:
            raise ValueError(f"Completion marker without accumulator: {marker}")
    for name in committed:
        open(os.path.join(DONE, name + ".done"), "a").close()
    q = mp.Queue()
    for i in order:
        q.put((i, tree[i]))
    for _ in range(NWORK):
        q.put(None)
    ps = [mp.Process(target=worker, args=(w, q, valid, midx, ntok, lines, tree)) for w in range(NWORK)]
    for p in ps: p.start()
    for p in ps: p.join()
    bad = [p.exitcode for p in ps if p.exitcode != 0]
    assert not bad, f"worker failures {bad}"
    ndone = len(os.listdir(DONE))
    print("shards done", ndone)
    assert ndone == 116
