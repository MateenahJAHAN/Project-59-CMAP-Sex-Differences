"""Step 2: read and reduce, vehicle row groups only.

For every row group listed in the index with vehicle rows: one HTTP range read of the row group's
bytes into memory, decode genes, expressions, drug, sample, BARCODE_SUB_LIB_ID, cell_line_id and plate,
keep rows with drug == 'DMSO_TF', drop the treated rows at once, free the raw bytes and log it.
Treated rows are never written and no statistic is computed on them.

Outputs (vehicle cells only):
  data/F2/cells_parts/shard_XXXXX.parquet  per-cell marker counts
  data/F2/pb_line_plate.npy                all-gene pseudobulk per line x plate (analysis set)
  data/F2/sparse/<line>.{genes,counts}.u16, .nnz.u32, .barcodes.txt  sparse per-cell counts
  logs/F2_reduce_rowgroups.tsv              bytes fetched and freed per row group
Resumable through data/F2/reduce_ckpt.json.
"""
import os, sys, time, json, threading, glob
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd
import pyarrow.parquet as pq, pyarrow.compute as pc
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hf_io
from common import runlog, MARKERS, NAMED_SPARSE_LINES, DEEP_SPARSE

D = "data/F2"
os.makedirs(f"{D}/cells_parts", exist_ok=True)
os.makedirs(f"{D}/sparse", exist_ok=True)
CK = f"{D}/reduce_ckpt.json"
RGLOG = "logs/F2_reduce_rowgroups.tsv"

gm = pd.read_parquet(f"{D}/meta/gene_metadata.parquet")
G = len(gm)
assert gm.token_id.is_unique and gm.gene_symbol.notna().all()
LUT = np.full(int(gm.token_id.max()) + 1, -1, dtype=np.int64)
LUT[gm.token_id.to_numpy()] = np.arange(G)
sym = gm.gene_symbol.to_numpy()
midx = {}
for m in MARKERS:
    hit = np.flatnonzero(sym == m)
    assert len(hit) == 1, (m, len(hit))
    midx[m] = int(hit[0])
MK = np.full(G, -1, dtype=np.int64)
for j, m in enumerate(MARKERS):
    MK[midx[m]] = j
MT = np.zeros(G, dtype=bool)
MT[np.flatnonzero(pd.Series(sym).str.startswith("MT-").to_numpy())] = True

obs = pd.read_parquet(f"{D}/obs_dmso.parquet")
full = obs[obs.pass_filter == "full"].copy()
assert full.BARCODE_SUB_LIB_ID.is_unique
LINES = sorted(full.cell_name.unique())
PLATES = sorted(full.plate.unique(), key=lambda p: int(p.replace("plate", "")))
assert len(LINES) == 50 and len(PLATES) == 14
L2I = {l: i for i, l in enumerate(LINES)}
P2I = {p: i for i, p in enumerate(PLATES)}
cl2name = full[["cell_line", "cell_name"]].drop_duplicates()
assert cl2name.cell_line.is_unique and cl2name.cell_name.is_unique
CVCL2I = {c: L2I[n] for c, n in zip(cl2name.cell_line, cl2name.cell_name)}
full_idx = pd.Index(full.BARCODE_SUB_LIB_ID.to_numpy())
full_line = full.cell_name.map(L2I).to_numpy()
full_plate = full.plate.to_numpy()
full_sample = full["sample"].to_numpy()
NAMED = {L2I[l] for l in NAMED_SPARSE_LINES}
NU = len(LINES) * len(PLATES)

veh = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{D}/index_parts/veh_*.parquet"))], ignore_index=True)
fl = pd.read_csv("results/F2/F2_filelist.tsv", sep="\t")
size_of = dict(zip(fl.path, fl["size"]))
COLS = ["genes", "expressions", "drug", "sample", "BARCODE_SUB_LIB_ID", "cell_line_id", "plate"]
lock = threading.Lock()


def process_rowgroup(pf, f, shard, rg, vexp):
    md = pf.metadata.row_group(rg)
    names = [pf.metadata.schema.column(j).path for j in range(pf.metadata.num_columns)]
    starts, ends = [], []
    for j in range(md.num_columns):
        c = md.column(j)
        if names[j].split(".")[0] in COLS:
            s = c.dictionary_page_offset if c.dictionary_page_offset is not None else c.data_page_offset
            starts.append(s); ends.append(s + c.total_compressed_size)
    f.prefetch([(min(starts), max(ends) - min(starts))])
    fetched = max(ends) - min(starts)
    t = pf.read_row_group(rg, columns=COLS)
    t = t.filter(pc.equal(t.column("drug"), "DMSO_TF"))  # treated rows dropped here
    freed = f.drop()  # raw bytes freed (row-group bytes and footer copy)
    n = t.num_rows
    bc = t.column("BARCODE_SUB_LIB_ID").to_numpy(zero_copy_only=False)
    assert n == len(vexp), (shard, rg, n, len(vexp))
    assert (bc == vexp.BARCODE_SUB_LIB_ID.to_numpy()).all(), (shard, rg)
    gl = t.column("genes").combine_chunks()
    el = t.column("expressions").combine_chunks()
    off = gl.offsets.to_numpy()
    assert (off == el.offsets.to_numpy()).all()
    tok = gl.values.to_numpy()[off[0]:off[-1]]
    ex = el.values.to_numpy()[off[0]:off[-1]]
    cell = np.repeat(np.arange(n), np.diff(off))
    gi = np.where(tok < len(LUT), LUT[np.minimum(tok, len(LUT) - 1)], -1)
    nongene = gi < 0
    zero = (~nongene) & (ex == 0)
    keep = (~nongene) & (ex != 0)
    nongene_tokens = np.unique(tok[nongene])
    cell_k, gi_k, ex_k = cell[keep], gi[keep], ex[keep]
    assert (ex_k > 0).all() and (ex_k == np.round(ex_k)).all(), (shard, rg)
    cnt = ex_k.astype(np.int64)
    key = cell_k * G + gi_k
    assert len(np.unique(key)) == len(key), ("duplicate gene in a cell", shard, rg)
    umis = np.bincount(cell_k, weights=cnt, minlength=n).astype(np.int64)
    ngenes = np.bincount(cell_k, minlength=n)
    mito = np.bincount(cell_k[MT[gi_k]], weights=cnt[MT[gi_k]], minlength=n).astype(np.int64)
    me = MK[gi_k]; s = me >= 0
    M = np.zeros((n, len(MARKERS)), dtype=np.int64)
    np.add.at(M, (cell_k[s], me[s]), cnt[s])
    # analysis-set membership
    pos = full_idx.get_indexer(bc)
    inset = pos >= 0
    line_shard = np.array([CVCL2I.get(x, -1) for x in t.column("cell_line_id").to_numpy(zero_copy_only=False)])
    plate_shard = t.column("plate").to_numpy(zero_copy_only=False)
    sample_shard = t.column("sample").to_numpy(zero_copy_only=False)
    pi = pos[inset]
    assert (full_line[pi] == line_shard[inset]).all(), ("line mismatch", shard, rg)
    assert (full_plate[pi] == plate_shard[inset]).all(), ("plate mismatch", shard, rg)
    assert (full_sample[pi] == sample_shard[inset]).all(), ("sample mismatch", shard, rg)
    cells = pd.DataFrame(dict(shard=shard, row_group=rg, barcode=bc, cell_line_id=t.column("cell_line_id").to_numpy(zero_copy_only=False),
                              plate=plate_shard, sample=sample_shard, in_analysis_set=inset, umis=umis, n_genes=ngenes,
                              mito_umis=mito, n_nongene_entries=np.bincount(cell[nongene], minlength=n),
                              n_zero_entries=np.bincount(cell[zero], minlength=n)))
    for j, m in enumerate(MARKERS):
        cells[m] = M[:, j]
    # pseudobulk (analysis set only)
    unit = np.full(n, -1)
    unit[inset] = line_shard[inset] * len(PLATES) + np.array([P2I[p] for p in plate_shard[inset]])
    ue = unit[cell_k]; se = ue >= 0
    uu, uinv = np.unique(ue[se], return_inverse=True)
    pb = np.bincount(uinv * G + gi_k[se], weights=cnt[se], minlength=len(uu) * G).reshape(len(uu), G).astype(np.int32)
    # sparse store: analysis set, named lines (all cells) or UMIs >= DEEP_SPARSE
    store = inset & ((umis >= DEEP_SPARSE) | np.isin(line_shard, list(NAMED)))
    sp = {}
    if store.any():
        assert cnt.max() < 65536
        st_e = store[cell_k]
        for li in np.unique(line_shard[store]):
            cs = np.flatnonzero(store & (line_shard == li))
            em = st_e & np.isin(cell_k, cs)
            sp[int(li)] = (bc[cs], ngenes[cs].astype(np.uint32), gi_k[em].astype(np.uint16), cnt[em].astype(np.uint16))
    del t, gl, el, tok, ex
    rgrow = dict(shard=shard, row_group=rg, bytes_fetched=fetched, bytes_freed=freed, vehicle_rows=n,
                 zero_entries=int(zero.sum()), nongene_entries=int(nongene.sum()),
                 nongene_tokens=",".join(map(str, nongene_tokens)))
    return cells, (uu, pb), sp, rgrow


def process_shard(shard):
    v = veh[veh.shard == shard]
    path = f"data/train-{shard:05d}-of-03388.parquet"
    f = hf_io.RangeFile(path, int(size_of[path]), tail=1 << 16)
    pf = pq.ParquetFile(f)
    tail_blocks = list(f.blocks)
    out_cells, out_pb, out_sp, rgrows = [], [], [], []
    for rg, vexp in v.groupby("row_group", sort=True):
        f.blocks = list(tail_blocks)
        c, p, s, r = process_rowgroup(pf, f, shard, int(rg), vexp.sort_values("row_in_group"))
        out_cells.append(c); out_pb.append(p); out_sp.append(s); rgrows.append(r)
    f.blocks = []
    return shard, pd.concat(out_cells, ignore_index=True), out_pb, out_sp, rgrows


def main():
    t0 = time.time()
    shards = sorted(veh.shard.unique().tolist())
    sp_files = {}
    if os.path.exists(CK):
        ck = json.load(open(CK))
        done = set(ck["done"])
        PB = np.load(f"{D}/" + ck["pb_snapshot"])
        with open(RGLOG, "ab") as stream:
            stream.truncate(ck["log_size"])
        for li in range(len(LINES)):
            for ext, sz in ck["sizes"][str(li)].items():
                fn = f"{D}/sparse/{li:02d}.{ext}"
                with open(fn, "ab") as g:
                    g.truncate(sz)
        print("resume from checkpoint:", len(done), "shards done", flush=True)
    else:
        done = set()
        PB = np.zeros((NU, G), dtype=np.int64)
        if any(os.path.getsize(f) for f in glob.glob(f"{D}/sparse/*")):
            raise RuntimeError("Incomplete sparse files without a checkpoint. Use a fresh work directory.")
        with open(RGLOG, "w") as g:
            g.write("shard\trow_group\tbytes_fetched\tbytes_freed\tvehicle_rows\tzero_entries\tnongene_entries\tnongene_tokens\n")
    for li in range(len(LINES)):
        sp_files[li] = {ext: open(f"{D}/sparse/{li:02d}.{ext}", "ab") for ext in ["genes.u16", "counts.u16", "nnz.u32", "barcodes.txt"]}
    todo = [s for s in shards if s not in done]
    print(f"{len(todo)} shards to do", flush=True)

    def ckpt():
        for li in sp_files:
            for fh in sp_files[li].values():
                fh.flush(); os.fsync(fh.fileno())
        snapshot = f"pb_checkpoint_{len(done)}.npy"
        np.save(f"{D}/pb_line_plate.tmp.npy", PB)
        os.replace(f"{D}/pb_line_plate.tmp.npy", f"{D}/{snapshot}")
        sizes = {str(li): {ext: os.path.getsize(f"{D}/sparse/{li:02d}.{ext}") for ext in sp_files[li]} for li in sp_files}
        with open(CK + ".tmp", "w") as g:
            json.dump(dict(done=sorted(done), sizes=sizes, lines=LINES, plates=PLATES,
                           pb_snapshot=snapshot, log_size=os.path.getsize(RGLOG)), g)
        os.replace(CK + ".tmp", CK)
        for prior in glob.glob(f"{D}/pb_checkpoint_*.npy"):
            if os.path.basename(prior) != snapshot:
                os.remove(prior)
        completed = f"{D}/pb_line_plate.npy"
        next_link = completed + ".link"
        if os.path.lexists(next_link):
            os.unlink(next_link)
        os.symlink(snapshot, next_link)
        os.replace(next_link, completed)

    if not os.path.exists(CK):
        ckpt()
    n_done = 0
    with ThreadPoolExecutor(max_workers=2) as ex:
        def ordered_shards():
            for start in range(0, len(todo), 4):
                yield from ex.map(process_shard, todo[start:start + 4])
        it = ordered_shards()
        for shard, cells, pbs, sps, rgrows in it:
            cells.to_parquet(f"{D}/cells_parts/shard_{shard:05d}.parquet", index=False)
            for uu, pb in pbs:
                PB[uu] += pb.astype(np.int64)
            for sp in sps:
                for li, (bcs, nnz, gg, cc) in sp.items():
                    fh = sp_files[li]
                    fh["barcodes.txt"].write(("\n".join(bcs) + "\n").encode())
                    fh["nnz.u32"].write(nnz.tobytes()); fh["genes.u16"].write(gg.tobytes()); fh["counts.u16"].write(cc.tobytes())
            with open(RGLOG, "a") as g:
                for r in rgrows:
                    g.write("\t".join(str(r[k]) for k in ["shard", "row_group", "bytes_fetched", "bytes_freed", "vehicle_rows", "zero_entries", "nongene_entries", "nongene_tokens"]) + "\n")
            done.add(shard); n_done += 1
            if n_done % 100 == 0:
                ckpt()
                print(f"{n_done}/{len(todo)} shards, {time.time()-t0:.0f}s, http {hf_io.STATS}", flush=True)
    ckpt()
    runlog("s03_reduce", t0, f"shards {len(todo)}; http {hf_io.STATS}")
    print("done", time.time() - t0, hf_io.STATS, flush=True)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test":
        s = int(sys.argv[2])
        t1 = time.time()
        shard, cells, pbs, sps, rgrows = process_shard(s)
        print(time.time() - t1, len(cells), cells.in_analysis_set.sum(), rgrows)
        print(cells[["umis", "n_genes", "XIST", "RPS4Y1"]].describe().to_string())
    else:
        main()
