"""Read row (gene) and column (instance/signature) ids of each gctx via HTTP range reads.
Writes data/E/ids/<key>_row_ids.txt and <key>_col_ids.txt."""
import sys, os, time
import fsspec, h5py, numpy as np
URLS = {
 "l3_trt": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level3/level3_beta_trt_cp_n1805898x12328.gctx",
 "l3_ctl": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level3/level3_beta_ctl_n188708x12328.gctx",
 "l5_trt": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/level5_beta_trt_cp_n720216x12328.gctx",
}
outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True)
fs = fsspec.filesystem("https")
for k, u in URLS.items():
    t = time.time()
    with fs.open(u, "rb", block_size=8*2**20, cache_type="readahead") as f:
        h = h5py.File(f, "r")
        rows = [x.decode() for x in h["0/META/ROW/id"][:]]
        cols = [x.decode() for x in h["0/META/COL/id"][:]]
        shp = h["0/DATA/0/matrix"].shape
    assert len(rows) == shp[1] and len(cols) == shp[0], (k, len(rows), len(cols), shp)
    assert len(set(rows)) == len(rows) and len(set(cols)) == len(cols)
    open(f"{outdir}/{k}_row_ids.txt", "w").write("\n".join(rows) + "\n")
    open(f"{outdir}/{k}_col_ids.txt", "w").write("\n".join(cols) + "\n")
    print(k, len(rows), len(cols), rows[:3], cols[:2], f"{time.time()-t:.1f}s", flush=True)
