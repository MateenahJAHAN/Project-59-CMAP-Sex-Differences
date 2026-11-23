"""Inspect LINCS2020 gctx files over HTTP range reads (no full download).
Prints dataset paths, shapes, dtypes, chunking and the matrix byte offset."""
import sys, json, time
import fsspec, h5py, numpy as np

URLS = {
 "l3_trt": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level3/level3_beta_trt_cp_n1805898x12328.gctx",
 "l3_ctl": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level3/level3_beta_ctl_n188708x12328.gctx",
 "l5_trt": "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/level5_beta_trt_cp_n720216x12328.gctx",
}
out = {}
fs = fsspec.filesystem("https", block_size=2**20)
for k, u in URLS.items():
    t = time.time()
    f = fs.open(u, "rb", block_size=2**20, cache_type="bytes")
    h = h5py.File(f, "r")
    info = {}
    def visit(name, obj):
        if isinstance(obj, h5py.Dataset):
            info[name] = dict(shape=list(obj.shape), dtype=str(obj.dtype), chunks=obj.chunks,
                              compression=obj.compression, offset=obj.id.get_offset())
    h.visititems(visit)
    out[k] = info
    print(k, json.dumps(info, indent=1), f"{time.time()-t:.1f}s", flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
