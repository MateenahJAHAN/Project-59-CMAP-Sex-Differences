"""Inspect the layout of the level-5 gctx files by HTTP range reads (no full download).
Writes results/D/D_gctx_layout.tsv."""
import sys, os, fsspec, h5py, datetime
ROOT = os.environ["P59_LINCS_WORK"]
B = "https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/"
FILES = ["level5_beta_trt_cp_n720216x12328.gctx", "level5_beta_ctl_n58022x12328.gctx"]
rows = []
fs = fsspec.filesystem("https", block_size=2 * 1024 * 1024)
for fn in FILES:
    with fs.open(B + fn, "rb") as fo, h5py.File(fo, "r") as h:
        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                lay = obj.id.get_create_plist().get_layout()
                off = obj.id.get_offset() if lay == h5py.h5d.CONTIGUOUS else None
                rows.append([fn, name, str(obj.shape), str(obj.dtype), str(obj.chunks),
                             {0: "compact", 1: "contiguous", 2: "chunked"}.get(lay, str(lay)), str(off),
                             str(obj.compression)])
                print(rows[-1], flush=True)
        h.visititems(visit)
out = os.path.join(ROOT, "results", "D", "D_gctx_layout.tsv")
with open(out, "w") as f:
    f.write("file\tdataset\tshape\tdtype\tchunks\tlayout\tbyte_offset\tcompression\n")
    for r in rows:
        f.write("\t".join(r) + "\n")
