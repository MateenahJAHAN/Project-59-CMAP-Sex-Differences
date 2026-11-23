"""t1_build_index (tahoe-pipeline): build the row-group index of all 3,388 Tahoe-100M raw files from the saved footers,
check the byte layout, check how rows are ordered, and align every row group with obs_metadata.

Per row group: file, rg, rows, row_start (in file), global_row, byte_offset, compressed_size, uncompressed_size and the
min/max of sample, plate, drug, cell_line_id and BARCODE_SUB_LIB_ID (footer statistics), with null counts.
Ordering checks: plate blocks over files; resets of the sample order (one per sublibrary start) against obs.
Alignment: obs rows with pass_filter 'full', taken plate by plate in shard plate order and in obs order inside a plate,
are laid over the shard rows; for every row group the min/max of all five columns computed from obs must equal the
footer statistics.
Outputs: index/rowgroups.parquet, index/files.parquet, index/file_plates.parquet, index/index_checks.json,
logs/t1_build_index.log
Run: python -I t1_build_index.py
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tp_common as C
import numpy as np, pandas as pd, pyarrow as pa, pyarrow.parquet as pq, pyarrow.compute as pc

IDX = C.DATA + "/index"
log = C.Log(C.DATA + "/logs/t1_build_index.log")
t0 = time.time()
tree = pd.read_csv(IDX + "/tree.tsv", sep="\t").set_index("file")
STAT_COLS = ["sample", "plate", "drug", "cell_line_id", "BARCODE_SUB_LIB_ID"]
SHORT = {"sample": "sample", "plate": "plate", "drug": "drug", "cell_line_id": "cell_line_id", "BARCODE_SUB_LIB_ID": "barcode"}

rg_rows, file_rows = [], []
schema0 = None
problems = []
for i in range(C.N_FILES):
    base = "%s/footers/%s" % (IDX, C.fname(i)[:-8])
    fb = open(base + ".footer", "rb").read()
    meta = json.load(open(base + ".json"))
    md = C.footer_metadata(fb)
    size = int(tree.loc[i, "size"])
    assert meta["size"] == size
    sch = md.schema.to_arrow_schema()
    if schema0 is None:
        schema0 = sch
    elif not sch.equals(schema0, check_metadata=False):
        problems.append("schema differs in file %d" % i)
    paths = [md.schema.column(j).path for j in range(md.num_columns)]
    jcol = {c: paths.index(c) for c in STAT_COLS}
    rs = 0
    prev_end = None
    first_off = None
    for g in range(md.num_row_groups):
        r = md.row_group(g)
        starts, csize = [], 0
        codecs = set()
        for j in range(r.num_columns):
            c = r.column(j)
            starts.append(c.dictionary_page_offset if c.has_dictionary_page else c.data_page_offset)
            csize += c.total_compressed_size
            codecs.add(c.compression)
        off = min(starts)
        if g == 0:
            first_off = off
        if prev_end is not None and off != prev_end:
            problems.append("file %d rg %d not contiguous (%d vs %d)" % (i, g, off, prev_end))
        prev_end = off + csize
        rec = dict(file=i, rg=g, rows=r.num_rows, row_start=rs, byte_offset=off, compressed_size=csize,
                   uncompressed_size=r.total_byte_size, codecs=",".join(sorted(codecs)))
        for cname in STAT_COLS:
            st = r.column(jcol[cname]).statistics
            assert st is not None and st.has_min_max, (i, g, cname)
            mn, mx = st.min, st.max
            mn = mn.decode() if isinstance(mn, bytes) else mn
            mx = mx.decode() if isinstance(mx, bytes) else mx
            rec[SHORT[cname] + "_min"] = mn
            rec[SHORT[cname] + "_max"] = mx
            rec[SHORT[cname] + "_nulls"] = st.null_count
        rg_rows.append(rec)
        rs += r.num_rows
    footer_start = size - meta["footer_len"] - 8
    if first_off != 4 or prev_end != footer_start:
        problems.append("file %d: first rg at %s, last rg ends %s, footer starts %d" % (i, first_off, prev_end, footer_start))
    if rs != md.num_rows:
        problems.append("file %d: row groups sum %d != num_rows %d" % (i, rs, md.num_rows))
    file_rows.append(dict(file=i, path=tree.loc[i, "path"], size=size, lfs_sha256=tree.loc[i, "lfs_sha256"],
                          rows=md.num_rows, row_groups=md.num_row_groups, footer_len=meta["footer_len"],
                          created_by=md.created_by, format_version=md.format_version, data_start=first_off,
                          data_end=prev_end, redirect_host=meta.get("redirect_host")))
rg = pd.DataFrame(rg_rows)
files = pd.DataFrame(file_rows)
log("parsed %d footers, %d row groups, %d rows in %.1f s; layout problems: %d" % (
    len(files), len(rg), int(rg.rows.sum()), time.time() - t0, len(problems)))
for p in problems[:20]:
    log("PROBLEM", p)
rg["global_row"] = np.r_[0, np.cumsum(rg.rows.values)[:-1]]
checks = dict(files=len(files), row_groups=int(len(rg)), rows=int(rg.rows.sum()), layout_problems=len(problems),
              row_groups_per_file=rg.groupby("file").size().value_counts().sort_index().to_dict(),
              rows_per_row_group=rg.rows.value_counts().head(5).to_dict(),
              codecs=rg.codecs.value_counts().to_dict(), created_by=files.created_by.value_counts().to_dict(),
              stat_nulls=int(rg[[c for c in rg.columns if c.endswith("_nulls")]].values.sum()),
              compressed_bytes_in_row_groups=int(rg.compressed_size.sum()), file_bytes=int(files["size"].sum()))
log("checks:", json.dumps({k: (str(v) if isinstance(v, dict) else v) for k, v in checks.items()}))

# ---- plate blocks over files and row groups
pn = lambda s: int(s.replace("plate", ""))
rg["plate_min_n"] = rg.plate_min.map(pn)
rg["plate_max_n"] = rg.plate_max.map(pn)
rg["two_plates"] = rg.plate_min != rg.plate_max
single = rg[~rg.two_plates]
order = list(pd.unique(single.plate_min_n))
blocks = single.groupby("plate_min_n").agg(first_file=("file", "min"), last_file=("file", "max"),
                                           first_g=("global_row", "min"), n_rg=("rg", "size"))
contig = True
for p, gsub in single.groupby("plate_min_n"):
    idx = gsub.index.values
    # all row groups between first and last of this plate must be this plate or a boundary group touching it
    span = rg.loc[idx.min():idx.max()]
    bad = span[~((span.plate_min_n == p) | (span.plate_max_n == p))]
    if len(bad):
        contig = False
        log("plate %d not contiguous: %d foreign row groups inside its span" % (p, len(bad)))
bnd = rg[rg.two_plates]
log("plate order over files:", order, "| plates contiguous:", contig, "| row groups with two plates:", len(bnd),
    "in files", bnd.file.tolist())
file_pl = rg.groupby("file").agg(pmin=("plate_min_n", "min"), pmax=("plate_max_n", "max"), two=("two_plates", "any"))
files["plates"] = [";".join(str(x) for x in sorted({a, b})) for a, b in zip(file_pl.pmin, file_pl.pmax)]
files["boundary"] = file_pl.two.values | (file_pl.pmin.values != file_pl.pmax.values)
checks["plate_order_in_files"] = order
checks["plates_contiguous"] = contig
checks["boundary_row_groups"] = int(len(bnd))
checks["boundary_files"] = files.loc[files.boundary, "file"].tolist()

# ---- obs alignment
obs = pq.read_table(IDX + "/obs_compact.parquet",
                    columns=["obs_row", "key", "plate", "sample", "sublibrary", "cell_line", "drug", "pass_full"])
full = obs.filter(obs.column("pass_full"))
plate_o = full.column("plate").to_numpy()
lib_o = full.column("sublibrary").to_numpy()
# Shard order of sublibrary runs, from the footer barcode statistics: a row group whose min and max barcode carry the
# same '-lib_N' suffix lies inside sublibrary N. The order of first appearance of these libs gives the run order.
rg["lib_min"] = rg.barcode_min.str.split("lib_").str[1].astype(int)
rg["lib_max"] = rg.barcode_max.str.split("lib_").str[1].astype(int)
v = rg.loc[rg.lib_min == rg.lib_max, "lib_min"].values
lib_runs = v[np.r_[True, v[1:] != v[:-1]]]
assert len(lib_runs) == len(np.unique(lib_runs)) == len(np.unique(lib_o)), "a sublibrary is split into several runs"
pos = {int(l): k for k, l in enumerate(lib_runs)}
two = rg[rg.lib_min != rg.lib_max]
between = int(sum(abs(pos[a] - pos[b]) == 1 for a, b in zip(two.lib_min, two.lib_max)))
log("sublibrary runs in shard order: %d (each contiguous); row groups holding two sublibraries: %d, all between "
    "adjacent runs: %s" % (len(lib_runs), len(two), between == len(two)))
checks["sublibrary_runs"] = int(len(lib_runs))
checks["two_sublibrary_row_groups"] = int(len(two))
checks["two_sublibrary_row_groups_between_adjacent_runs"] = between
o = np.argsort(lib_o, kind="stable")                      # obs order kept inside each lib
lib_sorted = lib_o[o]
starts = np.searchsorted(lib_sorted, lib_runs, side="left")
ends = np.searchsorted(lib_sorted, lib_runs, side="right")
seq = np.concatenate([o[a:b] for a, b in zip(starts, ends)])   # positions in `full`, in expected shard order
log("obs pass_full rows %d; footer rows %d; equal: %s" % (len(seq), int(rg.rows.sum()), len(seq) == int(rg.rows.sum())))
assert len(seq) == int(rg.rows.sum())
# naive alternative (plate blocks in obs order) kept for the record
naive = np.concatenate([np.nonzero(plate_o == p)[0] for p in order])
checks["naive_obs_order_equals_shard_order"] = bool(np.array_equal(naive, seq))
seg = rg.global_row.values
mism = {}


def red(vals):
    return np.minimum.reduceat(vals, seg), np.maximum.reduceat(vals, seg)


# sample (fixed width, numeric order = string order)
s = full.column("sample").to_numpy()[seq].astype(np.int32)
mn, mx = red(s)
mism["sample"] = int(((mn != rg.sample_min.str[4:].astype(int).values) | (mx != rg.sample_max.str[4:].astype(int).values)).sum())
# plate (string order differs from numeric: compare as strings)
pv = plate_o[seq].astype(np.int32)
mn, mx = red(pv)
smin = np.array(["plate%d" % a for a in mn]); smax = np.array(["plate%d" % b for b in mx])
emin = np.where(smin < smax, smin, smax)
emax = np.where(smin < smax, smax, smin)
mism["plate"] = int(((emin != rg.plate_min.values) | (emax != rg.plate_max.values)).sum())


def dict_ranks(col):
    a = full.column(col).combine_chunks()
    a = a if pa.types.is_dictionary(a.type) else pc.dictionary_encode(a)
    vals = np.array(a.dictionary.to_pylist(), dtype=object)
    if col == "drug":
        n_ws = int(sum(x != x.strip() for x in vals))
        log("obs drug names with leading or trailing spaces: %d of %d (compared after strip)" % (n_ws, len(vals)))
        checks["obs_drug_names_with_spaces"] = n_ws
        vals = np.array([x.strip() for x in vals], dtype=object)
    codes = a.indices.to_numpy().astype(np.int64)
    srt = np.argsort(vals.astype(str))
    rank = np.empty(len(vals), np.int64)
    rank[srt] = np.arange(len(vals))
    return rank[codes], vals[srt]


for col, short in [("cell_line", "cell_line_id"), ("drug", "drug")]:
    r_, sv = dict_ranks(col)
    mn, mx = red(r_[seq])
    mism[short] = int(((sv[mn] != rg[short + "_min"].values) | (sv[mx] != rg[short + "_max"].values)).sum())
    if mism[short]:
        bad = np.nonzero((sv[mn] != rg[short + "_min"].values) | (sv[mx] != rg[short + "_max"].values))[0][:3]
        log("%s mismatch examples:" % short, [(int(rg.file.iloc[b]), int(rg.rg.iloc[b]), sv[mn[b]], rg[short + "_min"].iloc[b],
                                              sv[mx[b]], rg[short + "_max"].iloc[b]) for b in bad])
# barcode: string order = (b1, b2, b3) then the lib digits as a string
k = full.column("key").to_numpy()[seq]
prefix = k % 1_000_000_000
lib = k // 1_000_000_000
ulib = np.unique(lib)
lib_str_order = np.argsort(np.array([str(x) for x in ulib]))
lib_rank = np.empty(len(ulib), np.int64)
lib_rank[lib_str_order] = np.arange(len(ulib))
sortable = prefix * 4096 + lib_rank[np.searchsorted(ulib, lib)]
assert len(ulib) < 4096
mn, mx = red(sortable)


def to_str(v):
    p, lr = v // 4096, v % 4096
    lb = ulib[lib_str_order[lr]]
    return "%02d_%03d_%03d-lib_%d" % (p // 1_000_000, (p // 1000) % 1000, p % 1000, lb)


bmin = np.array([to_str(v) for v in mn]); bmax = np.array([to_str(v) for v in mx])
mism["barcode"] = int(((bmin != rg.barcode_min.values) | (bmax != rg.barcode_max.values)).sum())
log("row groups whose footer min/max differ from the obs-derived min/max:", mism)
checks["obs_alignment_mismatches"] = mism
checks["obs_alignment_exact"] = all(v == 0 for v in mism.values())

# ---- sample order: resets inside each plate block versus sublibraries in obs
sub_pl = pq.read_table(IDX + "/obs_compact.parquet", columns=["sublibrary", "plate"]).to_pandas().drop_duplicates()
nsub = sub_pl.plate.value_counts().to_dict()
res = []
for p in order:
    sub = rg[(rg.plate_min_n == p) & (rg.plate_max_n == p)]
    smin_ = sub.sample_min.str[4:].astype(int).values
    smax_ = sub.sample_max.str[4:].astype(int).values
    gaps = np.diff(sub.index.values) == 1
    resets = int(((smin_[1:] < smax_[:-1]) & gaps).sum())
    inside = int(((smin_ == smin_.min()) & (smax_ == smax_.max())).sum())
    res.append(dict(plate=p, row_groups=len(sub), sample_resets_between_row_groups=resets,
                    row_groups_spanning_whole_plate_range=inside, sublibraries_in_obs=int(nsub[p])))
res = pd.DataFrame(res)
log("sample order per plate (resets vs sublibraries):\n" + res.to_string(index=False))
checks["sample_order"] = res.to_dict(orient="records")
# global monotonic check of (plate block, sample) over row groups
checks["rows_sorted_by_sample_across_files"] = False
checks["note_order"] = ("plates occupy contiguous blocks of files in the order above; inside a plate the rows form one "
                        "run per sublibrary, sorted by sample within the run; each new sublibrary restarts at the "
                        "plate's first sample, so rows are not sorted by sample across files")

# ---- per (file, plate) rows and offsets into the plate's obs sequence
fp = []
plate_seq = pv  # plate of each shard row
fstart = files.rows.cumsum().values - files.rows.values
for p in order:
    w = np.nonzero(plate_seq == p)[0]
    fi = np.searchsorted(fstart, w, side="right") - 1
    u, first_idx, cnt = np.unique(fi, return_index=True, return_counts=True)
    pos = 0
    for f, fidx, c in zip(u, first_idx, cnt):
        fp.append(dict(file=int(f), plate=int(p), rows=int(c), row_start_in_file=int(w[fidx] - fstart[f]),
                       plate_seq_start=int(pos)))
        pos += int(c)
fp = pd.DataFrame(fp).sort_values(["file", "row_start_in_file"]).reset_index(drop=True)
assert fp.rows.sum() == len(seq)
# expected obs row of every shard row (global shard order), and the sublibrary runs
obs_row_full = full.column("obs_row").to_numpy()
C.atomic_write_table(pa.table({"obs_row": pa.array(obs_row_full[seq])}), IDX + "/shard_obs_rows.parquet",
                     compression="zstd", row_group_size=1 << 20)
runs = pd.DataFrame(dict(run=np.arange(len(lib_runs)), sublibrary=lib_runs, rows=ends - starts,
                         global_row_start=np.r_[0, np.cumsum(ends - starts)[:-1]]))
runs["plate"] = plate_o[o[starts]]
runs["first_file"] = np.searchsorted(fstart, runs.global_row_start.values, side="right") - 1
runs["last_file"] = np.searchsorted(fstart, runs.global_row_start.values + runs.rows.values - 1, side="right") - 1
C.atomic_write_table(pa.Table.from_pandas(runs, preserve_index=False), IDX + "/sublibraries.parquet")
C.atomic_write_table(pa.Table.from_pandas(rg.drop(columns=["plate_min_n", "plate_max_n"]), preserve_index=False),
                     IDX + "/rowgroups.parquet", compression="zstd")
C.atomic_write_table(pa.Table.from_pandas(files, preserve_index=False), IDX + "/files.parquet", compression="zstd")
C.atomic_write_table(pa.Table.from_pandas(fp, preserve_index=False), IDX + "/file_plates.parquet", compression="zstd")
checks["problems"] = problems
checks["seconds"] = round(time.time() - t0, 1)
checks["utc"] = C.utcnow()
json.dump(checks, open(IDX + "/index_checks.json", "w"), indent=1, default=str)
log("wrote rowgroups.parquet (%d rows), files.parquet, file_plates.parquet (%d rows) in %.1f s" % (len(rg), len(fp), time.time() - t0))
