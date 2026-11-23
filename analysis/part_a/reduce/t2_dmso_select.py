"""t2_dmso_select (tahoe-pipeline): choose the row groups that can hold DMSO (vehicle) cells and estimate the bytes.

1. The 28 DMSO_TF wells from sample_metadata (2 per plate).
2. Row groups whose footer sample range [min, max] contains a DMSO sample (sample ids have fixed width, so string
   order is numeric order).
3. Exact expected DMSO rows per row group from the obs alignment (index/shard_obs_rows.parquet); every row group with
   expected DMSO rows must be selected.
4. Contiguous selected row groups of a file are merged into one byte range.
5. Expected DMSO cells per plate and line from obs_metadata (pass_filter 'full' = the rows in the shards).
Outputs: index/dmso_rowgroups.parquet, index/dmso_ranges.parquet, index/dmso_expected_counts.tsv,
index/dmso_select.json, logs/t2_dmso_select.log
Run: python -I t2_dmso_select.py
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tp_common as C
import numpy as np, pandas as pd, pyarrow as pa, pyarrow.parquet as pq

IDX = C.DATA + "/index"
log = C.Log(C.DATA + "/logs/t2_dmso_select.log")
sm = pd.read_parquet(C.META["sample_metadata"])
dm = sm[sm.drug == "DMSO_TF"][["sample", "plate"]].sort_values("sample")
log("DMSO_TF wells:", len(dm), "| per plate:", dm.plate.value_counts().unique().tolist(), "|", dm["sample"].tolist())
assert len(dm) == 28 and (dm.plate.value_counts() == 2).all()
assert sm["sample"].str.len().eq(8).all() and sm["sample"].str.startswith("smp_").all()
dms = np.array(sorted(int(s[4:]) for s in dm["sample"]))

rg = pd.read_parquet(IDX + "/rowgroups.parquet")
smin = rg.sample_min.str[4:].astype(int).values
smax = rg.sample_max.str[4:].astype(int).values
# a row group can hold a DMSO sample if any DMSO sample lies in [min, max]
k = np.searchsorted(dms, smin, side="left")
can = (k < len(dms)) & (dms[np.minimum(k, len(dms) - 1)] <= smax)
rg["can_hold_dmso"] = can

# exact expected DMSO rows per row group from the alignment
obs = pq.read_table(IDX + "/obs_compact.parquet", columns=["sample", "cell_line", "plate", "pass_full"])
osample = obs.column("sample").to_numpy()
srow = pq.read_table(IDX + "/shard_obs_rows.parquet").column("obs_row").to_numpy()
is_dmso_shard = np.isin(osample[srow], dms)
seg = rg.global_row.values
rg["expected_dmso_rows"] = np.add.reduceat(is_dmso_shard.astype(np.int64), seg)
missed = int(((rg.expected_dmso_rows > 0) & ~rg.can_hold_dmso).sum())
sel = rg[rg.can_hold_dmso].copy()
log("row groups selected by sample range: %d of %d (%.2f%%); with expected DMSO rows: %d; selected but expected 0: %d; "
    "expected DMSO but not selected: %d" % (len(sel), len(rg), 100 * len(sel) / len(rg),
                                             int((rg.expected_dmso_rows > 0).sum()),
                                             int((sel.expected_dmso_rows == 0).sum()), missed))
assert missed == 0
tot_dmso = int(rg.expected_dmso_rows.sum())
log("DMSO rows expected in shards: %d; rows in selected row groups: %d (DMSO share %.3f)" % (
    tot_dmso, int(sel.rows.sum()), tot_dmso / sel.rows.sum()))

# merge contiguous selected row groups of a file into byte ranges
ranges = []
for f, g in sel.groupby("file"):
    g = g.sort_values("rg")
    cur = None
    for r in g.itertuples():
        if cur is not None and r.rg == cur["rg_last"] + 1:
            cur["rg_last"] = r.rg
            cur["nbytes"] = r.byte_offset + r.compressed_size - cur["start"]
            cur["rows"] += r.rows
            cur["expected_dmso_rows"] += r.expected_dmso_rows
        else:
            if cur is not None:
                ranges.append(cur)
            cur = dict(file=f, rg_first=r.rg, rg_last=r.rg, start=r.byte_offset, nbytes=r.compressed_size,
                       rows=r.rows, expected_dmso_rows=r.expected_dmso_rows)
    ranges.append(cur)
ranges = pd.DataFrame(ranges)
assert ranges.nbytes.sum() == sel.compressed_size.sum()
files = pd.read_parquet(IDX + "/files.parquet").set_index("file")
summ = dict(dmso_wells=dm["sample"].tolist(), row_groups_selected=int(len(sel)), row_groups_total=int(len(rg)),
            files_touched=int(sel.file.nunique()), byte_ranges=int(len(ranges)),
            bytes_selected=int(sel.compressed_size.sum()), bytes_all_files=int(files["size"].sum()),
            share_of_bytes=round(float(sel.compressed_size.sum() / files["size"].sum()), 4),
            rows_selected=int(sel.rows.sum()), dmso_rows_expected=tot_dmso,
            selected_row_groups_without_dmso=int((sel.expected_dmso_rows == 0).sum()),
            hours_at_9_04_MBps=round(sel.compressed_size.sum() / 9.04e6 / 3600, 2),
            hours_at_7_11_MBps=round(sel.compressed_size.sum() / 7.11e6 / 3600, 2))
log(json.dumps(summ)[:2000])
# per plate
pp = sel.assign(plate=np.where(sel.plate_min == sel.plate_max, sel.plate_min, sel.plate_min + ";" + sel.plate_max)).groupby(
    "plate").agg(row_groups=("rg", "size"), bytes=("compressed_size", "sum"), rows=("rows", "sum"),
                 dmso_rows=("expected_dmso_rows", "sum"))
log("per plate:\n" + pp.to_string())

# expected DMSO cells per plate x line, from obs (pass_filter 'full' and all cells)
cl = obs.column("cell_line").to_pandas().astype(str).values
pl = obs.column("plate").to_numpy()
pf = obs.column("pass_full").to_numpy()
m = np.isin(osample, dms)
df = pd.DataFrame(dict(plate=pl[m], sample=osample[m], cell_line_id=cl[m], full=pf[m]))
exp = df.groupby(["plate", "cell_line_id"]).agg(cells_all=("full", "size"), cells_full=("full", "sum")).reset_index()
clm = pd.read_parquet(C.META["cell_line_metadata"])[["cell_name", "Cell_ID_Cellosaur"]].drop_duplicates()
exp = exp.merge(clm, left_on="cell_line_id", right_on="Cell_ID_Cellosaur", how="left").drop(columns="Cell_ID_Cellosaur")
exp.to_csv(IDX + "/dmso_expected_counts.tsv", sep="\t", index=False)
per_line = exp.groupby("cell_name")[["cells_all", "cells_full"]].sum()
log("DMSO cells, pass_filter 'full': total %d; A549 %d; SW48 %d | all pass_filter values: A549 %d; SW48 %d" % (
    int(exp.cells_full.sum()), int(per_line.loc["A549", "cells_full"]), int(per_line.loc["SW48", "cells_full"]),
    int(per_line.loc["A549", "cells_all"]), int(per_line.loc["SW48", "cells_all"])))
log("DMSO cells per line ('full'): min %d (%s), median %d, max %d" % (
    per_line.cells_full.min(), per_line.cells_full.idxmin(), per_line.cells_full.median(), per_line.cells_full.max()))
summ.update(dmso_cells_full=int(exp.cells_full.sum()), dmso_cells_all_pass_filter=int(exp.cells_all.sum()),
            A549_full=int(per_line.loc["A549", "cells_full"]), SW48_full=int(per_line.loc["SW48", "cells_full"]),
            A549_all=int(per_line.loc["A549", "cells_all"]), SW48_all=int(per_line.loc["SW48", "cells_all"]),
            lines=int(per_line.shape[0]), utc=C.utcnow())
assert summ["dmso_cells_full"] == tot_dmso
C.atomic_write_table(pa.Table.from_pandas(sel, preserve_index=False), IDX + "/dmso_rowgroups.parquet")
C.atomic_write_table(pa.Table.from_pandas(ranges, preserve_index=False), IDX + "/dmso_ranges.parquet")
json.dump(summ, open(IDX + "/dmso_select.json", "w"), indent=1)
