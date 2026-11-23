"""Step 1b: shard index. For every shard: read the footer; for every row group record rows,
compressed bytes and the min/max statistics of drug and sample. Row groups whose statistics admit
DMSO_TF are read for drug, sample, BARCODE_SUB_LIB_ID, cell_line_id and plate (no expression column);
the vehicle rows are kept (row position, barcode, sample, line, plate), the other rows are dropped
in memory at once and never written. Resumable: one output file per shard under data/F2/index_parts/."""
import os, sys, time, json, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import pyarrow.parquet as pq, pyarrow.compute as pc, pyarrow as pa
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hf_io
from common import runlog

OUT = "data/F2/index_parts"
os.makedirs(OUT, exist_ok=True)
fl = pd.read_csv("results/F2/F2_filelist.tsv", sep="\t")
shards = fl[fl.path.str.startswith("data/train-")].sort_values("path").reset_index(drop=True)
sm = pd.read_parquet("data/F2/meta/sample_metadata.parquet")
dmso_samples = sorted(sm.loc[sm.drug == "DMSO_TF", "sample"].tolist())
assert len(dmso_samples) == 28, len(dmso_samples)
SMALL = ["drug", "sample", "BARCODE_SUB_LIB_ID", "cell_line_id", "plate"]


def stat(c):
    st = c.statistics
    if st is None or not st.has_min_max:
        return None, None
    mn, mx = st.min, st.max
    if isinstance(mn, bytes):
        mn, mx = mn.decode(), mx.decode()
    return mn, mx


def do_shard(k):
    path, size = shards.path[k], int(shards["size"][k])
    sid = int(path.split("-")[1])
    fout = f"{OUT}/shard_{sid:05d}.json"
    if os.path.exists(fout):
        return sid, "skip"
    f = hf_io.RangeFile(path, size, tail=1 << 16)
    pf = pq.ParquetFile(f)
    md = pf.metadata
    names = [md.schema.column(j).path for j in range(md.num_columns)]
    idx = {n: j for j, n in enumerate(names)}
    rgs, veh = [], []
    for i in range(md.num_row_groups):
        rg = md.row_group(i)
        dmn, dmx = stat(rg.column(idx["drug"]))
        smn, smx = stat(rg.column(idx["sample"]))
        pmn, pmx = stat(rg.column(idx["plate"]))
        if dmn is None or smn is None:
            cand, why = True, "no_stats"
        else:
            cand = (dmn <= "DMSO_TF" <= dmx) and any(smn <= s <= smx for s in dmso_samples)
            why = "stats_admit" if cand else "stats_exclude"
        rec = dict(shard=sid, row_group=i, rows=rg.num_rows, compressed_bytes=sum(rg.column(j).total_compressed_size for j in range(rg.num_columns)),
                   drug_min=dmn, drug_max=dmx, sample_min=smn, sample_max=smx, plate_min=pmn, plate_max=pmx,
                   candidate=cand, reason=why, vehicle_rows=0)
        if cand:
            cols = [rg.column(idx[c]) for c in SMALL]
            st = min((c.dictionary_page_offset if c.dictionary_page_offset is not None else c.data_page_offset) for c in cols)
            en = max((c.dictionary_page_offset if c.dictionary_page_offset is not None else c.data_page_offset) + c.total_compressed_size for c in cols)
            f.prefetch([(st, en - st)])
            t = pf.read_row_group(i, columns=SMALL)
            m = pc.equal(t.column("drug"), "DMSO_TF")
            pos = pc.indices_nonzero(m)
            tv = t.filter(m).drop_columns(["drug"])
            del t
            rec["vehicle_rows"] = tv.num_rows
            if tv.num_rows:
                dv = tv.to_pandas()
                dv.insert(0, "row_in_group", pos.to_numpy())
                dv.insert(0, "row_group", i)
                dv.insert(0, "shard", sid)
                veh.append(dv)
        rgs.append(rec)
    f.drop()
    if veh:
        pd.concat(veh, ignore_index=True).to_parquet(f"{OUT}/veh_{sid:05d}.parquet", index=False)
    with open(fout + ".tmp", "w") as g:
        json.dump(dict(shard=sid, path=path, size=size, num_rows=md.num_rows, num_row_groups=md.num_row_groups,
                       footer_bytes=md.serialized_size, columns=names, row_groups=rgs, misses=f.misses), g)
    os.replace(fout + ".tmp", fout)
    return sid, sum(r["vehicle_rows"] for r in rgs)


if __name__ == "__main__":
    t0 = time.time()
    todo = list(range(len(shards)))
    done = 0
    with ThreadPoolExecutor(max_workers=2) as ex:
        futs = {ex.submit(do_shard, k): k for k in todo}
        for fu in as_completed(futs):
            sid, r = fu.result()
            done += 1
            if done % 100 == 0:
                print(f"{done}/{len(todo)} shards, {time.time()-t0:.0f}s, http {hf_io.STATS}", flush=True)
    runlog("s02_index", t0, f"shards {len(todo)}; http {hf_io.STATS}")
    print("done", time.time() - t0, hf_io.STATS)
