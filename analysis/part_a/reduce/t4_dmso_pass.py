"""t4_dmso_pass (tahoe-pipeline): selective pass over all DMSO (vehicle) cells of Tahoe-100M.

For every file with DMSO row groups (index/dmso_ranges.parquet), fetch only those byte ranges (HTTP Range, at most
4 requests in flight, backoff on 429/5xx), parse them with pyarrow using the saved footer (tp_common.SparseFile raises
on any read outside the fetched bytes), and run the shared reducer with keep_samples = the 28 DMSO wells.
Per-file partials go to dmso/work/ with a manifest (resumable). When all files are done, finalize():
  cells -> obs join and checks -> dmso/cells_markers.parquet
  partial pseudobulks -> one sum per (cell line, sample, stratum) -> dmso/pseudobulk_counts.parquet, groups, phases
  six-line rows -> dmso/six_lines/plateN.parquet (original columns, SNAPPY)
  QC: DMSO cells per line per plate against obs_metadata; bytes, time, requests.
No sex state is called and no drug-response statistic is computed.
Run: python -I t4_dmso_pass.py [--workers 4] [--limit N] [--finalize-only]
"""
import os, sys, json, time, glob, argparse, threading
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tp_common as C
import tp_reducer as R
import numpy as np, pandas as pd, pyarrow as pa, pyarrow.parquet as pq, pyarrow.compute as pc
from concurrent.futures import ThreadPoolExecutor, as_completed

ap = argparse.ArgumentParser()
ap.add_argument("--workers", type=int, default=4)
ap.add_argument("--limit", type=int, default=None)
ap.add_argument("--finalize-only", action="store_true")
args = ap.parse_args()
assert args.workers <= 4
IDX = C.DATA + "/index"
DM = C.DATA + "/dmso"
WORK = DM + "/work"
os.makedirs(WORK, exist_ok=True)
log = C.Log(C.DATA + "/logs/t4_dmso_pass.log")
MAN = WORK + "/manifest.jsonl"
mlock = threading.Lock()

ranges = pd.read_parquet(IDX + "/dmso_ranges.parquet")
selrg = pd.read_parquet(IDX + "/dmso_rowgroups.parquet")
files = pd.read_parquet(IDX + "/files.parquet").set_index("file")
sm = pd.read_parquet(C.META["sample_metadata"])
DMSO = sorted(int(s[4:]) for s in sm.loc[sm.drug == "DMSO_TF", "sample"])
assert len(DMSO) == 28


def load_manifest():
    done = {}
    if os.path.exists(MAN):
        for line in open(MAN):
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue
            done[r["file"]] = r
    return done


def outputs_exist(fi):
    b = "%s/%s" % (WORK, C.fname(fi)[:-8])
    return all(os.path.exists(b + s) for s in [".cells.parquet", ".pb.parquet", ".groups.parquet", ".six.parquet"])


red = R.Reducer()
http = C.Http(max_conc=args.workers, log=log)
stat = dict(bytes=0, files=0, t0=time.time())
slock = threading.Lock()


def one(fi):
    t0 = time.time()
    size = int(files.loc[fi, "size"])
    md = C.footer_metadata(open("%s/footers/%s.footer" % (IDX, C.fname(fi)[:-8]), "rb").read())
    rr = ranges[ranges.file == fi].sort_values("start")
    bufs = []
    for r in rr.itertuples():
        resp = http.get(C.resolve_url(fi), headers={"Range": "bytes=%d-%d" % (r.start, r.start + r.nbytes - 1)})
        if resp.status_code != 206:
            raise RuntimeError("expected 206, got %d" % resp.status_code)
        a, b, total = C.parse_content_range(resp.headers["Content-Range"])
        body = resp.content
        if not (a == r.start and b == r.start + r.nbytes - 1 and total == size and len(body) == r.nbytes):
            raise RuntimeError("range mismatch: %s, got %d bytes" % (resp.headers["Content-Range"], len(body)))
        bufs.append((int(r.start), body))
    nbytes = sum(len(b) for _, b in bufs)
    sf = C.SparseFile(size, bufs)
    rgs = sorted(selrg.loc[selrg.file == fi, "rg"].tolist())
    tb = pq.ParquetFile(sf, metadata=md).read_row_groups(rgs)
    rg_of, row_of, nexp = R.rows_meta_for(md, rgs)
    if tb.num_rows != nexp:
        raise RuntimeError("rows read %d != footer rows %d" % (tb.num_rows, nexp))
    out = red.reduce(tb, fi, rg_of, row_of, keep_samples=DMSO, want_six=True)
    exp_dmso = int(selrg.loc[selrg.file == fi, "expected_dmso_rows"].sum())
    if out["stats"]["rows_kept"] != exp_dmso:
        raise RuntimeError("DMSO rows %d != expected %d" % (out["stats"]["rows_kept"], exp_dmso))
    base = "%s/%s" % (WORK, C.fname(fi)[:-8])
    R.write_partials(out, base, six_path=base + ".six.parquet")
    rec = dict(file=int(fi), status="done", bytes=int(nbytes), requests=len(bufs), row_groups=len(rgs),
               rows_read=int(tb.num_rows), rows_footer=int(nexp), dmso_rows=out["stats"]["rows_kept"],
               dmso_expected=exp_dmso, sparse_reads=len(sf.reads), seconds=round(time.time() - t0, 2),
               utc=C.utcnow(), stats=out["stats"])
    with mlock:
        with open(MAN, "a") as f:
            f.write(json.dumps(rec) + "\n")
    with slock:
        stat["bytes"] += nbytes
        stat["files"] += 1
    return rec


def run_pass():
    done = load_manifest()
    todo = [int(f) for f in sorted(ranges.file.unique()) if not (done.get(int(f), {}).get("status") == "done" and outputs_exist(int(f)))]
    if args.limit:
        todo = todo[:args.limit]
    need = int(ranges[ranges.file.isin(todo)].nbytes.sum())
    log("DMSO pass: files with DMSO row groups %d; already done %d; to do now %d (%.2f GB)" % (
        ranges.file.nunique(), ranges.file.nunique() - len(todo) if not args.limit else len(done), len(todo), need / 1e9))
    t0 = time.time()
    failed = []
    with ThreadPoolExecutor(args.workers) as ex:
        futs = {ex.submit(one, f): f for f in todo}
        k = 0
        for fu in as_completed(futs):
            f = futs[fu]
            k += 1
            try:
                fu.result()
            except Exception as e:
                failed.append(f)
                log("FAILED file %d: %s" % (f, str(e)[:300]))
                with mlock:
                    with open(MAN, "a") as fh:
                        fh.write(json.dumps(dict(file=int(f), status="failed", error=str(e)[:300], utc=C.utcnow())) + "\n")
            if k % 50 == 0 or k == len(todo):
                el = time.time() - t0
                rate = stat["bytes"] / max(el, 1e-9)
                left = need - stat["bytes"]
                log("progress %d/%d files, %.2f GB, %.2f MB/s, ETA %.1f min, http %s" % (
                    k, len(todo), stat["bytes"] / 1e9, rate / 1e6, left / max(rate, 1) / 60, http.stats))
    el = time.time() - t0
    log("pass done: %d files, %d failed, %.2f GB in %.0f s (%.2f MB/s); http %s" % (
        len(todo) - len(failed), len(failed), stat["bytes"] / 1e9, el, stat["bytes"] / max(el, 1e-9) / 1e6, http.stats))
    return failed, el


def finalize():
    t0 = time.time()
    done = load_manifest()
    allf = sorted(int(f) for f in ranges.file.unique())
    missing = [f for f in allf if not (done.get(f, {}).get("status") == "done" and outputs_exist(f))]
    if missing:
        log("finalize skipped: %d files not done, e.g. %s" % (len(missing), missing[:10]))
        return False
    recs = pd.DataFrame([done[f] for f in allf])
    lines = red.lines
    # ---- cells
    cells = pd.concat([pq.read_table("%s/%s.cells.parquet" % (WORK, C.fname(f)[:-8])).to_pandas() for f in allf],
                      ignore_index=True)
    ob = R.ObsIndex()
    cells, jchk = R.join_cells(cells, ob)
    del ob
    cells["cell_name"] = cells.cell_line_id.astype(str).map(dict(zip(lines.cell_line_id, lines.cell_name)))
    order = ["plate", "sample", "drug", "cell_line_id", "cell_name", "BARCODE_SUB_LIB_ID", "key", "total_umi", "n_genes"] + \
            C.MARKERS + ["y8", "stratum", "phase", "obs_gene_count", "obs_tscp_count", "pcnt_mito", "file", "rg", "row", "obs_row"]
    cells = cells[order]
    cells["_p"] = C.plate_num(pa.array(cells.plate.astype(str).values))
    cells = cells.sort_values(["_p", "sample", "cell_line_id", "file", "row"]).drop(columns="_p").reset_index(drop=True)
    # ---- pseudobulk
    grp = pd.concat([pq.read_table("%s/%s.groups.parquet" % (WORK, C.fname(f)[:-8])).to_pandas() for f in allf],
                    ignore_index=True)
    dmg = R.DenseMerge(grp.gkey.values, dtype=np.int64)
    for f in allf:
        dmg.add(pq.read_table("%s/%s.pb.parquet" % (WORK, C.fname(f)[:-8])))
    merged = dmg.result()
    del dmg
    fin = R.finalize_tables(merged, grp, cells, lines)
    pbt = fin["pseudobulk"]
    plate_of_sample = dict(zip(sm["sample"], sm["plate"]))
    samp = pbt.column("sample").combine_chunks()
    pl = pa.array([plate_of_sample[s] for s in samp.dictionary.to_pylist()])
    pbt = pbt.add_column(0, "plate", pa.DictionaryArray.from_arrays(samp.indices, pl))
    # ---- QC counts against obs
    exp = pd.read_csv(IDX + "/dmso_expected_counts.tsv", sep="\t")
    got = cells.groupby(["plate", "cell_line_id"], observed=True).agg(
        dmso_cells=("key", "size"), median_total_umi=("total_umi", "median"), median_n_genes=("n_genes", "median")).reset_index()
    exp["plate"] = "plate" + exp.plate.astype(str)
    qc = exp.merge(got, on=["plate", "cell_line_id"], how="outer")
    qc["dmso_cells"] = qc.dmso_cells.fillna(0).astype(int)
    qc["equal_obs_full"] = qc.dmso_cells == qc.cells_full
    per_line = qc.groupby("cell_name")[["dmso_cells", "cells_full", "cells_all"]].sum().sort_values("dmso_cells")
    checks = dict(jchk)
    checks.update(fin["checks"])
    checks.update(files=len(allf), rows_read=int(recs.rows_read.sum()), rows_footer=int(recs.rows_footer.sum()),
                  rows_read_equal_footer=bool((recs.rows_read == recs.rows_footer).all()),
                  dmso_rows=int(recs.dmso_rows.sum()), dmso_rows_expected=int(recs.dmso_expected.sum()),
                  line_plate_pairs=int(len(qc)), line_plate_pairs_equal_obs=int(qc.equal_obs_full.sum()),
                  A549_dmso_cells=int(per_line.loc["A549", "dmso_cells"]), SW48_dmso_cells=int(per_line.loc["SW48", "dmso_cells"]))
    # ---- six lines per plate
    os.makedirs(DM + "/six_lines", exist_ok=True)
    six_rows = 0
    six_tabs = {}
    for f in allf:
        t = pq.read_table("%s/%s.six.parquet" % (WORK, C.fname(f)[:-8]))
        if t.num_rows == 0:
            continue
        for p in pc.unique(t.column("plate")).to_pylist():
            six_tabs.setdefault(p, []).append(t.filter(pc.equal(t.column("plate"), p)))
    for p, tabs in six_tabs.items():
        t = pa.concat_tables(tabs)
        six_rows += t.num_rows
        C.atomic_write_table(t, "%s/six_lines/%s.parquet" % (DM, p), compression="snappy", row_group_size=10000)
    six_ids = set(lines.loc[lines.six_line, "cell_line_id"])
    checks["six_line_rows"] = int(six_rows)
    checks["six_line_rows_expected"] = int(cells.cell_line_id.astype(str).isin(six_ids).sum())
    # ---- write
    C.atomic_write_table(pa.Table.from_pandas(cells, preserve_index=False), DM + "/cells_markers.parquet", compression="zstd")
    C.atomic_write_table(pbt, DM + "/pseudobulk_counts.parquet", compression="zstd", row_group_size=1 << 20)
    fin["groups"].to_csv(DM + "/pseudobulk_groups.tsv", sep="\t", index=False)
    fin["phase"].to_csv(DM + "/phase_counts.tsv", sep="\t", index=False)
    qc.to_csv(DM + "/qc_dmso_cells_per_line_plate.tsv", sep="\t", index=False)
    gm = pd.read_parquet(C.META["gene_metadata"])
    gm.to_csv(DM + "/genes.tsv", sep="\t", index=False)
    lines.to_csv(DM + "/lines.tsv", sep="\t", index=False)
    summ = dict(checks=checks, bytes_fetched=int(recs.bytes.sum()), http_requests_data=int(recs.requests.sum()),
                row_groups_read=int(recs.row_groups.sum()), files=len(allf), pass_seconds_sum_per_file=float(recs.seconds.sum()),
                finalize_seconds=round(time.time() - t0, 1), utc=C.utcnow())
    json.dump(summ, open(DM + "/dmso_pass_summary.json", "w"), indent=1, default=str)
    log("finalize checks:", json.dumps(checks, default=str))
    log("DMSO cells per line (read / obs full):\n" + per_line.to_string())
    return True


if __name__ == "__main__":
    if not args.finalize_only:
        failed, el = run_pass()
        json.dump(dict(failed=failed, seconds=round(el, 1), http=http.stats, utc=C.utcnow()),
                  open(WORK + "/last_run.json", "w"))
        if failed:
            raise SystemExit(f"{len(failed)} Tahoe expression shards failed. Re-run this stage to retry them.")
    if args.limit is None and not finalize():
        raise SystemExit("Could not finalize the complete Tahoe reduction")
