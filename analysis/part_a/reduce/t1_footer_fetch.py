"""t1_footer_fetch (tahoe-pipeline): read the parquet footer of every Tahoe-100M raw-count file with HTTP range
requests and save it, so the index can be built offline and later passes need no footer reads.

- Repo revision pinned (tp_common.REVISION); the HF API is asked for the current sha and data/ tree at the pinned revision. File sizes and LFS SHA-256 are saved to index/tree.tsv.
- At most 4 requests in flight; exponential backoff with jitter on 429, 5xx and connection errors (tp_common.Http).
- Resumable: a file whose footer is already saved and parses is skipped.
Outputs: index/tree.tsv, index/footers/train-XXXXX.footer (+ .json), logs/t1_footer_fetch.log
Run: python -I t1_footer_fetch.py [--limit N]
"""
import os, sys, json, time, struct, argparse
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tp_common as C
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed

ap = argparse.ArgumentParser()
ap.add_argument("--limit", type=int, default=None)
ap.add_argument("--conc", type=int, default=4)
args = ap.parse_args()
assert args.conc <= 4
IDX = C.DATA + "/index"
FD = IDX + "/footers"
os.makedirs(FD, exist_ok=True)
log = C.Log(C.DATA + "/logs/t1_footer_fetch.log")
http = C.Http(max_conc=args.conc, log=log)

# ---- 1. repo info and tree at the pinned revision
tree_path = IDX + "/tree.tsv"
if not os.path.exists(tree_path):
    r = http.get("https://huggingface.co/api/datasets/%s/revision/%s" % (C.REPO, C.REVISION))
    info = r.json()
    log("dataset revision:", info.get("sha"), "| lastModified:", info.get("lastModified"), "| pinned:", C.REVISION,
        "| same:", info.get("sha") == C.REVISION)
    C.atomic_write_bytes(IDX + "/repo_info.json", r.content)
    url = "https://huggingface.co/api/datasets/%s/tree/%s/data" % (C.REPO, C.REVISION)
    rows, page = [], 0
    while url:
        r = http.get(url)
        page += 1
        for e in r.json():
            if e.get("type") == "file":
                rows.append(dict(path=e["path"], size=e["size"], lfs_sha256=(e.get("lfs") or {}).get("oid"),
                                 lfs_size=(e.get("lfs") or {}).get("size"), xet_hash=e.get("xetHash")))
        url = r.links.get("next", {}).get("url")
        time.sleep(0.5)
    tree = pd.DataFrame(rows)
    tree["file"] = tree.path.str.extract(r"train-(\d{5})-of-03388")[0].astype(int)
    tree = tree.sort_values("file").reset_index(drop=True)
    assert len(tree) == C.N_FILES and tree.file.tolist() == list(range(C.N_FILES))
    assert (tree["size"] == tree["lfs_size"]).all()
    tree.to_csv(tree_path, sep="\t", index=False)
tree = pd.read_csv(tree_path, sep="\t").set_index("file")


# ---- 2. footers
def have(i):
    p = "%s/%s.footer" % (FD, C.fname(i)[:-8])
    if not os.path.exists(p):
        return False
    try:
        meta = C.footer_metadata(open(p, "rb").read())
        with open(p.removesuffix(".footer") + ".json") as f:
            sidecar = json.load(f)
        return sidecar["file"] == i and sidecar["rows"] == meta.num_rows
    except Exception:
        return False


def fetch(i):
    size = int(tree.loc[i, "size"])
    url = C.resolve_url(i)
    t0 = time.time()
    r = http.get(url, headers={"Range": "bytes=-65536"}, stream=True)
    try:
        if r.status_code != 206:
            raise RuntimeError("expected 206, got %d" % r.status_code)
        a, b, total = C.parse_content_range(r.headers["Content-Range"])
        tail = r.content
        host = requests_host(r.url)
    finally:
        C.Http.release(r)
        r.close()
    assert total == size, (i, total, size)
    assert len(tail) == b - a + 1 == min(65536, size), (i, len(tail), a, b)
    assert tail[-4:] == b"PAR1", "bad magic"
    flen = struct.unpack("<I", tail[-8:-4])[0]
    if flen + 8 > len(tail):
        r2 = http.get(url, headers={"Range": "bytes=%d-%d" % (size - flen - 8, size - 1)})
        tail = r2.content
        assert len(tail) == flen + 8
    footer = tail[-(flen + 8):]
    md = C.footer_metadata(footer)
    base = "%s/%s" % (FD, C.fname(i)[:-8])
    C.atomic_write_bytes(base + ".footer", footer)
    C.atomic_write_bytes(base + ".json", json.dumps(dict(
        file=i, size=size, footer_len=flen, rows=md.num_rows, row_groups=md.num_row_groups, redirect_host=host,
        utc=C.utcnow(), seconds=round(time.time() - t0, 2))).encode())
    return i, md.num_rows, md.num_row_groups


def requests_host(u):
    from urllib.parse import urlparse
    return urlparse(u).netloc


todo = [i for i in range(C.N_FILES) if not have(i)]
if args.limit:
    todo = todo[:args.limit]
log("footers already saved: %d; to fetch now: %d" % (C.N_FILES - len([i for i in range(C.N_FILES) if not have(i)]), len(todo)))
t0 = time.time()
done, failed = 0, []
with ThreadPoolExecutor(args.conc) as ex:
    futs = {ex.submit(fetch, i): i for i in todo}
    for f in as_completed(futs):
        i = futs[f]
        try:
            f.result()
            done += 1
        except Exception as e:
            failed.append(i)
            log("FAILED footer %d: %s" % (i, str(e)[:300]))
        if done % 200 == 0 and done:
            log("progress %d/%d in %.0f s; http stats %s" % (done, len(todo), time.time() - t0, http.stats))
log("fetched %d footers in %.0f s; failed %d %s; http stats %s" % (done, time.time() - t0, len(failed), failed[:20], http.stats))
missing = [i for i in range(C.N_FILES) if not have(i)]
log("footers saved now: %d of %d; missing: %s" % (C.N_FILES - len(missing), C.N_FILES, missing[:30]))
sys.exit(0 if not missing else 2)

if failed:
    raise SystemExit(f"{len(failed)} Tahoe footers failed. Re-run to retry them.")
