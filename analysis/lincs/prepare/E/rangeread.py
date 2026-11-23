"""HTTP range reads of whole rows from a contiguous, uncompressed HDF5 matrix (LINCS gctx).

The matrix layout (offset, shape, dtype) is read from the file with h5py over fsspec, then rows are
fetched with plain HTTP Range requests, merging nearby rows into one request.
"""
import time, threading
from concurrent.futures import ThreadPoolExecutor
import numpy as np, requests, fsspec, h5py

_local = threading.local()


def session():
    if not hasattr(_local, "s"):
        _local.s = requests.Session()
    return _local.s


def layout(url, path="0/DATA/0/matrix"):
    fs = fsspec.filesystem("https")
    with fs.open(url, "rb", block_size=2**20, cache_type="bytes") as f:
        h = h5py.File(f, "r")
        d = h[path]
        assert d.chunks is None and d.compression is None, "matrix is not contiguous"
        return dict(offset=int(d.id.get_offset()), shape=tuple(d.shape), dtype=d.dtype.newbyteorder("<") if d.dtype.byteorder == "=" else d.dtype)


def h5py_rows(url, rows, path="0/DATA/0/matrix"):
    """Reference read through h5py (used only for the equality check)."""
    fs = fsspec.filesystem("https")
    with fs.open(url, "rb", block_size=2**20, cache_type="bytes") as f:
        h = h5py.File(f, "r")
        return np.stack([h[path][int(r), :] for r in rows])


def fetch(url, start, end, tries=6):
    for k in range(tries):
        try:
            r = session().get(url, headers={"Range": f"bytes={start}-{end}"}, timeout=120)
            if r.status_code == 206 and len(r.content) == end - start + 1:
                return r.content
            err = f"status {r.status_code} len {len(r.content)}"
        except Exception as e:  # network error: retry with backoff
            err = repr(e)
        time.sleep(2 ** k)
    raise RuntimeError(f"range read failed {start}-{end}: {err}")


def plan_runs(rows, max_gap=3, max_len=64):
    rows = np.unique(np.asarray(rows, dtype=np.int64))
    runs = []
    s = p = rows[0]
    for r in rows[1:]:
        if r - p <= max_gap + 1 and r - s < max_len:
            p = r
        else:
            runs.append((s, p)); s = p = r
    runs.append((s, p))
    return runs


def read_rows(url, lay, rows, cols=None, threads=24, log=None, label=""):
    """Return array (len(rows), len(cols)) in the order of `rows`."""
    rows = np.asarray(rows, dtype=np.int64)
    n_genes = lay["shape"][1]
    dt = np.dtype(lay["dtype"])
    rb = n_genes * dt.itemsize
    runs = plan_runs(rows)
    out_cols = n_genes if cols is None else len(cols)
    pos = {int(r): i for i, r in enumerate(rows)}
    assert len(pos) == len(rows), "duplicate rows requested"
    out = np.empty((len(rows), out_cols), dtype=dt)
    t0 = time.time(); done = [0, 0]
    lock = threading.Lock()

    def job(run):
        s, e = run
        buf = fetch(url, lay["offset"] + s * rb, lay["offset"] + (e + 1) * rb - 1)
        block = np.frombuffer(buf, dtype=dt).reshape(e - s + 1, n_genes)
        for r in range(s, e + 1):
            i = pos.get(r)
            if i is not None:
                out[i] = block[r - s] if cols is None else block[r - s, cols]
        with lock:
            done[0] += 1; done[1] += len(buf)
            if log and done[0] % 2000 == 0:
                el = time.time() - t0
                log(f"{label} {done[0]}/{len(runs)} requests, {done[1]/1e9:.2f} GB, {done[1]/1e6/el:.1f} MB/s")

    with ThreadPoolExecutor(threads) as ex:
        list(ex.map(job, runs))
    if log:
        el = time.time() - t0
        log(f"{label} done: {len(rows)} rows, {len(runs)} requests, {done[1]/1e9:.3f} GB in {el:.0f}s")
    return out, done[1]
