"""tp_common: shared constants and helpers for the Tahoe-100M blind pipeline (label tahoe-pipeline).

Nothing here calls a sex state or computes a drug-response statistic.
Every script imports this module after putting its own directory on sys.path (needed under python -I).
"""
import os, sys, time, json, random, struct, threading, hashlib, shutil
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.compute as pc
import requests

# ---------------------------------------------------------------- locations
DATA = os.path.join(os.environ["P59_DATA"], "tahoe")
META = {k: DATA + "/index/meta/%s.parquet" % k for k in
        ["sample_metadata", "cell_line_metadata", "gene_metadata", "drug_metadata"]}
META["obs_metadata"] = DATA + "/index/meta/obs_metadata.parquet"


# ---------------------------------------------------------------- dataset
REPO = "tahoebio/Tahoe-100M"
REVISION = "349e2a55f5538a3e73b32a06c6bdd8e9e87011e8"   # repo sha seen by the pilot; re-checked by t1 against the API
N_FILES = 3388
UA = "generic-research-client/1.0"


def fname(i):
    return "train-%05d-of-03388.parquet" % int(i)


def resolve_url(i, rev=REVISION):
    return "https://huggingface.co/datasets/%s/resolve/%s/data/%s" % (REPO, rev, fname(i))


# ---------------------------------------------------------------- genes
Y8 = ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "ZFY", "USP9Y", "NLGN4Y"]
YX = ["TXLNGY", "PRKY"]
XG = ["RPS4X", "DDX3X", "KDM5C", "KDM6A", "EIF1AX", "ZFX", "USP9X", "NLGN4X"]
HK = ["ACTB", "GAPDH"]
MARKERS = ["XIST", "TSIX"] + Y8 + YX + XG + HK          # TSIX is extra (the pilot kept it too)
N_SPECIAL = 3          # token ids 0, 1, 2 are not in gene_metadata (special tokens); genes are 3..62712
CLS_TOKEN = 1          # first entry of every cell, value -2.0
N_TOKENS = 62713       # max token id + 1
Y_MIN, XIST_MIN = 3, 3  # stratum thresholds (proposal 3a, misassigned-cell filter)

# six lines kept as full per-cell rows: pilot mosaic lines and homogeneous reference lines
SIX_NAMES = ["A549", "SW48", "NCI-H460", "LoVo", "RKO", "HT-29"]

PHASES = ["G1", "S", "G2M"]
PASS_FILTERS = ["full", "minimal"]


def stratum_code(y8, xist):
    """bit 1 (value 2): Y UMIs over the 8 Y genes >= 3; bit 0 (value 1): XIST UMIs >= 3."""
    return (2 * (np.asarray(y8) >= Y_MIN) + (np.asarray(xist) >= XIST_MIN)).astype(np.int8)


STRATUM_LABEL = {0: "y_lt3_xist_lt3", 1: "y_lt3_xist_ge3", 2: "y_ge3_xist_lt3", 3: "y_ge3_xist_ge3"}


def marker_lut(gene_metadata_path=META["gene_metadata"]):
    gm = pd.read_parquet(gene_metadata_path)
    assert gm.token_id.is_unique and gm.gene_symbol.is_unique
    assert int(gm.token_id.min()) == N_SPECIAL and int(gm.token_id.max()) == N_TOKENS - 1
    lut = np.full(N_TOKENS, -1, dtype=np.int32)
    for k, g in enumerate(MARKERS):
        tok = gm.loc[gm.gene_symbol == g, "token_id"].astype(int).tolist()
        assert len(tok) == 1, (g, tok)
        lut[tok[0]] = k
    return lut


def six_line_ids(clm_path=META["cell_line_metadata"]):
    clm = pd.read_parquet(clm_path)[["cell_name", "Cell_ID_Cellosaur"]].drop_duplicates()
    out = {}
    for nm in SIX_NAMES:
        ids = clm.loc[clm.cell_name == nm, "Cell_ID_Cellosaur"].unique().tolist()
        assert len(ids) == 1, (nm, ids)
        out[nm] = ids[0]
    return out


# ---------------------------------------------------------------- keys
def encode_keys(arr):
    """BARCODE_SUB_LIB_ID 'BB_BBB_BBB-lib_N' -> int64 lib*1e9 + b1*1e6 + b2*1e3 + b3 (lossless for this format)."""
    arr = arr.combine_chunks() if isinstance(arr, pa.ChunkedArray) else arr
    ok = pc.match_substring_regex(arr, r"^\d{2}_\d{3}_\d{3}-lib_\d+$")
    assert pc.all(ok).as_py(), "unexpected BARCODE_SUB_LIB_ID format"
    i64 = pa.int64()
    b1 = pc.cast(pc.utf8_slice_codeunits(arr, 0, 2), i64).to_numpy()
    b2 = pc.cast(pc.utf8_slice_codeunits(arr, 3, 6), i64).to_numpy()
    b3 = pc.cast(pc.utf8_slice_codeunits(arr, 7, 10), i64).to_numpy()
    lib = pc.cast(pc.utf8_slice_codeunits(arr, 15), i64).to_numpy()
    return lib * 1_000_000_000 + b1 * 1_000_000 + b2 * 1_000 + b3


def sample_num(arr):
    arr = arr.combine_chunks() if isinstance(arr, pa.ChunkedArray) else arr
    return pc.cast(pc.utf8_slice_codeunits(arr, 4), pa.int32()).to_numpy()


def plate_num(arr):
    arr = arr.combine_chunks() if isinstance(arr, pa.ChunkedArray) else arr
    return pc.cast(pc.utf8_slice_codeunits(arr, 5), pa.int32()).to_numpy()


# ---------------------------------------------------------------- logging and files
def utcnow():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


class Log:
    def __init__(self, path, also_print=True):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        self.f = open(path, "a")
        self.p = also_print
        self.lock = threading.Lock()

    def __call__(self, *a):
        s = "%s %s" % (utcnow(), " ".join(str(x) for x in a))
        with self.lock:
            if self.p:
                print(s, flush=True)
            self.f.write(s + "\n")
            self.f.flush()


def atomic_write_bytes(path, data):
    tmp = path + ".tmp%d" % os.getpid()
    with open(tmp, "wb") as f:
        f.write(data)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def atomic_write_table(table, path, **kw):
    tmp = path + ".tmp%d" % os.getpid()
    pq.write_table(table, tmp, **kw)
    os.replace(tmp, path)


def sha256_file(path, bufsize=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(bufsize), b""):
            h.update(chunk)
    return h.hexdigest()


def free_gb(path):
    return shutil.disk_usage(path).free / 1e9


# ---------------------------------------------------------------- HTTP with backoff
class Http:
    """GET with at most `max_conc` requests in flight, exponential backoff with jitter on 429, 5xx and
    connection errors, Retry-After honoured, and a shared pause for all threads after a 429."""

    def __init__(self, max_conc=4, log=None, base=2.0, cap=300.0, max_attempts=10):
        self.sem = threading.BoundedSemaphore(max_conc)
        self.local = threading.local()
        self.lock = threading.Lock()
        self.pause_until = 0.0
        self.base, self.cap, self.max_attempts = base, cap, max_attempts
        self.log = log or (lambda *a: None)
        self.stats = dict(requests=0, http_429=0, http_5xx=0, conn_errors=0, other_status=0)

    def session(self):
        s = getattr(self.local, "s", None)
        if s is None:
            s = requests.Session()
            s.headers["User-Agent"] = UA
            self.local.s = s
        return s

    def _bump(self, k):
        with self.lock:
            self.stats[k] += 1

    def _delay(self, attempt, retry_after=None):
        d = min(self.cap, self.base * (2 ** attempt)) * (0.5 + random.random())
        if retry_after:
            try:
                d = max(d, float(retry_after))
            except ValueError:
                pass
        return d

    def get(self, url, headers=None, stream=False, timeout=(30, 120), ok=(200, 206)):
        last = None
        for attempt in range(self.max_attempts):
            w = self.pause_until - time.time()
            if w > 0:
                time.sleep(w)
            self.sem.acquire()
            try:
                self._bump("requests")
                ranged_url = url + "?range=" + headers["Range"].removeprefix("bytes=") if headers and "Range" in headers else url
                r = self.session().get(ranged_url, headers=headers, stream=stream, timeout=timeout, allow_redirects=True)
            except (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError) as e:
                self.sem.release()
                self._bump("conn_errors")
                last = "conn: %s" % str(e)[:200]
                d = self._delay(attempt)
                self.log("GET retry %d after %s; sleep %.1f s" % (attempt + 1, last, d))
                time.sleep(d)
                continue
            if stream and r.status_code in ok:
                # caller reads the body and must call release() when done
                r._tp_sem = self.sem
                return r
            self.sem.release()
            if r.status_code in ok:
                return r
            code = r.status_code
            ra = r.headers.get("Retry-After")
            r.close()
            if code == 429 or code >= 500:
                self._bump("http_429" if code == 429 else "http_5xx")
                d = self._delay(attempt, ra)
                if code == 429:
                    with self.lock:
                        self.pause_until = max(self.pause_until, time.time() + d)
                last = "HTTP %d" % code
                self.log("GET retry %d after %s (Retry-After %s); sleep %.1f s" % (attempt + 1, last, ra, d))
                time.sleep(d)
                continue
            self._bump("other_status")
            raise RuntimeError("HTTP %d for %s" % (code, url))
        raise RuntimeError("gave up after %d attempts: %s" % (self.max_attempts, last))

    @staticmethod
    def release(r):
        s = getattr(r, "_tp_sem", None)
        if s is not None:
            r._tp_sem = None
            s.release()


def parse_content_range(h):
    # 'bytes 0-99/1234'
    unit, rest = h.split(" ", 1)
    rng, total = rest.split("/")
    a, b = rng.split("-")
    return int(a), int(b), int(total)


def footer_metadata(footer_bytes):
    """footer_bytes = thrift FileMetaData + 4-byte length + b'PAR1' (the last flen+8 bytes of the file)."""
    return pq.read_metadata(pa.BufferReader(b"PAR1" + footer_bytes))


# ---------------------------------------------------------------- sparse remote file for pyarrow
class SparseFile:
    """Read-only file object that holds only some byte ranges of a remote file.
    Any read outside the held ranges raises, which proves that pyarrow touched only those bytes."""

    def __init__(self, size, ranges):
        self.size = size
        self.ranges = sorted(ranges)  # list of (start, bytes)
        self.pos = 0
        self.closed = False
        self.reads = []

    def seekable(self):
        return True

    def readable(self):
        return True

    def writable(self):
        return False

    def tell(self):
        return self.pos

    def seek(self, off, whence=0):
        if whence == 0:
            self.pos = off
        elif whence == 1:
            self.pos += off
        else:
            self.pos = self.size + off
        return self.pos

    def read(self, n=-1):
        if n is None or n < 0:
            n = self.size - self.pos
        a, b = self.pos, min(self.pos + n, self.size)
        for s, buf in self.ranges:
            if s <= a and b <= s + len(buf):
                self.pos = b
                self.reads.append((a, b - a))
                return buf[a - s:b - s]
        # a read may span several fetched ranges that touch each other
        parts, cur = [], a
        for s, buf in self.ranges:
            if s <= cur < s + len(buf):
                take = min(b, s + len(buf)) - cur
                parts.append(buf[cur - s:cur - s + take])
                cur += take
                if cur == b:
                    self.pos = b
                    self.reads.append((a, b - a))
                    return b"".join(parts)
        raise IOError("read outside fetched ranges: %d-%d" % (a, b))

    def close(self):
        self.closed = True
