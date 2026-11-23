"""HTTP range access to Hugging Face dataset files (Tahoe-100M), one connection per
caller, with redirect caching, retry with exponential backoff on HTTP 429/5xx and
connection errors, and an in-memory range cache. Nothing read here is written to disk."""
import io
import time
import random
import threading
import requests

BASE = "https://huggingface.co/datasets/tahoebio/Tahoe-100M/resolve/349e2a55f5538a3e73b32a06c6bdd8e9e87011e8/"
_tls = threading.local()
STATS = {"requests": 0, "bytes": 0, "retries": 0, "http429": 0}
_lock = threading.Lock()


def session():
    s = getattr(_tls, "s", None)
    if s is None:
        s = requests.Session()
        s.headers["User-Agent"] = "p59-analysis-client/1.0"
        _tls.s = s
        _tls.redir = {}
    return s


def get_range(path, start, end_incl):
    """Return bytes [start, end_incl] of a dataset file."""
    for attempt in range(12):
        url = BASE + path + f"?range={start}-{end_incl}"
        try:
            r = session().get(url, headers={"Range": f"bytes={start}-{end_incl}"}, timeout=180)
            if r.status_code == 206:
                b = r.content
                if (len(b) != end_incl - start + 1 or
                        not r.headers.get("Content-Range", "").startswith(f"bytes {start}-{end_incl}/")):
                    raise requests.RequestException("short read")
                with _lock:
                    STATS["requests"] += 1
                    STATS["bytes"] += len(b)
                return b
            if r.status_code in (403, 404, 410) and attempt < 3:
                continue  # expired signed URL: re-resolve
            if r.status_code == 429 or r.status_code >= 500:
                with _lock:
                    STATS["retries"] += 1
                    STATS["http429"] += r.status_code == 429
                time.sleep(min(180, 2 ** attempt + random.random()))
                continue
            raise RuntimeError(f"HTTP {r.status_code} for {path} {start}-{end_incl}")
        except requests.RequestException:
            with _lock:
                STATS["retries"] += 1
            time.sleep(min(180, 2 ** attempt + random.random()))
    raise RuntimeError(f"giving up on {path} {start}-{end_incl}")


class RangeFile(io.RawIOBase):
    """Seekable read-only file over HTTP ranges with an in-memory cache of fetched blocks."""

    def __init__(self, path, size, tail=1 << 20):
        self.path, self.size, self.pos = path, size, 0
        self.blocks = []  # list of (start, bytes)
        self.misses = 0
        t0 = max(0, size - tail)
        self.add(t0, get_range(path, t0, size - 1))

    def add(self, start, data):
        self.blocks.append((start, data))

    def prefetch(self, ranges, merge_gap=1 << 16):
        """ranges: list of (start, length). Merge near ranges and fetch them."""
        rs = sorted((s, s + n) for s, n in ranges if n > 0)
        merged = []
        for s, e in rs:
            if merged and s - merged[-1][1] <= merge_gap:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        for s, e in merged:
            if self._find(s, e - s) is None:
                self.add(s, get_range(self.path, s, e - 1))
        return len(merged)

    def drop(self):
        n = sum(len(b) for _, b in self.blocks)
        self.blocks = []
        return n

    def _find(self, s, n):
        for bs, b in self.blocks:
            if bs <= s and s + n <= bs + len(b):
                return b[s - bs:s - bs + n]
        return None

    def readable(self):
        return True

    def seekable(self):
        return True

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

    def readinto(self, buf):
        n = min(len(buf), self.size - self.pos)
        if n <= 0:
            return 0
        d = self._find(self.pos, n)
        if d is None:
            self.misses += 1
            d = get_range(self.path, self.pos, self.pos + n - 1)
            self.add(self.pos, d)
        buf[:n] = d
        self.pos += n
        return n
