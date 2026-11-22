"""Download helpers: streamed download with SHA-256 check and provenance logging."""
import hashlib, os, time, datetime, fcntl
import requests

def log_download(root, url, local_path, nbytes, sha, notes):
    p = os.path.join(root, "provenance/downloads.tsv")
    ts = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    rel = os.path.relpath(local_path, root)
    with open(p, "a") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        fh.write(f"F1\t{url}\t{rel}\t{nbytes}\t{sha}\t{ts}\t{notes}\n")
        fcntl.flock(fh, fcntl.LOCK_UN)

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()

def download(url, path, expected_sha=None, expected_size=None, retries=6):
    """Stream url to path. Returns (bytes, sha256). Raises if checksum fails after retries."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path) and expected_sha and os.path.getsize(path) == expected_size:
        s = sha256_file(path)
        if s == expected_sha:
            return os.path.getsize(path), s, "cached"
    tmp = path + ".part"
    for attempt in range(retries):
        try:
            h = hashlib.sha256(); n = 0
            with requests.get(url, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(tmp, "wb") as fh:
                    for chunk in r.iter_content(1 << 22):
                        fh.write(chunk); h.update(chunk); n += len(chunk)
            s = h.hexdigest()
            if expected_size is not None and n != expected_size:
                raise IOError(f"size mismatch {n} vs {expected_size}")
            if expected_sha is not None and s != expected_sha:
                raise IOError(f"sha mismatch {s} vs {expected_sha}")
            os.replace(tmp, path)
            return n, s, "downloaded"
        except Exception as ex:
            print(f"download retry {attempt} {url}: {ex}", flush=True)
            time.sleep(10 * (attempt + 1))
    raise RuntimeError(f"failed to download {url}")
