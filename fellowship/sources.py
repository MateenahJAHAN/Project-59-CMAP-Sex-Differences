"""Acquire the exact public LINCS and GTEx inputs used in the fellowship screens."""
import csv
import gzip
import hashlib
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil
import tarfile
import tempfile
import time
from urllib.request import Request, urlopen

REPO = Path(__file__).resolve().parents[1]
PROOF = REPO / 'verification/source_replay'
SOURCES = json.loads((PROOF / 'phases2_5/source_manifest.json').read_text())
FILES = {r['filename']: r['sha256'] for r in SOURCES.values() if isinstance(r, dict)}
MATRIX = json.loads((PROOF / 'lincs_inputs/comparison.json').read_text())
GTEX = json.loads((PROOF / 'lincs_inputs/gtex_comparison.json').read_text())
LINCS = 'https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/'


def digest(path):
    with path.open('rb') as file:
        return hashlib.file_digest(file, 'sha256').hexdigest()


def required(path):
    if not path.exists():
        return False
    if digest(path) != FILES[path.name]:
        raise ValueError(f'Checksum mismatch for {path}. Move the mismatched file and rerun.')
    return True


def install(path, writer):
    """A failed fetch never leaves an input with a final filename."""
    with tempfile.TemporaryDirectory(dir=path.parent) as dirname:
        tmp = Path(dirname) / path.name
        writer(tmp)
        required(tmp)
        os.replace(tmp, path)


def download(url, dest):
    with urlopen(url, timeout=120) as stream, dest.open('wb') as out:
        shutil.copyfileobj(stream, out, 1024 * 1024)


def fetch_metadata(data):
    for name, archived in [('cellinfo_beta.txt', 'cellinfo_beta (2).txt'),
                           ('geneinfo_beta.txt', 'geneinfo_beta (1).txt')]:
        if not required(data / name):
            install(data / name, lambda tmp: shutil.copyfile(REPO / 'cmap' / archived, tmp))
    if not required(data / 'siginfo_beta.txt'):
        install(data / 'siginfo_beta.txt', lambda tmp: download(LINCS + 'siginfo_beta.txt', tmp))
    if not required(data / 'signif.sbgenes.txt'):
        archive = data / 'GTEx_Analysis_v8_sbgenes.tar.gz'
        if not required(archive):
            install(archive, lambda tmp: download(GTEX['url'], tmp))
        def extract(tmp):
            with tarfile.open(archive) as tar, tar.extractfile(GTEX['member']) as src, tmp.open('wb') as out:
                shutil.copyfileobj(src, out)
        install(data / 'signif.sbgenes.txt', extract)


def fetch_matrix(data):
    if all([required(data / name) for name in ['row_ids.txt', 'sig50.sigs.txt', 'sig50.npy']]):
        return
    import h5py
    import fsspec
    import numpy as np

    with gzip.open(PROOF / 'lincs_inputs/signature_sha256.tsv.gz', 'rt') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    indices = [int(r['provider_row']) for r in rows]
    assert len(indices) == len(set(indices)) == 36523 and indices == sorted(indices)
    url = MATRIX['source']
    info = urlopen(Request(url, method='HEAD'), timeout=60).headers
    if info['ETag'] != MATRIX['source_etag'] or int(info['Content-Length']) != MATRIX['source_file_bytes']:
        raise ValueError('The provider matrix changed. Check the recorded version before using it.')
    with fsspec.open(url, block_size=4*1024**2) as stream, h5py.File(stream, 'r') as hdf:
        columns = hdf['0/META/COL/id'][:]
        genes = hdf['0/META/ROW/id'][:]
        mat = hdf['0/DATA/0/matrix']
        assert mat.chunks is None and mat.dtype == np.dtype('<f4')
        assert list(mat.shape) == MATRIX['provider_shape'] and mat.id.get_offset() == MATRIX['matrix_offset']
    if any(columns[i].decode() != r['sig_id'] for i, r in zip(indices, rows)):
        raise ValueError('Source signature order differs from the study')
    install(data/'row_ids.txt', lambda p: p.write_bytes(b'\n'.join(genes) + b'\n'))
    install(data/'sig50.sigs.txt', lambda p: p.write_bytes(('\n'.join(r['sig_id'] for r in rows)+'\n').encode()))
    stride = len(genes) * 4
    groups = []
    for pos, i in enumerate(indices):
        if not groups or i-indices[groups[-1][0]] >= 128 or i-indices[groups[-1][-1]] > 7:
            groups.append([])
        groups[-1].append(pos)

    def build(tmp):
        arr = np.lib.format.open_memmap(tmp, mode='w+', dtype='<f4', shape=(len(indices), len(genes)))
        # Use persistent connections for adjacent byte ranges. Each worker writes disjoint rows.
        import requests
        import threading
        local = threading.local()
        def fetch(group):
            first, last = indices[group[0]], indices[group[-1]]
            lo = MATRIX['matrix_offset'] + first * stride
            hi = MATRIX['matrix_offset'] + (last+1) * stride - 1
            for attempt in range(5):
                try:
                    if not hasattr(local, 'http'):
                        local.http = requests.Session()
                    with local.http.get(url, headers={'Range':f'bytes={lo}-{hi}', 'If-Match':MATRIX['source_etag']}, stream=True, timeout=(20,90)) as resp:
                        if resp.status_code != 206 or resp.headers.get('Content-Range') != f"bytes {lo}-{hi}/{MATRIX['source_file_bytes']}":
                            raise IOError('Provider did not return the requested byte range')
                        payload = resp.content
                    if len(payload) != hi-lo+1:
                        raise IOError('Truncated matrix response')
                    for pos in group:
                        start = (indices[pos]-first)*stride
                        row = memoryview(payload)[start:start+stride]
                        if hashlib.sha256(row).hexdigest() != rows[pos]['provider_float32_sha256']:
                            raise IOError(f'Provider checksum mismatch for signature {rows[pos]["sig_id"]}')
                        arr[pos] = np.frombuffer(row, dtype='<f4')
                    return
                except (requests.RequestException,OSError):
                    if attempt == 4:
                        raise
                    time.sleep(2**attempt)
        with ThreadPoolExecutor(max_workers=12) as pool:
            for i, _ in enumerate(pool.map(fetch, groups), 1):
                if i % 1000 == 0:
                    print(f'Fetched {i}/{len(groups)} matrix ranges', flush=True)
        arr.flush()
    install(data/'sig50.npy', build)


def prepare(data):
    data.mkdir(parents=True, exist_ok=True)
    fetch_metadata(data)
    fetch_matrix(data)
