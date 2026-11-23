"""Load the per-line sparse vehicle-cell counts written by s03_reduce.py."""
import gzip
from pathlib import Path
import numpy as np, scipy.sparse as sps

D = "data/F2/sparse"
G = 62710


def read_bytes(path):
    file = Path(path)
    if file.exists():
        return file.read_bytes()
    with gzip.open(str(file) + ".gz", "rb") as stream:
        return stream.read()


def load_line(li):
    prefix = f"{D}/{li:02d}."
    nnz = np.frombuffer(read_bytes(prefix + "nnz.u32"), dtype=np.uint32).astype(np.int64)
    genes = np.frombuffer(read_bytes(prefix + "genes.u16"), dtype=np.uint16)
    counts = np.frombuffer(read_bytes(prefix + "counts.u16"), dtype=np.uint16)
    bcs = read_bytes(prefix + "barcodes.txt").decode().splitlines()
    assert len(bcs) == len(nnz) and nnz.sum() == len(genes) == len(counts), (li, len(bcs), len(nnz), nnz.sum(), len(genes))
    indptr = np.concatenate([[0], np.cumsum(nnz)])
    csr = sps.csr_matrix((counts.astype(np.float32), genes.astype(np.int32), indptr), shape=(len(nnz), G))
    return dict(barcodes=np.array(bcs), csr=csr)
