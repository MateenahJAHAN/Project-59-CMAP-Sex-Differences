"""Rebuild well-level marker counts directly from the fourteen sealed Tahoe cell tables."""
import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

OUT = Path(__file__).resolve().parent
MARKERS = ["XIST","RPS4Y1","DDX3Y","KDM5D","UTY","EIF1AY","ZFY","USP9Y","NLGN4Y"]
KEYS = ["plate","sample","drug","cell_line_id"]
SUMS = ["total_umi"] + MARKERS

def run(root):
    manifest = {line.split()[1]: line.split()[0] for line in (root/"stream_SHA256SUMS.txt").read_text().splitlines() if len(line.split()) == 2}
    groups, sources = [], []
    for plate in range(1,15):
        rel = f"plates/plate{plate}/cells_markers.parquet"
        source = root/"stream"/rel
        with source.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        assert digest == manifest[rel]
        parquet = pq.ParquetFile(source)
        chunks = []
        for batch in parquet.iter_batches(columns=KEYS + SUMS, batch_size=400000):
            counts = batch.to_pandas().groupby(KEYS, sort=False, observed=True)
            sums = counts[SUMS].sum()
            sums["cells"] = counts.size()
            chunks.append(sums)
        merged = pd.concat(chunks).groupby(level=list(range(len(KEYS)))).sum()
        expected = json.loads((source.parent/"MERGE_OK.json").read_text())["cells"]
        assert merged.cells.sum() == parquet.metadata.num_rows == expected
        groups.append(merged)
        sources.append({"table":rel, "sha256":digest, "rows":expected})
        print(f"plate {plate}: {expected:,} cells",flush=True)
    grouped = pd.concat(groups).reset_index()
    assert grouped.cells.sum() == 95624334 and not grouped.duplicated(KEYS).any()
    grouped.to_csv(OUT/"well_counts.tsv.gz", sep="\t", index=False, compression={"method":"gzip","mtime":0})
    (OUT/"cell_table_manifest.json").write_text(json.dumps(sources,indent=2)+"\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=Path.home()/"p59_data/tahoe")
    run(parser.parse_args().root)
