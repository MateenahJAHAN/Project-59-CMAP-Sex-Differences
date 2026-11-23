"""Step 1a: vehicle (DMSO_TF) rows of obs_metadata. Only rows with drug == 'DMSO_TF' are kept;
other rows are dropped per row group right after reading and never written."""
import sys, time
import pyarrow.parquet as pq, pyarrow.compute as pc, pyarrow as pa
import pandas as pd
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import runlog
t0 = time.time()
pf = pq.ParquetFile("data/F2/meta/obs_metadata.parquet")
cols = ["plate", "BARCODE_SUB_LIB_ID", "sample", "gene_count", "tscp_count", "pcnt_mito", "drug",
        "cell_line", "cell_name", "phase", "pass_filter"]
parts, n_all = [], 0
drug_vals = set()
for i in range(pf.metadata.num_row_groups):
    t = pf.read_row_group(i, columns=cols)
    n_all += t.num_rows
    d = t.column("drug")
    drug_vals.update(v for v in pc.unique(d).to_pylist() if v and "DMSO" in v)
    t = t.filter(pc.equal(d, "DMSO_TF"))
    parts.append(t.to_pandas())
    del t, d
obs = pd.concat(parts, ignore_index=True)
for c in ["cell_line", "cell_name", "phase", "pass_filter"]:
    obs[c] = obs[c].astype(str)
assert n_all == 100648790, n_all
assert (obs.drug == "DMSO_TF").all()
print("drug values containing 'DMSO':", sorted(drug_vals))
print("vehicle rows (all pass_filter):", len(obs))
print(obs.pass_filter.value_counts())
print("phase values:", obs.phase.value_counts().to_dict())
obs.to_parquet("data/F2/obs_dmso.parquet", index=False)
runlog("s01_obs", t0, f"obs rows {n_all}; DMSO_TF rows {len(obs)}")
