"""Preserve the sample identities and fixed vehicle definitions used in the analysis."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import pandas as pd

OUT = Path(__file__).resolve().parent
p = argparse.ArgumentParser()
p.add_argument("--root",type=Path,default=Path.home()/"p59_data/tahoe")
root = p.parse_args().root
cutoffs = root/"partA_local/partA_local_cutoffs.tsv"
states = root/"partA_local/partA_local_states.tsv"
samples = root/"index/meta/sample_metadata.parquet"
manifest = {str(f.relative_to(root)): hashlib.sha256(f.read_bytes()).hexdigest() for f in [cutoffs,states,samples]}
df = pd.read_csv(states, sep="\t")
df[["cellosaurus_id","cell_name","annotation","organ","xist_high","y_high","confidence"]].to_csv(OUT/"reference_lines.tsv",sep="\t",index=False)
pd.read_csv(cutoffs, sep="\t")[["score","cutoff"]].to_csv(OUT/"reference_cutoffs.tsv",sep="\t",index=False)
meta = pd.read_parquet(samples)
entry = [ast.literal_eval(s) for s in meta.drugname_drugconc]
assert all(len(row)==1 and row[0][2]=="uM" and row[0][0]==name for row,name in zip(entry,meta.drug))
meta["dose_um"] = [row[0][1] for row in entry]
meta[["sample","plate","drug","dose_um"]].to_csv(OUT/"sample_key.tsv",sep="\t",index=False)
(OUT/"source_tables.json").write_text(json.dumps(manifest,indent=2)+"\n")
