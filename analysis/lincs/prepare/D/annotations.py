"""Rebuild landmark chromosomes from Ensembl 110 and historical drug overlap from the Phase 5 selection."""
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sys

import pandas as pd

work = Path(sys.argv[1]).resolve()
here = Path(__file__).resolve().parents[2]
repo = here.parents[1]
out = work / 'results/references'
out.mkdir(parents=True, exist_ok=True)
lm_path = work / 'results/D/D_landmarks.tsv'
gi_path = work / 'sources/geneinfo_beta.txt'
gtf_path = work / 'sources/Homo_sapiens.GRCh38.110.gtf.gz'

# ---------- 1. landmark chromosomes ----------
lm = pd.read_csv(lm_path, sep="\t", dtype={"gene_id": str})
assert len(lm) == 978 and lm.gene_id.is_unique and lm.gene_symbol.is_unique
gi = pd.read_csv(gi_path, sep="\t", dtype=str, usecols=["gene_id", "gene_symbol", "ensembl_id", "feature_space"])
assert gi.gene_id.is_unique
gi_lm = gi[gi.feature_space == "landmark"]
assert len(gi_lm) == 978 and set(gi_lm.gene_id) == set(lm.gene_id)
chrom_by_id, ids_by_name = {}, {}
with gzip.open(gtf_path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if f[2] != "gene":
            continue
        attr = f[8]
        gid = re.search(r'gene_id "([^"]+)"', attr).group(1)
        m = re.search(r'gene_name "([^"]+)"', attr)
        assert gid not in chrom_by_id
        chrom_by_id[gid] = f[0]
        if m:
            ids_by_name.setdefault(m.group(1), []).append(gid)
t = lm.merge(gi_lm[["gene_id", "ensembl_id"]], on="gene_id", how="left", validate="1:1")
assert len(t) == 978
chrom, key = [], []
for _, r in t.iterrows():
    e = r.ensembl_id if isinstance(r.ensembl_id, str) else ""
    if e in chrom_by_id:
        chrom.append(chrom_by_id[e]); key.append("geneinfo ensembl_id")
    elif len(ids_by_name.get(r.gene_symbol, [])) == 1:
        chrom.append(chrom_by_id[ids_by_name[r.gene_symbol][0]]); key.append("unique GTF gene_name")
    else:
        chrom.append("unmatched"); key.append("none")
t["chrom"] = chrom
t["match_key"] = key
t["is_X"] = t.chrom == "X"
t = t[["gene_id", "gene_symbol", "ensembl_id", "chrom", "is_X", "match_key"]]
t.to_csv(os.path.join(out, "landmark_chromosomes.tsv"), sep="\t", index=False)
print("landmarks:", len(t), "match keys:", t.match_key.value_counts().to_dict(), "X-linked:", int(t.is_X.sum()),
      "Y:", int((t.chrom == "Y").sum()), "unmatched:", int((t.chrom == "unmatched").sum()))


legacy = set(json.loads((repo / 'Phase5_PerSignature_Categories/reference/top50_names.json').read_text()))
panel = pd.read_csv(work / 'results/D/D_panel_definitions.tsv', sep='\t')
eligible = panel.included & panel.cmap_name.isin(legacy) & panel.panel.isin(['P1', 'P2'])
panel.loc[eligible, ['panel', 'cmap_name']].sort_values(['panel','cmap_name']).to_csv(out/'old_list_overlap.tsv', sep='\t', index=False)
for name in ('landmark_chromosomes.tsv', 'old_list_overlap.tsv'):
    reference = here / 'inputs/G' / name
    if (out / name).read_bytes() != reference.read_bytes():
        raise ValueError(f'Rebuilt annotation differs: {name}')
    print(name, hashlib.sha256(reference.read_bytes()).hexdigest(), flush=True)
