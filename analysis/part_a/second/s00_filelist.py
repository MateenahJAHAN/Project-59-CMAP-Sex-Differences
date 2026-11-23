"""Step 0: list Tahoe-100M shard and metadata files with LFS SHA-256 (HF tree API, paginated)."""
import json, re, sys, time, requests
import pandas as pd
OUT = sys.argv[1] if len(sys.argv) > 1 else "results/F2/F2_filelist.tsv"
rows = []
for sub in ["data", "metadata"]:
    url = f"https://huggingface.co/api/datasets/tahoebio/Tahoe-100M/tree/349e2a55f5538a3e73b32a06c6bdd8e9e87011e8/{sub}"
    while url:
        r = requests.get(url, timeout=120); r.raise_for_status()
        for x in r.json():
            if x["type"] == "file":
                rows.append(dict(path=x["path"], size=x["size"], lfs_sha256=(x.get("lfs") or {}).get("oid", "")))
        m = re.search(r'<([^>]+)>;\s*rel="next"', r.headers.get("link", ""))
        url = m.group(1) if m else None
df = pd.DataFrame(rows)
shards = df[df.path.str.match(r"^data/train-\d{5}-of-03388\.parquet$")]
assert len(shards) == 3388, len(shards)
assert shards.path.is_unique
print("shards", len(shards), "bytes", int(shards["size"].sum()))
df.to_csv(OUT, sep="\t", index=False)
