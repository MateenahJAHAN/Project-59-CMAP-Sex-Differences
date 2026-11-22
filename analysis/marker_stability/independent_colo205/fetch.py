"""Fetch missing inputs and confirm the recorded SHA-256 for every source."""
import csv
from hashlib import file_digest
from pathlib import Path
import requests

ROOT = Path(__file__).resolve().parent

def check_sources():
    dest=ROOT/"inputs"
    dest.mkdir(exist_ok=True)
    for row in csv.DictReader((ROOT/"sources.tsv").open(), delimiter="\t"):
        p=dest/row["file"]
        if not p.exists():
            tmp=p.with_suffix(p.suffix+".tmp")
            with requests.get(row["url"],stream=True,timeout=(25,180)) as r:
                r.raise_for_status()
                with tmp.open("wb") as f:
                    for block in r.iter_content(2**20):
                        f.write(block)
            with tmp.open("rb") as f:
                if file_digest(f,"sha256").hexdigest()!=row["sha256"]:
                    raise ValueError("Source has changed since this analysis: "+row["file"])
            tmp.replace(p)
        with p.open("rb") as f:
            if file_digest(f,"sha256").hexdigest()!=row["sha256"]:
                raise ValueError("Checksum mismatch: "+str(p))
if __name__=="__main__":
    check_sources()
