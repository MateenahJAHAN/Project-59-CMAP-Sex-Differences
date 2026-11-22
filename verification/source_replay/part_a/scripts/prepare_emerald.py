"""Prepare the three pinned EmeraldBay metadata files and all 116 original shards' source checksums."""
import argparse
import csv
import json
from pathlib import Path

import sys
sys.dont_write_bytecode=True
sys.path.insert(0,str(Path(__file__).resolve().parent))
from fetch import download

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('run', type=Path, help='scratch directory for this run')
args = parser.parse_args()
base = Path(__file__).resolve().parents[1]
spec = json.loads((base/'receipts/emerald_metadata.json').read_text())
root = args.run.resolve()
folder = root/'data/F1/meta'
folder.mkdir(parents=True, exist_ok=True)
(root/'results/F1').mkdir(parents=True, exist_ok=True)
(root/'provenance').mkdir(parents=True, exist_ok=True)
url = 'https://huggingface.co/datasets/tahoebio/EmeraldBay/resolve/'+spec['revision']+'/metadata/'
for name, check in spec['metadata'].items():
    download(url+name, str(folder/name), check['sha256'], check['size'])
with (base/'receipts/emerald_all_condition_coverage.tsv').open(newline='') as file:
    coverage = list(csv.DictReader(file, delimiter='\t'))
assert len(coverage) == len({row['shard'] for row in coverage}) == 116
tree = [{'path':'expression_data/'+r['shard'],'lfs':{'oid':r['sha256'],'size':int(r['bytes'])}} for r in coverage]
(folder/'tree_expression.json').write_text(json.dumps(tree)+'\n')
print(f'Prepared {len(tree)} pinned shards and {len(spec["metadata"])} metadata files in {root}')
