"""Compare the exact GTEx v8 gene-set table against the member of the provider release."""
import argparse
import hashlib,io,json,tarfile
from pathlib import Path
import requests
root=Path(__file__).resolve().parent
urls=['https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/sb-eqtl/GTEx_Analysis_v8_sbgenes.tar.gz',
      'https://storage.googleapis.com/gtex_analysis_v8/sex_biased_genes_and_sbeqtl_data/GTEx_Analysis_v8_sbgenes.tar.gz']
for url in urls:
    response=requests.get(url,timeout=100)
    if response.ok:break
else:response.raise_for_status()
with tarfile.open(fileobj=io.BytesIO(response.content),mode='r:gz') as t:
    hits=[n for n in t.getnames() if n.split('/')[-1]=='signif.sbgenes.txt']
    assert len(hits)==1,hits
    raw=t.extractfile(hits[0]).read()
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('gtex_table', type=Path, help='Path to the local signif.sbgenes.txt file')
local = parser.parse_args().gtex_table.read_bytes()
report={'url':url,'archive_sha256':hashlib.sha256(response.content).hexdigest(),'archive_bytes':len(response.content),
        'member':hits[0],'member_sha256':hashlib.sha256(raw).hexdigest(),'member_bytes':len(raw),
        'table_identical':local==raw}
(root/'gtex_comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print(report)
