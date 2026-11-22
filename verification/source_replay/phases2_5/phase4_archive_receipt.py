"""Compare each original notebook run log with each historical CSV at drug–lineage level."""
import hashlib,json,re
from pathlib import Path
import pandas as pd
import compare

out=Path(__file__).resolve().parent/'comparison'
src=next(compare.REPO.glob('Phase4_*'))
keys=['drug','lineage','tissue','sex_direction']
csvs=sorted(src.glob('*.csv'))
tables=[pd.read_csv(p).sort_values(keys).reset_index(drop=True) for p in csvs]
assert tables[0].equals(tables[1])
ref=tables[0]
ref['significant']=ref.category!='C'
truth=ref.groupby(keys[:2]).significant.sum()
results=[]
for n in sorted(src.glob('*.ipynb')):
    cells=json.loads(n.read_text())['cells']
    k=next(i for i,c in enumerate(cells) if 'all_classifications = []' in ''.join(c['source']))
    log='\n'.join(''.join(o.get('text',[])) for o in cells[k]['outputs'])
    pat=r'\[\s*\d+/50\]\s+(.*?)\s+×\s+(.*?)\s+→\s+\d+s · (\d+) significant rows'
    rows=[(d.strip(),l.strip(),int(count)) for d,l,count in re.findall(pat,log)]
    captured=pd.DataFrame(rows,columns=['drug','lineage','count']).set_index(keys[:2])['count'].sort_index()
    assert len(captured)==captured.index.nunique() and captured.index.equals(truth.index)
    srcfn=''.join(cells[49]['source']); pipelinesrc=''.join(cells[k]['source'])
    results.append({'notebook':n.name,'source_cells':[49,k],'function_cell_sha256':hashlib.sha256(srcfn.encode()).hexdigest(),'batch_cell_sha256':hashlib.sha256(pipelinesrc.encode()).hexdigest(),'logged_pairs':len(captured),'pair_counts_different_from_either_csv':int((captured!=truth).sum())})
result={'both_csvs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in csvs},'rows_each':len(ref),'columns':list(tables[1]),'all_keyed_rows_and_all_fields_equal':True,'captured_notebook_runs':results,'batch_logic':'Cell 51 initializes an empty classification list and calls the analysis and classification functions from cell 49 once per pair. Classification takes NES, FDR, set name and category together from a single significant GSEA row. Cell 52 reads the completed CSV for plotting and copies it; it does not merge numeric columns.'}
(out/'phase4_historical_exports.json').write_text(json.dumps(result,indent=2)+'\n')
print(result,flush=True)
