from pathlib import Path
import pandas as pd
b=Path(__file__).resolve().parent
keys=['drug','lineage','tissue']
a=pd.read_csv(b/'comparison/phase5_tested_comparison.csv')
c=pd.read_csv(b/'sensitivity_gseapy_1.3.1/comparison/phase5_tested_comparison.csv')
cols=['NES_replay','FDR_replay','category_replay']
joined=a.merge(c[keys+cols],on=keys,suffixes=('_1.2.1','_1.3.1'),validate='one_to_one')
assert len(joined)==620
joined=joined.rename(columns={'NES_old':'NES_archive','FDR_old':'FDR_archive','category_old':'category_archive'})
joined.to_csv(b/'comparison/phase5_software_versions.csv',index=False)
