"""Read both raw Infinium colors and compare COLO 205 to the same-study controls."""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from analyze import IN, OUT, write
sys.path.insert(0,str(Path(__file__).resolve().parent/"nci_dna"))
from check_idats import read_means
from analyse import probe_manifest

def gdsc():
    ann=probe_manifest(IN)
    cols=["probeID","address_A","address_B","designType","channel"]
    chip=pd.read_csv(IN/"manifest.gz",sep="\t",usecols=cols,low_memory=False).set_index("probeID").loc[ann.index]
    rows=[]; probes=[]
    for gsm,cell in [("GSM1669683","COLO205"),("GSM1669868","HCT116"),("GSM1669869","HCT15"),("GSM1669913","HT29")]:
        means={c:read_means(next(IN.glob(f"{gsm}*{c}.idat.gz"))) for c in ["Red","Grn"]}
        total=chip.address_A.map(means["Grn"])+chip.address_A.map(means["Red"])
        for color in ["Red","Grn"]:
            idx=(chip.designType=="I")&(chip.channel==color)
            total.loc[idx]=chip.loc[idx,"address_A"].map(means[color])+chip.loc[idx,"address_B"].map(means[color])
        assert total.notna().all()
        sig=np.log2(1+total)
        base=sig[ann.arm=="auto"].median()
        for arm in ["Yp","Yq"]:
            v=sig[ann.arm==arm]
            rows.append(dict(cell=cell,accession=gsm,arm=arm,n_probes=len(v),autosome_median=base,arm_median=v.median(),centered=v.median()-base))
            co=ann.loc[v.index,["arm","CpG_beg"]].copy()
            co["gsm"]=gsm;co["cell"]=cell;co["centered"]=v-base
            probes.append(co.reset_index())
    write(rows,"gdsc_arm_medians")
    write(pd.concat(probes),"gdsc_Y_probes")

if __name__=="__main__":
    gdsc()
