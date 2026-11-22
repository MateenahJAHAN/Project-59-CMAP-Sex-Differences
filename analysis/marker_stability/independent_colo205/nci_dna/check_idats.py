"""Recover total signal from both color IDATs, including probes omitted from the series table."""
import argparse
import gzip
import hashlib
from pathlib import Path
import struct

import numpy as np
import pandas as pd
import requests

from analyse import OUT, probe_manifest, write_plots

SAMPLES = [
    ("GSE79185", "COLO205", "GSM2087430", "9376537080_R05C01"),
    ("GSE86078", "COLO205", "GSM2293091", "200486890054_R01C02"),
    ("GSE86078", "HCC2998", "GSM2293106", "3998426003_R05C02"),
    ("GSE86078", "HCT116", "GSM2293107", "200486890054_R04C02"),
    ("GSE86078", "HCT15", "GSM2293108", "200486890054_R05C02"),
    ("GSE86078", "HT29", "GSM2293128", "200693590019_R01C01"),
    ("GSE86078", "SW620", "GSM2293194", "200693590019_R06C02"),
]


def read_means(path):
    raw = gzip.decompress(path.read_bytes())
    assert raw[:4] == b"IDAT" and struct.unpack_from("<Q",raw,4)[0] == 3
    fields = struct.unpack_from("<I",raw,12)[0]
    offsets = dict(struct.unpack_from("<HQ",raw,16+i*10) for i in range(fields))
    n = struct.unpack_from("<I",raw,offsets[1000])[0]
    addresses = np.frombuffer(raw,"<i4",n,offsets[102])
    assert len(set(addresses)) == n
    return pd.Series(np.frombuffer(raw,"<u2",n,offsets[104]).astype(float),index=addresses)


def main(cache):
    ann = probe_manifest(cache)
    cols = ["probeID","address_A","address_B","designType","channel"]
    chip = pd.read_csv(cache/"manifest.gz",sep="\t",usecols=cols,low_memory=False).set_index("probeID").loc[ann.index]
    matrices = pd.read_csv(OUT/"y_probe_intensities.tsv.gz",sep="\t")
    centers = pd.read_csv(OUT/"autosomal_reference.tsv",sep="\t").set_index(["study","cell_key"])
    rows, meta, files = [], [], []
    for study,cell,gsm,chipid in SAMPLES:
        means={}
        for color in ("Red","Grn"):
            name=f"{gsm}_{chipid}_{color}.idat.gz"
            url=f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{name}"
            path=cache/name
            if not path.exists():
                r=requests.get(url,timeout=60)
                r.raise_for_status()
                path.write_bytes(r.content)
            files.append(dict(name=name,url=url,sha256=hashlib.sha256(path.read_bytes()).hexdigest()))
            means[color]=read_means(path)
        # Infinium II uses both channels at A; Infinium I uses A + B in its designated channel.
        total=pd.Series(index=chip.index,dtype=float)
        for channel, idx in chip.groupby("channel").groups.items():
            probes=chip.loc[idx]
            if channel=="Both":
                total.loc[idx]=probes.address_A.map(means["Red"])+probes.address_A.map(means["Grn"])
            else:
                total.loc[idx]=probes.address_A.map(means[channel])+probes.address_B.map(means[channel])
        assert total.notna().all()
        intensity=np.log2(total+1)
        autosomal_center=intensity[ann.arm=="auto"].median()
        values=ann.loc[ann.arm!="auto",["arm","CpG_beg"]].copy()
        values["log2_total"]=intensity
        values["log2_centered"]=intensity-autosomal_center
        values["study"]=study
        values["gsm"]=gsm
        values["cell_key"]=cell
        values=values.join(matrices.loc[(matrices.study==study) & (matrices.cell_key==cell)].set_index("probeID").log2_centered.rename("submitted_centered"))
        comparable=values.submitted_centered.notna()
        delta=(values.loc[comparable,"log2_total"] -
               (values.loc[comparable,"submitted_centered"] + centers.loc[(study,cell),"autosomal_median"]))
        assert delta.abs().max() < 2e-6, (gsm, delta.abs().max())
        for arm,sub in values.groupby("arm"):
            meta.append(dict(study=study,cell=cell,accession=gsm,arm=arm,n_total=len(sub),n_deposited=int(sub.submitted_centered.count()),
                             raw_idat_median=sub.log2_centered.median(),autosomal_center=autosomal_center,
                             deposited_median=sub.submitted_centered.median(),max_abs_raw_log_difference_on_published_probes=float(delta.loc[sub.index.intersection(delta.index)].abs().max())))
        rows.append(values.reset_index())
    raw = pd.concat(rows,ignore_index=True)
    raw.to_csv(OUT/"idat_probes.tsv",sep="\t",index=False,float_format="%.7f")
    pd.DataFrame(files).to_csv(OUT/"idat_sources.tsv",sep="\t",index=False)
    pd.DataFrame(meta).to_csv(OUT/"idat_summary.tsv",sep="\t",index=False,float_format="%.7f")
    lookup = matrices[["study","cell_key","reference"]].drop_duplicates()
    raw = raw.merge(lookup, on=["study","cell_key"], validate="many_to_one")
    combined = pd.concat([matrices[matrices.study=="GSE79185"],raw[raw.study=="GSE86078"]],ignore_index=True)
    assert combined.log2_centered.notna().all()
    stats = combined.groupby(["study","cell_key","reference","arm"]).log2_centered.agg(["count","median"]).reset_index()
    stats.to_csv(OUT/"complete_arm_medians.tsv",sep="\t",index=False,float_format="%.7f")
    write_plots(combined,stats,"complete_dna_intensity")
    print(pd.DataFrame(meta).round(5).to_string(index=False))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("cache",type=Path)
    main(p.parse_args().cache)
