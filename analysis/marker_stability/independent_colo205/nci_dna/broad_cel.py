"""Extract chromosome Y nonpolymorphic probe intensities from three CCLE SNP 6.0 CEL files."""
import argparse
import gzip
import hashlib
from pathlib import Path
import shutil
import sqlite3
import struct
import tarfile

import numpy as np
import pandas as pd
import requests

OUT=Path(__file__).resolve().parent
CHIP_URL="https://bioconductor.org/packages/3.23/data/annotation/src/contrib/pd.genomewidesnp.6_3.14.1.tar.gz"
DB="pd.genomewidesnp.6/inst/extdata/pd.genomewidesnp.6.sqlite"
CELLS={"COLO205":"GSM888007","HCT15":"GSM888133","HT29":"GSM888212"}


def obtain(path,url):
    if not path.exists():
        tmp=path.with_suffix(".part")
        with requests.get(url,stream=True,timeout=120) as r, tmp.open("wb") as f:
            r.raise_for_status()
            r.raw.decode_content = path.suffix == ".txt"
            shutil.copyfileobj(r.raw,f)
        tmp.replace(path)
    with path.open("rb") as f:
        return {"file":path.name,"url":url,"sha256":hashlib.file_digest(f,"sha256").hexdigest()}


def cel(path):
    with gzip.open(path,"rb") as f:
        magic,version,nrows,ncols,n,h_size=struct.unpack("<6i",f.read(24))
        assert (magic,version)==(64,4) and n==nrows*ncols
        header=f.read(h_size)
        assert f"Rows={nrows}\n".encode() in header and f"Cols={ncols}\n".encode() in header
        for _ in range(2):
            f.read(struct.unpack("<I",f.read(4))[0])
        _,outliers,masks,subgrids=struct.unpack("<iIIi",f.read(16))
        assert subgrids==0
        records=np.frombuffer(f.read(n*10),dtype=[("mean","<f4"),("sd","<f4"),("pixels","<i2")])
        assert len(records)==n and (records["mean"]>=0).all()
        flags=np.frombuffer(f.read(4*(masks+outliers)),dtype="<u2").reshape(-1,2)
        assert len(flags)==masks+outliers and (flags[:,0]<ncols).all() and (flags[:,1]<nrows).all()
        flagged=flags[:,1].astype("int64")*ncols+flags[:,0]
        return records["mean"],flagged,ncols


def main(cache):
    manifest=[obtain(cache/"snp6.tar.gz",CHIP_URL)]
    db=cache/DB
    if not db.exists():
        db.parent.mkdir(parents=True,exist_ok=True)
        with tarfile.open(cache/"snp6.tar.gz") as tar, db.open("wb") as f:
            shutil.copyfileobj(tar.extractfile(DB),f)
    with sqlite3.connect(db) as conn:
        # Exclude centromeric q11.1, distal q12 and the terminal p11.31 band
        # to ensure that both sets are non-PAR irrespective of small build-boundary shifts.
        probes=pd.read_sql_query("""select s.fsetid,s.man_fsetid,s.chrom,s.chrom_start,s.cytoband,p.fid,p.x,p.y
           from featureSetCNV s join pmfeatureCNV p on s.fsetid=p.fsetid
           where s.chrom in ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22')
           or (s.chrom='Y' and s.xpar=0 and s.cytoband in ('p11.2','q11.21','q11.221','q11.222','q11.223','q11.23'))""",conn)
    assert probes.fid.is_unique and (probes.chrom_start>=0).all()
    probes["arm"]=np.where(probes.chrom!="Y","autosome",np.where(probes.cytoband=="p11.2","Yp","Yq"))
    arms=probes.set_index("fsetid")[["man_fsetid","arm","chrom_start","cytoband"]].drop_duplicates()
    assert arms.index.is_unique
    out=[]; yvalues=[]
    for cell,gsm in CELLS.items():
        path=cache/(gsm+".CEL.gz")
        manifest.append(obtain(path,f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{path.name}"))
        metaname=gsm+".brief.txt"
        manifest.append(obtain(cache/metaname,f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gsm}&targ=self&form=text&view=brief"))
        descriptor=(cache/metaname).read_text()
        assert f"^SAMPLE = {gsm}" in descriptor and "!Sample_molecule_ch1 = genomic DNA" in descriptor
        mean,flagged,cols=cel(path)
        index=probes.fid.to_numpy()-1
        assert (index==probes.y*cols+probes.x).all()
        vals=pd.Series(np.log2(1+mean[index].astype(float)),index=probes.fsetid)
        flag_mask=np.isin(index,flagged)
        for label,excluded in [("all",False),("excluding_flagged",True)]:
            locus=vals.mask(flag_mask) if excluded else vals
            locus=locus.groupby(level=0).median()
            vals_locus=arms.join(locus.rename("log_intensity"))
            center=vals_locus.loc[vals_locus.arm=="autosome","log_intensity"].median()
            vals_locus["centered"]=vals_locus.log_intensity-center
            if label=="all":
                yvalues.append(vals_locus.loc[vals_locus.arm!="autosome"].assign(gsm=gsm,cell=cell))
            for arm,g in vals_locus.groupby("arm"):
                out.append(dict(cell=cell,sample=gsm,probes=label,arm=arm,n_sets=g.centered.count(),autosome_center=center,median=g.centered.median(),
                                q25=g.centered.quantile(.25),q75=g.centered.quantile(.75)))
    pd.concat(yvalues).to_csv(OUT/"broad_y_probes.tsv.gz",sep="\t",index=False,float_format="%.7f",compression={"method":"gzip","mtime":0})
    summary=pd.DataFrame(out)
    summary.to_csv(OUT/"broad_summary.tsv",sep="\t",index=False,float_format="%.7f")
    pd.DataFrame(manifest).to_csv(OUT/"broad_sources.tsv",sep="\t",index=False)
    print(summary.query('arm!="autosome"')[["cell","probes","arm","n_sets","median"]].round(4).to_string(index=False))

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("cache",type=Path)
    main(p.parse_args().cache)
