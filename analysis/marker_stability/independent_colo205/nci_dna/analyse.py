"""Extract non-pseudoautosomal Y fluorescence from two deposited 450K studies."""
import argparse
import hashlib
from pathlib import Path
import re
from shutil import copyfileobj

from bs4 import BeautifulSoup
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests

OUT = Path(__file__).resolve().parent
SOURCES = {
    "GSE79185.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79185/suppl/GSE79185_NCI60_raw_int.txt.gz",
    "GSE86078.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE86nnn/GSE86078/suppl/GSE86078_non_normalized.txt.gz",
    "manifest.gz": "https://zhouserver.research.chop.edu/InfiniumAnnotation/HM450/HM450.hg38.manifest.tsv.gz",
    "bands.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
    "cellminer.html": "https://discover.nci.nih.gov/cellminer/celllineMetadata.do",
}
TARGET = "COLO205"
# Spelling variants used by the NCI intensity file.
ALIASES = {"IGR0V1": "IGROV1", "HL60": "HL60TB", "U031": "UO31", "A549": "A549ATCC"}


def norm(name):
    name = re.sub(r"[^A-Z0-9]", "", name.upper())
    return ALIASES.get(name, name)


def get_inputs(cache):
    cache.mkdir(exist_ok=True, parents=True)
    manifest = []
    for filename, url in SOURCES.items():
        dest = cache / filename
        if not dest.exists():
            part = dest.with_suffix(".part")
            with requests.get(url, stream=True, timeout=120) as response, part.open("wb") as f:
                response.raise_for_status()
                response.raw.decode_content = filename.endswith('.html')
                copyfileobj(response.raw, f)
            part.replace(dest)
        with dest.open("rb") as f:
            digest = hashlib.file_digest(f, "sha256").hexdigest()
        manifest.append(dict(file=filename, url=url, size=dest.stat().st_size, sha256=digest))
    pd.DataFrame(manifest).to_csv(OUT / "sources.tsv", sep="\t", index=False)


def amel_labels(cache):
    soup = BeautifulSoup((cache / "cellminer.html").read_bytes(), "html.parser")
    labels = {}
    for row in soup.select_one("table").select("tr"):
        link = row.find("a", attrs={"data-open-modal": True})
        if link is None:
            continue
        heading = row.find("th")
        for footnote in heading.select("sup"):
            footnote.decompose()
        cell = heading.get_text(" ", strip=True).split(":", 1)[1]
        modal = soup.find(id=link["aria-controls"])
        obs = modal.find(string=lambda x: x and x.strip() == "AMEL_1").find_parent("tr")
        key1, a, key2, b = [t.get_text(strip=True) for t in obs.select("td")]
        assert (key1, key2) == ("AMEL_1", "AMEL_2")
        if set((a,b)) <= {"X","Y"}:
            labels[norm(cell)] = "X/Y" if "Y" in (a,b) else "X-only"
    assert labels[TARGET] == "X-only" and labels["HCT15"] == "X/Y"
    return labels


def probe_manifest(cache):
    cols = ["probeID", "CpG_chrm", "CpG_beg", "CpG_end", "MASK_general", "MASK_mapping"]
    frame = pd.read_csv(cache/"manifest.gz", sep="\t", usecols=cols, low_memory=False)
    frame = frame.loc[frame.probeID.str.startswith("cg") & ~frame.MASK_general & ~frame.MASK_mapping].copy()
    assert frame.probeID.is_unique
    frame["arm"] = ""
    frame.loc[frame.CpG_chrm.str.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2])"), "arm"] = "auto"
    bands = pd.read_csv(cache/"bands.gz", sep="\t", header=None, names=["chrom","start","end","band","stain"])
    for b in bands.loc[(bands.chrom=="chrY") & (bands.stain!="acen")].itertuples():
        # GRCh38 PARs on chrY: 10,001–2,781,479 and 56,887,903–57,217,415 (1-based).
        keep = ((frame.CpG_chrm=="chrY") & (frame.CpG_beg>=2_781_479) &
                (frame.CpG_end<=56_887_902) & (frame.CpG_beg>=b.start) & (frame.CpG_end<=b.end))
        frame.loc[keep,"arm"] = "Y"+b.band[0]
    return frame.loc[frame.arm!=""].set_index("probeID")


def extract(study, cache, annotation, labels):
    path = cache/(study+".gz")
    nci = study == "GSE79185"
    header = 0 if nci else 4
    cols = list(pd.read_csv(path, sep="\t", skiprows=header, nrows=0).columns)
    suffix = ".Methylated" if nci else " Methylated Signal"
    other = ".UnMethylated" if nci else " Unmethylated Signal"
    cells = [k.removesuffix(suffix) for k in cols if k.endswith(suffix) and norm(k.removesuffix(suffix)) in labels]
    assert len(set(map(norm,cells))) == len(cells)
    meth = [k+suffix for k in cells]
    unmeth = [k+other for k in cells]
    autos, ys = [], []
    for piece in pd.read_csv(path, sep="\t", skiprows=header, usecols=[cols[0], *meth, *unmeth], chunksize=12_000, na_values=["null"]):
        piece = piece.set_index(cols[0])
        ann = annotation.reindex(piece.index)
        valid = ann.arm.notna()
        vals = np.log2(1.0 + piece.loc[valid, meth].to_numpy(dtype=np.float64) +
                       piece.loc[valid, unmeth].to_numpy(dtype=np.float64))
        groups = ann.loc[valid,"arm"]
        is_auto = groups.values=="auto"
        autos.append(vals[is_auto])
        ys.append(pd.DataFrame(vals[~is_auto],index=groups.index[~is_auto],columns=cells))
    auto = np.concatenate(autos)
    center = np.nanmedian(auto,axis=0)
    y = pd.concat(ys)
    assert y.index.is_unique and np.isfinite(center).all()
    y = y.subtract(pd.Series(center,index=cells),axis="columns")
    data = y.rename_axis("probeID").reset_index().melt(id_vars="probeID",var_name="line",value_name="log2_centered")
    data = data.join(annotation[["arm","CpG_beg"]],on="probeID")
    data["cell_key"] = data.line.map(norm)
    data["reference"] = data.cell_key.map(labels).where(data.cell_key!=TARGET, "COLO 205")
    data["study"] = study
    stat = data.groupby(["study","cell_key","reference","arm"]).log2_centered.agg(["count","median"]).reset_index()
    centers = pd.DataFrame({"study":study,"line":cells,"cell_key":list(map(norm,cells)),
                            "autosomal_median":center,"autosomal_probes":np.isfinite(auto).sum(axis=0)})
    return data, stat, centers


def write_plots(data, stat, stem="dna_intensity"):
    fig, axes = plt.subplots(2,2,figsize=(9.3,6.35),gridspec_kw={"width_ratios":[1,1.55]})
    colors={"X-only":"#77939b","X/Y":"#9368a0","COLO 205":"#c25326"}
    for index,(study,d) in enumerate(data.groupby("study")):
        ax,bx=axes[index]
        sub=stat.loc[stat.study==study]
        for i,arm in enumerate(["Yp","Yq"]):
            for offset,(g,color) in zip([-.2,.03,.25],colors.items()):
                vals=sub.loc[(sub.arm==arm) & (sub.reference==g)].sort_values("cell_key")
                jitter=np.random.default_rng(71).uniform(-.05,.05,len(vals)) if g!="COLO 205" else [0]
                ax.scatter(i+offset+np.array(jitter),vals["median"],s=13 if g!="COLO 205" else 53,
                           marker="o" if g!="COLO 205" else "D",color=color,alpha=.76,label=g if i==0 else None)
        counts = sub.loc[sub.arm=="Yp"].groupby("reference").size()
        ax.set(title=f"{study}  ·  {counts['X-only']} X-only, {counts['X/Y']} X/Y",xticks=[0,1],xticklabels=["Yp","Yq"],ylabel="Centered log₂ DNA intensity")
        for arm in ("Yp","Yq"):
            vals=d.loc[d.arm==arm].groupby(["probeID","CpG_beg","reference"]).log2_centered.median().unstack().reset_index().sort_values("CpG_beg")
            for g,color in colors.items():
                bx.scatter(vals.CpG_beg/1e6,vals[g],s=6 if g!="COLO 205" else 10,alpha=.45 if g!="COLO 205" else .85,c=color)
        bx.set(title="Chromosome Y probe positions",xlabel="GRCh38 position (Mb)",ylabel="Centered log₂ DNA intensity")
        bx.axvspan(10.3,10.6,color="#dddddd",zorder=-1)
        for a in (ax,bx):
            a.spines[["top","right"]].set_visible(False)
            a.tick_params(labelsize=8)
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="upper center",ncol=3,frameon=False)
    fig.tight_layout(rect=(0,0,1,.94),w_pad=1.7,h_pad=1)
    for ext in ["png","pdf"]:
        fig.savefig(OUT/(stem+"."+ext),dpi=185,bbox_inches="tight")
    plt.close(fig)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("cache",type=Path,help="directory for the public raw intensity files")
    cache=p.parse_args().cache
    get_inputs(cache)
    ann=probe_manifest(cache)
    labels=amel_labels(cache)
    runs=[extract(s,cache,ann,labels) for s in ("GSE79185","GSE86078")]
    data,stats,centers=(pd.concat(parts,ignore_index=True) for parts in zip(*runs))
    assert set(data.arm)=={"Yp","Yq"} and set(stats.cell_key[stats.reference=="COLO 205"])=={TARGET}
    data.to_csv(OUT/"y_probe_intensities.tsv.gz",sep="\t",index=False,float_format="%.7f",compression={"method":"gzip","mtime":0})
    for obj,stem in [(stats,"arm_medians"),(centers,"autosomal_reference")]:
        obj.to_csv(OUT/(stem+".tsv"),sep="\t",index=False,float_format="%.7f")
    summary=stats.groupby(["study","arm","reference"])["median"].agg(["size","median","min","max"])
    summary.to_csv(OUT/"arm_summary.tsv",sep="\t",float_format="%.7f")
    write_plots(data,stats)
    print(summary.to_string(float_format=lambda x:f"{x:.3f}"))
    print(data.drop_duplicates("probeID").groupby("arm").CpG_beg.agg(["count","min","max"]))

if __name__ == "__main__":
    main()
