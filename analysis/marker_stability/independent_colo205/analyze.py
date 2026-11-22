"""Rebuild the COLO 205 analysis from the inputs listed in sources.tsv."""
from pathlib import Path
from zipfile import ZipFile
import csv
import gzip
import re

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
IN = ROOT / "inputs"
OUT = ROOT / "results"
OUT.mkdir(exist_ok=True)
MARKERS = {
    "Yp": ["RPS4Y1", "ZFY"],
    "Yq": ["DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "NLGN4Y"],
    "X": ["XIST"],
    "autosome": ["ACTB", "GAPDH", "RPLP0"],
}
GENES = {g: region for region, genes in MARKERS.items() for g in genes}
LINES = ["COLO205", "HCT116", "HT29", "HCT15"]
CELL_LABELS = dict(zip(["CO:COLO 205", "CO:HCT-116", "CO:HT29", "CO:HCT-15"], LINES))


def write(rows, name):
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / (name + ".tsv"), sep="\t", index=False, float_format="%.8g", na_rep="NA")
    return frame


def rna_panel():
    # UCSC's comprehensive GENCODE v33 genePred includes transcript versions. Limit
    # gene assignment to primary chromosomes and reject ambiguous transcript IDs.
    lookup = {}
    with gzip.open(IN / "wgEncodeGencodeCompV33.txt.gz", "rt") as f:
        for line in f:
            rec = line.rstrip().split("\t")
            gene, chromosome, tx = rec[12], rec[2], rec[1]
            if gene in GENES and re.fullmatch(r"chr([0-9]+|X|Y)", chromosome):
                assert tx not in lookup or lookup[tx] == gene
                lookup[tx] = gene
    rna = []
    for acc, cell in [("GSM5857280", "COLO205"), ("GSM5857281", "HCT116")]:
        dat = pd.read_csv(next(IN.glob(f"{acc}*tsv.gz")), sep="\t")
        dat["gene"] = dat["target_id"].map(lookup)
        assert dat.target_id.is_unique
        for gene in GENES:
            g = dat[dat.gene == gene]
            tpm = g.tpm.sum(min_count=1)
            rna.append(dict(study="IRIC", cell=cell, gene=gene, arm=GENES[gene],
                            source=acc, n_transcripts=len(g), est_counts=g.est_counts.sum(min_count=1),
                            value=tpm, units="TPM", log2_one_plus=np.log2(1+tpm)))
    with ZipFile(IN / "nci60_RNA__RNA_seq_composite_expression.zip") as z:
        with z.open(next(x for x in z.namelist() if x.endswith(".html"))) as raw:
            for line in raw:
                if b"<th" in line and b"Gene name" in line:
                    header = re.findall(r"<th[^>]*>(.*?)</th>", line.decode())
                elif b'<td scope="row">' in line:
                    if not any((f">{g}</td>").encode() in line[:90] for g in GENES):
                        continue
                    vals = re.findall(r"<td[^>]*>(.*?)</td>", line.decode())
                    gene = vals[0]
                    if gene not in GENES:
                        continue
                    for cell_header, cell in CELL_LABELS.items():
                        v = float(vals[header.index(cell_header)])
                        rna.append(dict(study="NCI", cell=cell, gene=gene, arm=GENES[gene],
                                        source="CellMiner v2.16", n_transcripts=None, est_counts=None,
                                        value=v, units="log2(FPKM+1)", log2_one_plus=v))
    df = write(rna, "rna_genes")
    assert len(df) == 6 * len(GENES) and not df.duplicated(["study", "cell", "gene"]).any()
    a = df[df.arm.isin(["Yp","Yq"])].pivot_table(index=["study","cell"], columns="arm", values="log2_one_plus", aggfunc="median")
    a["Yq_minus_Yp"] = a.Yq - a.Yp
    write(a.reset_index(), "rna_arms")


def agilent_groups(chrom, start, end):
    # Inclusive hg17 PARs from the UCSC genome assembly page, and centromeric acen bands.
    if chrom.isdigit():
        return "autosome"
    if chrom == "X":
        return "X"
    if chrom == "Y":
        if start > 2692881 and end <= 11200000:
            return "Yp"
        if start > 12400000 and end < 57372174:
            return "Yq"
        return "Y_excluded"
    return "other"


def raw_array(acc):
    with gzip.open(next(IN.glob(f"{acc}*44k*")), "rt") as f:
        for line in f:
            if line.startswith("FEATURES"):
                break
        use = ["ControlType", "ProbeName", "SystematicName", "GeneName",
               *[c+t for c in "rg" for t in ("ProcessedSignal","IsWellAboveBG","IsSaturated","IsFeatNonUnifOL")]]
        d = pd.read_csv(f, sep="\t", names=line.rstrip().split("\t"), usecols=use, low_memory=False)
    d = d[d.ControlType == 0].copy()
    m = d.SystematicName.str.extract(r"^chr(\d+|X|Y):(\d+)-(\d+)$")
    d[["chrom","start","end"]] = m
    d = d.dropna(subset=["chrom"]).copy()
    d[["start","end"]] = d[["start","end"]].astype(int)
    d["group"] = [agilent_groups(c,s,e) for c,s,e in zip(d.chrom,d.start,d.end)]
    # If a probe name maps to conflicting positions, drop it entirely.
    nmap = d.groupby("ProbeName").SystematicName.nunique()
    return d[d.ProbeName.isin(nmap[nmap == 1].index)].copy()


def agilent():
    first = raw_array("GSM1181248")
    features = first.drop_duplicates("ProbeName").set_index("ProbeName")
    p = IN / "GSE48568_series_matrix.txt.gz"
    with gzip.open(p, "rt") as f:
        for line in f:
            if line.startswith("!Sample_title"):
                names=next(csv.reader([line],delimiter="\t"))[1:]
            if line.startswith("!series_matrix_table_begin"):
                break
        m=pd.read_csv(f,sep="\t",index_col=0,comment="!",na_values=["null"])
    assert m.index.is_unique
    assert m.index.isin(features.index).mean() > 0.95
    m = m.loc[m.index.intersection(features.index)]
    meta = dict(zip(m.columns, names))
    grouped = m.groupby(features.group).median()
    result=[]
    for acc,title in meta.items():
        background=grouped.loc["autosome",acc]
        for region in ["autosome","X","Yp","Yq","Y_excluded"]:
            values=m.loc[features.reindex(m.index).group==region,acc].dropna()
            result.append(dict(accession=acc, sample=title, region=region, n_probes=len(values),
                               median_log2_ratio=values.median(), centered=values.median()-background))
    tab=write(result,"agilent_all_arrays")
    centers=tab.pivot(index=["accession","sample"],columns="region",values="centered").reset_index()
    centers["Yq_minus_Yp"]=centers.Yq-centers.Yp
    write(centers, "agilent_arm_contrasts")
    sel={acc:title for acc,title in meta.items() if re.search(r"CO:(COLO205|HCT_116|HT29|HCT_15)_",title)}
    fY=features[features.chrom=="Y"][["chrom","start","end","group","GeneName"]]
    write(fY.join(m[list(sel)],how="inner").reset_index(names="probe"),"agilent_Y_probes")
    raw=[]
    for acc,title in sel.items():
        file = IN / (acc+".metadata.txt")
        if not file.exists():
            continue
        d=raw_array(acc)
        # These GEO samples always put Cy5 in channel 1; test DNA can be channel 1 or 2.
        info=file.read_text()
        assert "!Sample_label_ch1 = Cy5" in info and "!Sample_label_ch2 = Cy3" in info
        color="g" if "!Sample_source_name_ch1 = Reference" in info else "r"
        val=color+"ProcessedSignal"
        d["detected"]=d[color+"IsWellAboveBG"].astype(bool)
        d["acceptable"]=(d[color+"IsSaturated"]==0)&(d[color+"IsFeatNonUnifOL"]==0)
        d=d[d.acceptable & (d[val]>0)]
        # Collapse repeated spots, so every probe counts once.
        d=d.groupby(["ProbeName","group"]).agg(signal=(val,"median"), detected=("detected","mean")).reset_index()
        auto=np.log2(d.loc[d.group=="autosome","signal"]).median()
        for region in ["autosome","X","Yp","Yq"]:
            rows=d[d.group==region]
            raw.append(dict(accession=acc,sample=title,channel=color,region=region,n_probes=len(rows),
                            fraction_above_background=rows.detected.mean(), median_log2_signal=np.log2(rows.signal).median(),
                            centered=np.log2(rows.signal).median()-auto))
    write(raw,"agilent_raw_sensitivity")


def snp_array():
    # The deposited rows omit the platform's sex-chromosome probes. Count each
    # join directly in BOTH archives. Zero is coverage, never a dosage estimate.
    with (IN/"GPL8887.platform.txt").open() as f:
        for l in f:
            if l.startswith("!platform_table_begin"):
                break
        next(f)
        annotation={}
        for l in f:
            if l.startswith("!platform_table_end"):
                break
            address, chromosome, *_ = l.rstrip().split("\t")
            assert address not in annotation
            annotation[address] = chromosome
    mapper=pd.Series(annotation)
    rows=[]
    for name, ids in [("GPL8887", mapper.index), *[(file.name, pd.read_csv(file,sep="\t",usecols=[0],dtype=str).iloc[:,0]) for file in sorted(IN.glob("GSE55832_*.gz"))]]:
        chrom=mapper.reindex(list(ids))
        assert not pd.Index(ids).has_duplicates
        for group in ["autosome","X","Y","XY","Mt","unmapped"]:
            if group=="autosome":
                n=chrom.fillna("").str.isnumeric().sum()
            elif group=="unmapped":
                n=chrom.isna().sum()
            else:
                n=(chrom==group).sum()
            rows.append(dict(source=name,region=group,n=int(n),total=len(chrom)))
    write(rows,"wehi_coverage")

if __name__=="__main__":
    rna_panel()
    agilent()
    snp_array()
