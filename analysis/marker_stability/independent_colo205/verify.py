"""Independent checks using the original GTF, individual probe counts and source XML."""
from pathlib import Path
from statistics import median
from collections import defaultdict
import csv
import gzip
import json
import re
import xml.etree.ElementTree as ET
import pandas as pd
from analyze import IN, OUT, GENES, write

def verify():
    # Cross-check transcript sums using GENCODE's own GTF, separate from UCSC genePred.
    tx={}
    with gzip.open(IN/"gencode.v33.annotation.gtf.gz","rt") as f:
        for row in f:
            if row.startswith("#") or "\ttranscript\t" not in row:
                continue
            parts=row.split("\t")
            if not re.fullmatch(r"chr(?:[1-9]|1\d|2[012]|X|Y)",parts[0]):
                continue
            attr=parts[8]; gene=re.search(r'gene_name "([^"]+)"',attr).group(1)
            if gene in GENES:
                e=re.search(r'transcript_id "([^"]+)"',attr).group(1)
                tx[e]=gene
    expected=pd.read_csv(OUT/"rna_genes.tsv",sep="\t")
    for gsm,cell in [("GSM5857280","COLO205"),("GSM5857281","HCT116")]:
        measured=defaultdict(list)
        with gzip.open(next(IN.glob(gsm+"*tsv.gz")),"rt") as f:
            for r in csv.DictReader(f,delimiter="\t"):
                g=tx.get(r["target_id"])
                if g:measured[g].append(float(r["tpm"]))
        for r in expected[(expected.study=="IRIC") & (expected.cell==cell)].itertuples():
            assert len(measured[r.gene])==int(r.n_transcripts)
            assert abs(sum(measured[r.gene])-r.value)<max(1e-6,abs(r.value)*1e-6) if measured[r.gene] else pd.isna(r.value)
    # Compare ratio calculations with a stand-alone per-probe median from the saved table.
    a=pd.read_csv(OUT/"agilent_Y_probes.tsv",sep="\t")
    b=pd.read_csv(OUT/"agilent_arm_contrasts.tsv",sep="\t")
    ranks=[]
    for acc in ["GSM1181248","GSM1181264"]:
        contrast=median(a.loc[a.group=="Yq",acc].dropna())-median(a.loc[a.group=="Yp",acc].dropna())
        assert abs(contrast-b.loc[b.accession==acc,"Yq_minus_Yp"].item())<1e-7
        ranks.append(1+int((b.Yq_minus_Yp>contrast+1e-6).sum()))
    # Compare the same oligo across lines, pairing the orientation of the test dye.
    contrasts=[]
    for target,refs in [("GSM1181248",["GSM1181250","GSM1181251","GSM1181252"]),
                        ("GSM1181264",["GSM1181367","GSM1181363","GSM1181366"])]:
        for ref,label in zip(refs,["HCT116","HCT15","HT29"]):
            d=a[target]-a[ref]
            arms={arm:d[a.group==arm].median() for arm in ["Yp","Yq"]}
            contrasts.append(dict(target=target,control=ref,control_cell=label,
                                  matched_Yp=int(d[a.group=="Yp"].count()),matched_Yq=int(d[a.group=="Yq"].count()),
                                  median_Yp=arms["Yp"],median_Yq=arms["Yq"],Yq_minus_Yp=arms["Yq"]-arms["Yp"]))
    write(contrasts,"agilent_same_probe_line_contrasts")
    # Repeat 450K summaries separately in each chemistry as a cross-hybridization/assay check.
    probes=pd.read_csv(OUT/"gdsc_Y_probes.tsv",sep="\t")
    design=pd.read_csv(IN/"manifest.gz",usecols=["probeID","designType"],sep="\t")
    by=probes.merge(design,on="probeID",validate="many_to_one").groupby(["cell","designType","arm"]).centered.agg(["median","size"]).reset_index()
    write(by,"gdsc_by_probe_design")
    sky=ET.parse(IN/"nci_SKY.xml")
    case=next(c for c in sky.iter("SkyCase") if c.findtext("SkyCase_caseId")=="2801")
    assert case.findtext(".//CaseHeader_name")=="Colo205 (submitter:Tonon)"
    iscn=case.findtext(".//SkyCellHeader_karyotype")
    assert iscn.startswith("69~75,XX,-Y[9],")
    result={"sky_case":2801,"pubmed":case.findtext(".//SkyCitation_pubMedId"),
            "karyotype":iscn,"agilent_target_ranks_out_of_129_in_rep_order":ranks,
            "IRIC_independent_GENCODE_GTF_all_panel_values_match":True}
    (OUT/"verification.json").write_text(json.dumps(result,indent=2)+"\n")
    (OUT/"sky_case_2801.xml").write_bytes(ET.tostring(case,encoding="utf-8"))
if __name__=="__main__":
    verify()
