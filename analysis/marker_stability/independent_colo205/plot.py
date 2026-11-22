"""Plot values from the analyzed public data; dotted intervals mark centromeres."""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
mpl.rcParams["svg.hashsalt"]="colo205-external-arrays"
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from analyze import GENES, OUT
ROOT=Path(__file__).resolve().parent
COLOR={"COLO205":"#145f82","HCT15":"#8c5c1a","HCT116":"#677384","HT29":"#8e6880"}
LABEL={"COLO205":"COLO 205","HCT15":"HCT-15","HCT116":"HCT 116","HT29":"HT-29"}

def draw():
    rna=pd.read_csv(OUT/"rna_genes.tsv",sep="\t")
    ag=pd.read_csv(OUT/"agilent_arm_contrasts.tsv",sep="\t")
    gd=pd.read_csv(OUT/"gdsc_Y_probes.tsv",sep="\t")
    fig=plt.figure(figsize=(10.6,8.3),layout="constrained")
    gs=fig.add_gridspec(2,2,height_ratios=[1,1.05],width_ratios=[1.3,1])
    axs=[fig.add_subplot(gs[0,0]),fig.add_subplot(gs[0,1]),fig.add_subplot(gs[1,:])]
    fig.suptitle("COLO 205 shows stronger Yq than Yp signal across RNA and DNA assays",fontsize=15, fontweight="bold")
    ax=axs[0]
    target=rna[(rna.cell=="COLO205") & rna.arm.isin(["Yp","Yq"])]
    genes=list(GENES)[:8]
    for st,lab,shape,off,col in [("IRIC","Montréal · log₂(1+TPM)","o",-0.1,"#145f82"),("NCI","NCI · log₂(1+FPKM)","s",0.1,"#bf7641")]:
        sub=target[target.study==st].set_index("gene").loc[genes]
        ax.scatter(np.arange(len(sub))+off,sub.log2_one_plus,s=40,c=col,marker=shape,label=lab,zorder=3)
    ax.axvline(1.5,color="#aab3b8",ls=":",lw=1)
    ax.set_xticks(range(8),genes,rotation=44,ha="right")
    ax.set_ylabel("Reported RNA abundance, log scale")
    ax.set_title("A  |  Two bulk RNA studies",loc="left")
    ax.annotate("Yp",(0.45,2.03),ha="center",color="#57636c");ax.annotate("Yq",(4.5,2.03),ha="center",color="#57636c")
    ax.set_ylim(-.12,2.42);ax.legend(loc="upper left",bbox_to_anchor=(0,-.29),ncol=2,fontsize=8,frameon=False)
    ax=axs[1]
    a=ag.copy();a["line"]=a["sample"].str.replace(r"^.*:|_rep\d+$","",regex=True).str.replace("_","")
    a=a[a.line.isin(COLOR)]
    for i,line in enumerate(COLOR):
        vals=sorted(a.loc[a.line==line,"Yq_minus_Yp"]);assert len(vals)==2
        ax.plot(vals,[i,i],color=COLOR[line],lw=2)
        ax.scatter(vals,[i,i],c=COLOR[line],s=38,edgecolor="white",linewidth=.5,zorder=3)
    ax.set_yticks(range(len(COLOR)),[LABEL[k] for k in COLOR]);ax.invert_yaxis()
    ax.axvline(0,ls=":",lw=1,color="#889098");ax.set_xlim(-.19,1.04)
    ax.set_xlabel("Median Yq − median Yp, log₂ ratio");ax.set_title("B  |  NCI genomic array CGH",loc="left")
    ax.text(0,-.27,"Two dye swaps per line · 45 Yp and 90 Yq probes",transform=ax.transAxes,fontsize=8,color="#4f5962")
    ax=axs[2]
    for line in ["HCT15","HT29","COLO205"]:
        d=gd.loc[gd.cell==line].copy()
        if line=="COLO205":
            ax.scatter(d.CpG_beg/1e6,d.centered,s=9,alpha=.45,color=COLOR[line],lw=0)
        for arm in ["Yp","Yq"]:
            t=d[d.arm==arm]
            ax.hlines(t.centered.median(),t.CpG_beg.min()/1e6,t.CpG_beg.max()/1e6,colors=COLOR[line],ls="-" if line=="COLO205" else "--",lw=2.1)
    ax.set_title("C  |  GDSC genomic DNA signal at individual probes",loc="left")
    ax.set_xlabel("Position on chromosome Y, GRCh38 (Mb)");ax.set_ylabel("Log₂ DNA intensity relative\nto autosomal median")
    ax.axhline(0,color="#aab3b8",lw=.7);ax.axvspan(10.3,10.6,color="#d6dce0",alpha=.8)
    ax.text(.96,.045,"81 Yp probes · 161 Yq probes",transform=ax.transAxes,ha="right",fontsize=9)
    ax.legend(handles=[Line2D([0],[0],color=COLOR[x],lw=2,ls="-" if x=="COLO205" else "--",label=LABEL[x]+" arm medians") for x in ["COLO205","HT29","HCT15"]],ncol=3,fontsize=8,loc="lower left",bbox_to_anchor=(0,-.32),frameon=False)
    for ax in axs:
        ax.spines[["top","right"]].set_visible(False)
        ax.grid(axis="y",alpha=.15,lw=.6)
    fig.savefig(ROOT/"colo205_rna_dna.png",dpi=210,bbox_inches="tight")
    fig.savefig(ROOT/"colo205_rna_dna.svg",bbox_inches="tight",metadata={"Date":None})
if __name__=="__main__":
    draw()
