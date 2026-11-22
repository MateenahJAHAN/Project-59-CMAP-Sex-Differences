"""Plot arm differences on each array's own scale and the Broad positional sensitivity."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

P=Path(__file__).resolve().parent
a=pd.read_csv(P/"complete_arm_medians.tsv",sep="\t").rename(columns={"cell_key":"cell"})
b=pd.read_csv(P/"broad_summary.tsv",sep="\t").query('probes=="all" and arm!="autosome"').copy()
b["study"]="CCLE"
rows=pd.concat([a[["study","cell","arm","median"]],b[["study","cell","arm","median"]]])
groups=rows.pivot(index=["study","cell"],columns="arm",values="median").reset_index()
groups["difference"]=groups.Yq-groups.Yp
colors={"COLO205":"#ae451e","HCT15":"#3d7f80","HT29":"#56557e"}
fig=plt.figure(figsize=(9.2,4.8),layout="constrained")
left,right=fig.subplots(1,2,width_ratios=[1,1.55])
for i,study in enumerate(["GSE79185","GSE86078","CCLE"]):
    for r in groups.loc[groups.study==study].itertuples():
        name=r.cell
        jitter=(np.random.default_rng(sum(name.encode())).uniform(-.09,.09) if study!="CCLE" else {"COLO205":0,"HCT15":-.09,"HT29":.09}[name])
        left.scatter(i+jitter,r.difference,color=colors.get(name,"#a4afb0"),marker="D" if name=="COLO205" else "o",
                     s=50 if name=="COLO205" else 15,zorder=4 if name in colors else 1)
left.set(title="A   Within-array arm difference",xticks=range(3),xticklabels=["NCI\n450K","Candiolo\n450K","Broad\nSNP 6.0"],
         ylabel="Median log₂ intensity, Yq minus Yp",xlim=(-.45,2.45))
left.axhline(0,c="#777",lw=.6,zorder=0)
points=pd.read_csv(P/"broad_y_probes.tsv.gz",sep="\t")
# Equal-width 0.5 Mb bins. Locations are in the platform annotation's native coordinate system.
points["bin"]=(points.chrom_start//500_000)*.5+.25
for j,(cell,label) in enumerate([("COLO205","COLO 205"),("HCT15","HCT-15"),("HT29","HT-29")]):
    for arm,sub in points.loc[points.cell==cell].groupby("arm"):
        summary=sub.groupby("bin").centered.agg(["median","count"])
        assert summary["count"].sum()==len(sub)
        right.scatter(summary.index,summary["median"],label=label if arm=="Yp" else None,c=colors[cell],marker=["D","o","s"][j],s=13,
                      zorder=3 if cell=="COLO205" else 2)
right.set(title="B   Broad genomic DNA by position",xlabel="Chromosome Y position in array annotation (Mb)",
          ylabel="Centered log₂ intensity per 0.5 Mb bin")
right.legend(loc="lower center",bbox_to_anchor=(.54,-.32),ncol=3,frameon=False,handletextpad=.2,columnspacing=.9,fontsize=8)
for ax in (left,right):
    ax.spines[["top","right"]].set_visible(False)
    ax.tick_params(labelsize=9)
    ax.title.set_fontsize(12)
    ax.title.set_ha("left")
    ax.title.set_position((0,1))
for ext in ("png","pdf"):
    fig.savefig(P/("dna_platforms."+ext),dpi=220,bbox_inches="tight")
