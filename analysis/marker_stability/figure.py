"""Render the observed line-level concentration contrasts and cross-plate repeats."""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
plt.rcParams.update({"font.size":9.0, "axes.spines.top":False, "axes.spines.right":False,
                     "pdf.fonttype":42, "savefig.facecolor":"white"})
lines = pd.read_csv(RESULTS/"lines.tsv",sep="\t").query("in_primary")
rep = pd.read_csv(RESULTS/"repeated_conditions.tsv",sep="\t")
fig,(a,b) = plt.subplots(1,2,figsize=(8.4,3.65),layout="constrained",width_ratios=(1.05,1))
colors=["#a34a1b","#14636a"]
a.axhline(0,color="#999999",linewidth=.8)
for x, organ in enumerate(["Bowel","Lung","Skin"]):
    for grp,offset,clr,label in zip([False,True],[-.19,.19],colors,["Neither marker","Y genes or XIST"]):
        obs=lines.loc[(lines.organ==organ)&(lines.positive==grp),"mean_high_minus_low_log2_share"]
        a.scatter(x+offset+np.linspace(-.07,.07,len(obs)),obs,s=23,color=clr,label=label if x==0 else None,zorder=3)
        a.plot([x+offset-.1,x+offset+.1],[obs.mean()]*2,color=clr,linewidth=2.5)
a.set_xticks(range(3),["Bowel","Lung","Skin"])
a.set_ylabel("Mean high − low concentration\n(log₂ relative cell recovery)")
a.set_xlim(-.6,2.6)
a.set_title("A   Responses by cell line",loc="left",fontweight="bold")
a.legend(loc="upper right",fontsize=8,frameon=False)
for flag,color,label in [(False,"#9c6b42","Other plate sets"),(True,"#236f76","Plates 6 and 14 only")]:
    pts=rep.loc[rep.only_6_14==flag]
    b.scatter(pts["first"],pts.others,s=14,alpha=.66,color=color,label=label,edgecolors="none")
limits=[min(rep["first"].min(),rep.others.min())-.025,max(rep["first"].max(),rep.others.max())+.025]
b.plot(limits,limits,linewidth=.9,ls="--",color="0.5")
b.set_xlim(limits)
b.set_ylim(limits)
b.set_aspect("equal",adjustable="box")
b.set_title("B   Repeated drug–dose conditions",loc="left",fontweight="bold")
b.set_xlabel("Marker-group contrast on first plate (log₂)")
b.set_ylabel("Contrast on other plates (log₂)")
b.legend(loc="upper left",fontsize=7.5,frameon=False)
for suffix in ["png","pdf"]:
    fig.savefig(HERE/f"figure_1.{suffix}",dpi=240,metadata={"Creator":"Project 59 follow-up"})
