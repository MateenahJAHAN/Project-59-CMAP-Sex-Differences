from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

P=Path(__file__).resolve().parent
w=pd.read_csv(P/"well_counts.tsv.gz",sep="\t")
w=w[(w.cell_line_id=="CVCL_0218") & (w.cells>=200)].copy()
genes=["RPS4Y1","DDX3Y","KDM5D","UTY","EIF1AY","ZFY","USP9Y","NLGN4Y"]
w["score"]=np.log2(1+w[genes].multiply(1e6/w.total_umi, axis=0)).mean(axis=1)
cut=pd.read_csv(P/"reference_cutoffs.tsv",sep="\t").set_index("score").cutoff.y_score
line=pd.read_csv(P/"by_line.tsv",sep="\t").query("treated").sort_values("joint_agreement")
rep=pd.read_csv(P/"replicate_pairs.tsv.gz",sep="\t").query("cell_line_id == 'CVCL_0218'")
scores=w.merge(pd.read_csv(P/"sample_key.tsv",sep="\t")[["sample","dose_um"]],on="sample")
wide=scores[(scores.drug!="DMSO_TF") & (scores.plate.isin(["plate6","plate14"]))].pivot(index=["drug","dose_um"],columns="plate",values="score").dropna()
assert len(wide)==len(rep)
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False})
fig,axes=plt.subplots(1,2,figsize=(9.5,3.35),gridspec_kw={"width_ratios":[1,1.18]},constrained_layout=True)
ax=axes[0]
ax.scatter(np.arange(2,len(line)+1),100*line.joint_agreement.iloc[1:],s=12,color="#376f84",alpha=.78)
ax.scatter([1],100*line.joint_agreement.iloc[:1],s=28,color="#b8562c")
ax.annotate("COLO 205",(1,100*line.joint_agreement.iloc[0]),xytext=(5,49),color="#9d461f",arrowprops={"arrowstyle":"-","color":"#aa572e"})
ax.set(xlabel="Cell lines ordered by agreement",ylabel="Treated wells matching baseline (%)",ylim=(32,104),xticks=[1,12,24,36,48]);ax.set_title("A   Agreement across the panel",loc="left",fontsize=11)
ax=axes[1]
color=np.where((wide.plate6>=cut)&(wide.plate14>=cut),"#aa512b","#38738a")
ax.scatter(wide.plate6,wide.plate14,s=16,c=color,alpha=.8,edgecolors="white",linewidths=.25)
extent=[min(wide.min())-.12,max(wide.max())+.12]
ax.plot(extent,extent,":",lw=.8,color="#9ca6ad");ax.axvline(cut,color="#6b747b",lw=.9);ax.axhline(cut,color="#6b747b",lw=.9)
ax.set(xlabel="Y score on plate 6",ylabel="Y score on plate 14",xlim=extent,ylim=extent)
ax.set_title("B   COLO 205 across repeat plates",loc="left",fontsize=11)
fig.savefig(P/"marker_agreement.png",dpi=220,bbox_inches="tight");plt.close(fig)
diag=w.groupby(w.drug.ne("DMSO_TF")).score.describe().to_dict(orient="index")
print("COLO 205 score distributions",diag,"cut",cut)
