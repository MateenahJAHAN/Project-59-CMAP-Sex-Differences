"""Figures for the original Phases 3–4 screens using the complete replay."""
from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
ROOT=Path(os.environ.get('P59_PHASE_RUN', Path(__file__).resolve().parent))
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'axes.spines.left':False,'pdf.fonttype':42})
colors={'A':'#08775F','B':'#A44F13','C':'#586367'}
def save(fig,stem):
    for ext in ['png','pdf']:
        fig.savefig(FIG/(stem+'.'+ext),dpi=180,bbox_inches='tight',facecolor='white')
    plt.close(fig)
def base(n):
    fig,ax=plt.subplots(figsize=(8.3,12.3));fig.subplots_adjust(left=.30,right=.97,bottom=.07,top=.88)
    ax.grid(axis='x',color='#DADADA',linewidth=.6);ax.set_axisbelow(True);ax.tick_params(axis='y',length=0,labelsize=8.5)
    return fig,ax
def phase3():
    t=pd.read_csv(ROOT/'comparison/phase3_compounds.csv').copy()
    t['oldhits']=t.original_A+t.original_B;t['newhits']=t.replay_A+t.replay_B
    t=t.sort_values(['newhits','drug'],ascending=[True,False]).reset_index(drop=True)
    fig,ax=base(len(t))
    for y,r in t.iterrows():
        ax.plot([r.oldhits,r.newhits],[y,y],c='#899393',linewidth=1,zorder=1)
        ax.scatter(r.oldhits,y,marker='|',s=95,color='#111111',linewidth=1,zorder=2)
        ax.scatter(r.newhits,y,s=23,facecolor=colors[r.replay] if r.tested_lineages else 'white',edgecolor=colors[r.replay],zorder=3)
    ax.set_yticks(np.arange(len(t)),t.drug);ax.set_xlim(-.8,max(t.newhits.max(),t.oldhits.max())*1.08+.5)
    ax.set_xlabel('Significant drug–lineage–tissue comparisons per compound (FDR < 0.05)')
    fig.text(.03,.97,'Phase 3   Enrichment within cell lineages',fontsize=17,fontweight='bold')
    fig.text(.03,.942,'Dots show replay counts. Ticks show the archived counts.',fontsize=9)
    handles=[Line2D([],[],marker='o',linestyle='',color=v,label=k) for k,v in [('Own tissue',colors['A']),('Other tissue only',colors['B']),('Tested, no enrichment',colors['C'])]]
    handles.append(Line2D([], [], marker='o', linestyle='', markerfacecolor='white', markeredgecolor=colors['C'], label='No eligible lineage'))
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.03,.932),frameon=False,ncol=4)
    save(fig,'phase3_within_lineage_counts')
def phase4():
    all_t=pd.concat([pd.read_csv(f) for f in sorted((ROOT/'outputs/phase4').glob('*.csv.gz'))],ignore_index=True)
    all_t['significant']=all_t.FDR<.05
    g=all_t.groupby('drug').agg(sets_tested=('FDR','size'),sets_significant=('significant','sum')).reset_index()
    a=pd.read_csv(ROOT/'comparison/phase4_pair_counts.csv').groupby('drug')[['A','B','C','replay_A','replay_B','replay_C']].sum().reset_index()
    g=g.merge(a,on='drug',validate='one_to_one')
    g['replay_percent']=100*g.sets_significant/g.sets_tested
    g['archive_count_over_replay_denominator']=100*(g.A+g.B)/g.sets_tested
    assert (g.sets_significant==g.replay_A+g.replay_B).all()
    g.sort_values(['replay_percent','drug'],ascending=[True,False],inplace=True,ignore_index=True)
    g.to_csv(ROOT/'comparison/phase4_figure_denominators.csv',index=False)
    fig,ax=base(len(g))
    for y,r in g.iterrows():
        ax.plot([r.archive_count_over_replay_denominator,r.replay_percent],[y,y],c='#71807F',linewidth=1)
        ax.scatter(r.archive_count_over_replay_denominator,y,marker='|',s=85,color='black',linewidth=1,zorder=2)
        ax.scatter(r.replay_percent,y,s=23,color='#08775F',zorder=3)
    ax.set_yticks(np.arange(len(g)),g.drug);ax.set_xlim(-.5,max(g.replay_percent.max(),g.archive_count_over_replay_denominator.max())*1.07+1)
    ax.set_xlabel('Significant split reference sets / tested sets (%)')
    fig.text(.03,.97,'Phase 4   How often did a split gene set reach the threshold?',fontsize=16,fontweight='bold')
    fig.text(.03,.94,'Every tested combination of drug, lineage and reference set contributes to the denominator.\nDots use replayed calls. Ticks use archived calls divided by the same denominator (FDR < 0.05).',fontsize=9)
    save(fig,'phase4_tested_reference_sets')
def phase5():
    frame=pd.concat([pd.read_csv(f) for f in sorted((ROOT/'outputs/phase5').glob('*.csv.gz'))],ignore_index=True)
    counts=frame.groupby(['drug','lineage','category']).size().unstack(fill_value=0).reset_index()
    active=counts[counts.D<44].copy()
    active=active.sort_values(['A','B','C','drug','lineage'],ignore_index=True)
    fig,ax=plt.subplots(figsize=(9,6.6));fig.subplots_adjust(left=.39,right=.97,bottom=.10,top=.8)
    y=np.arange(len(active));left=np.zeros(len(active))
    for field,color,label in [('A','#08775F','Own-tissue enrichment'),('B','#BA6D2B','Other-tissue enrichment'),('C','#9DB1BB','Tested, FDR ≥ 0.05'),('D','#E1E4E5','Untested tissue')]:
        ax.barh(y,active[field],left=left,color=color,label=label,height=.70);left+=active[field].to_numpy()
    assert np.all(left==44)
    labs=[r.drug+' · '+r.lineage.replace('haematopoietic_and_lymphoid_tissue','blood / lymphoid').replace('_',' ') for r in active.itertuples()]
    ax.set_yticks(y,labs);ax.tick_params(axis='y',length=0);ax.set_xlim(0,44);ax.set_xticks([0,10,20,30,40,44]);ax.set_xlabel('GTEx reference tissues per drug–lineage pair')
    fig.text(.02,.963,'Phase 5   Tissue enrichment after pooling drug-response signatures',fontsize=15,fontweight='bold')
    fig.text(.02,.921,f'{len(active)} of {len(counts)} drug–lineage pairs yielded a test. {int((frame.category!="D").sum()):,} of {len(frame):,} tissue comparisons were tested.')
    handles,labels=ax.get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper left',bbox_to_anchor=(.01,.898),ncol=2,frameon=False,fontsize=9)
    save(fig,'phase5_testable_tissues')
if __name__=='__main__':
    phase3();phase4();phase5()
