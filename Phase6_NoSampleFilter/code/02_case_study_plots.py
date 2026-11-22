"""
Phase 6 - step 2: the per-case functions used for the single-pair figures.

IMPORTANT - this uses a DIFFERENT gene-set definition from 01_batch_abcd.py:

  01_batch_abcd.py   run_gsea   ->  44 sets, one per tissue,
                                    all sex-biased genes of that tissue together
  this file          _gsea      ->  88 sets, each tissue split into a
                                    FEMALE half and a MALE half by effsize sign

So NES and FDR from the batch table and from these plots are NOT comparable.
They are two different tests on the same data.

Verbatim from Phase6_NoSampleFilter/notebook/phase6_no_sample_filter.ipynb
(cells 34 and 36).
"""

# ------- CELL A: gctx + reference + core functions -------
# ============================================================
#  CELL A  —  restore: gctx + reference + core functions
#  download me 15-20 min lagenge
# ============================================================
!pip install cmapPy gseapy -q

from google.colab import drive
drive.mount('/content/drive', force_remount=True)

import os, json, glob, numpy as np, pandas as pd, gseapy
from cmapPy.pandasGEXpress.parse import parse

BASE = '/content/drive/MyDrive/Project59'
REF  = f'{BASE}/data/reference'

# ---------- gctx ----------
os.makedirs('/content/cmap', exist_ok=True)
gctx_file = '/content/cmap/level5_beta_trt_cp_n720216x12328.gctx'
if os.path.exists(gctx_file):
    print(f"gctx already here: {os.path.getsize(gctx_file)/1e9:.1f} GB")
else:
    print("downloading 33 GB — tab band mat karna\n")
    !cd /content/cmap && curl -O https://cmapreg.s3.amazonaws.com/galenAI/level5_beta_trt_cp_n720216x12328.gctx

gctx_sigs = set(parse(gctx_file, col_meta_only=True).index)
print(f"gctx signatures : {len(gctx_sigs):,}")

# ---------- reference ----------
cmap_to_gtex = {k: set(v) for k, v in json.load(open(f'{REF}/cmap_to_gtex.json')).items()}
gid_to_ensg  = json.load(open(f'{REF}/gid_to_ensg.json'))
top50_names  = json.load(open(f'{REF}/top50_names.json'))
genes_info   = pd.read_csv(f'{BASE}/data/cmap/geneinfo.csv')
sig_known    = pd.read_csv(f'{BASE}/results/sig_known/sig_known.csv', low_memory=False)

# ---------- sbg (GTEx significant sex-biased genes) ----------
cands = ([f'{BASE}/data/gtex/sbgenes_raw.csv'] +
         glob.glob(f'{BASE}/data/gtex/*sbgene*') +
         glob.glob(f'{BASE}/data/**/*signif*sbgene*', recursive=True))
sbg = None
for c in cands:
    if os.path.exists(c):
        sbg = pd.read_csv(c, sep='\t' if c.endswith(('.txt','.tsv')) else ',')
        print(f"sbg loaded from : {c}  ({len(sbg):,} rows)"); break
if sbg is None:
    raise FileNotFoundError("sbg nahi mila — data/gtex folder ki listing bhejo")

# ---------- ENS2SYM ----------
gi = genes_info[['gene_id','gene_symbol']].astype({'gene_id':str}).copy()
gi['ens'] = gi['gene_id'].map(gid_to_ensg)
ENS2SYM = dict(gi.dropna(subset=['ens']).drop_duplicates('ens')[['ens','gene_symbol']].values)
print(f"ENS2SYM         : {len(ENS2SYM):,} genes")

# ---------- palette ----------
F_COL, M_COL = '#C2185B', '#1565C0'
INK, GRID, GREY, MUTED = '#1F2427', '#D8DEE1', '#BDC5C9', '#6E777C'
MIN_SET = 15

# ============================================================
#  CORE FUNCTIONS
# ============================================================
def build_V1(drug, lineage, mode='full', sig_id=None, cell_line=None, quiet=False):
    s = sig_known[(sig_known['cmap_name'] == drug) &
                  (sig_known['cell_lineage'] == lineage)]
    s = s[s['sig_id'].isin(gctx_sigs)]
    if sig_id    is not None: s = s[s['sig_id']     == sig_id]
    elif cell_line is not None: s = s[s['cell_iname'] == cell_line]
    ids = s['sig_id'].tolist()
    if not ids:
        print("   koi signature nahi mila"); return None, None
    mz = parse(gctx_file, cid=ids).data_df
    mz = mz.iloc[:, 0] if len(ids) == 1 else mz.mean(axis=1)
    n0 = len(mz)
    if mode != 'full':
        mz = mz[mz.abs() > 2]
        if not quiet: print(f"   |z| > 2 : {n0:,} -> {len(mz):,} genes")
        if len(mz) == 0:
            print("   filter ke baad kuch nahi bacha"); return None, None
    cl = sorted(s['cell_iname'].dropna().unique())
    sx = sorted(s['donor_sex'].dropna().unique())
    src = (f"{len(ids)} signature{'s' if len(ids)>1 else ''}  ·  "
           f"cell line {cl[0] if len(cl)==1 else f'{len(cl)} lines'}  ·  "
           f"donor {sx[0] if len(sx)==1 else '/'.join(sx)}")
    d = pd.DataFrame({'gene_id': mz.index.astype(str), 'z': mz.values})
    d['ensembl_id'] = d['gene_id'].map(gid_to_ensg)
    V1 = (d.dropna(subset=['ensembl_id'])
            .groupby('ensembl_id')['z'].mean().sort_values(ascending=False))
    return V1, src


def split_sets(tissue, effsize_min=0.0):
    t = sbg[sbg['tissue'] == tissue].copy()
    t['ens'] = t['gene'].astype(str).str.split('.').str[0]
    eff = t.groupby('ens')['effsize'].mean()
    eff = eff[eff.abs() >= effsize_min]
    return {f'{tissue}__FEMALE': set(eff[eff > 0].index),
            f'{tissue}__MALE'  : set(eff[eff < 0].index)}


def _running_es(rnk, gene_set):
    hits = np.isin(rnk.index.values, list(gene_set))
    w    = np.abs(rnk.values)
    if hits.sum() == 0 or (~hits).sum() == 0:
        return np.zeros(len(rnk)), hits
    step = np.where(hits, w / w[hits].sum(), -1.0 / (~hits).sum())
    return np.cumsum(step), hits


def _leading_edge(cum, hits):
    i   = int(np.argmax(np.abs(cum)))
    es  = float(cum[i])
    idx = np.where(hits)[0]
    return i, es, (idx[idx >= i] if es < 0 else idx[idx <= i])


def _gsea(V1, tissue, n_perm=1000, quiet=False):
    halves = split_sets(tissue)
    Fk, Mk = f'{tissue}__FEMALE', f'{tissue}__MALE'
    idx = set(V1.index)
    use = {k: v for k, v in halves.items() if len(v & idx) >= MIN_SET}
    if not quiet:
        for k in halves:
            n = len(halves[k] & idx)
            print(f"   {k.split('__')[1]:<7} matched {n}"
                  + ("" if k in use else f"   SKIP (< {MIN_SET})"))
    if not use:
        print(f"   dono halves {MIN_SET} se chhote"); return None
    res = gseapy.prerank(rnk=V1, gene_sets=use, permutation_num=n_perm,
                         min_size=MIN_SET, max_size=100000, seed=42,
                         no_plot=True, verbose=False)
    r, stat = res.res2d.set_index('Term'), {}
    for key, lab in [(Fk, 'F'), (Mk, 'M')]:
        stat[lab] = ((float(r.loc[key, 'NES']), float(r.loc[key, 'NOM p-val']),
                      float(r.loc[key, 'FDR q-val'])) if key in r.index
                     else (np.nan, np.nan, np.nan))
    gs = V1.index.values
    def le(S):
        cum, hit = _running_es(V1, S)
        return np.array([], int) if hit.sum() == 0 else _leading_edge(cum, hit)[2]
    return dict(gs=gs,
                kF=np.isin(gs, list(halves[Fk])),
                kM=np.isin(gs, list(halves[Mk])),
                stat=stat, leF=le(halves[Fk]), leM=le(halves[Mk]))


def check_sizes(drug, lineage, tissue, sig_id=None, cell_line=None):
    h = split_sets(tissue)
    F, M = h[f'{tissue}__FEMALE'], h[f'{tissue}__MALE']
    home = set(cmap_to_gtex.get(lineage, []))
    print(f"\n{tissue}  library {len(F)+len(M):,}  (F {len(F):,} / M {len(M):,})")
    print(f"category : {'A  HOME' if tissue in home else 'B  CROSS'}")
    print(f"\n{drug} × {lineage} → {tissue}")
    print(f"{'mode':<12}{'ranked list':>14}{'F match':>10}{'M match':>10}   verdict")
    print("-"*62)
    for m in ['full', 'filtered']:
        V1, _ = build_V1(drug, lineage, m, sig_id, cell_line, quiet=True)
        if V1 is None: continue
        i = set(V1.index); nf, nm = len(F & i), len(M & i)
        ok = ('both halves' if min(nf, nm) >= MIN_SET else
              'one half'    if max(nf, nm) >= MIN_SET else 'NEITHER — skip')
        print(f"{m:<12}{len(V1):>14,}{nf:>10}{nm:>10}   {ok}")
    print(f"\n(ek half ko kam se kam {MIN_SET} genes chahiye)")

print("\nCELL A ready")

# ------- CELL B: sex_panel + sex_enrichment -------
# ============================================================
#  CELL B  —  sex_panel + sex_enrichment
# ============================================================
import numpy as np, pandas as pd, gseapy
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def sex_panel(drug, lineage, tissue, mode='filtered', sig_id=None,
              cell_line=None, n_perm=1000, save=True):
    V1, src = build_V1(drug, lineage, mode, sig_id, cell_line, quiet=True)
    if V1 is None: return None
    G = _gsea(V1, tissue, n_perm, quiet=True)
    if G is None: return None

    gs, stat = G['gs'], G['stat']
    Fset, Mset = set(gs[G['kF']]), set(gs[G['kM']])
    le = ({V1.index[i] for i in G['leF']} | {V1.index[i] for i in G['leM']})

    t = sbg[sbg['tissue'] == tissue].copy()
    t['ens'] = t['gene'].astype(str).str.split('.').str[0]
    eff = t.groupby('ens')['effsize'].mean()

    D = pd.DataFrame({'ens': V1.index, 'z': V1.values})
    D['rank'] = np.arange(1, len(D) + 1)
    D['gene'] = D['ens'].map(lambda g: ENS2SYM.get(g, g))
    D['sex']  = np.where(D['ens'].isin(Fset), 'F',
                np.where(D['ens'].isin(Mset), 'M', ''))
    D['eff']  = D['ens'].map(eff)
    D['LE']   = D['ens'].isin(le)

    S = D[D['sex'] != ''].copy()
    if len(S) == 0: print("koi sex-biased gene match nahi hua"); return None
    if len(S) > 400: print(f"WARNING: {len(S)} genes — figure bahut lamba hoga")

    blocks = [('UP by drug  \u2640  female-biased',
               S[(S.sex=='F') & (S.z>0)].sort_values('z', ascending=False), F_COL),
              ('UP by drug  \u2642  male-biased',
               S[(S.sex=='M') & (S.z>0)].sort_values('z', ascending=False), M_COL),
              ('DOWN by drug  \u2640  female-biased',
               S[(S.sex=='F') & (S.z<0)].sort_values('z'), F_COL),
              ('DOWN by drug  \u2642  male-biased',
               S[(S.sex=='M') & (S.z<0)].sort_values('z'), M_COL)]
    upF, upM, dnF, dnM = [len(b[1]) for b in blocks]
    nF, nM = int((S.sex=='F').sum()), int((S.sex=='M').sum())
    base = 100 * nF / len(S)
    rows = max(1, max(len(b[1]) for b in blocks))

    fig = plt.figure(figsize=(18.5, max(8.5, .175*rows + 2.6)))
    gsp = fig.add_gridspec(2, 4, height_ratios=[1, max(6, .175*rows)],
                           hspace=.10, wspace=.50)
    axS = fig.add_subplot(gsp[0, :])
    for i, (lab, f, m) in enumerate([('DOWN-regulated', dnF, dnM),
                                     ('UP-regulated',   upF, upM)]):
        tot = f + m; pf = 100*f/tot if tot else 0
        axS.barh(i, pf, .5, color=F_COL, edgecolor='white')
        axS.barh(i, 100-pf, .5, left=pf, color=M_COL, edgecolor='white')
        if f: axS.text(pf/2, i, f'{f}  ({pf:.0f}%)', ha='center', va='center',
                       color='white', fontsize=10.5, fontweight='bold')
        if m: axS.text(pf+(100-pf)/2, i, f'{m}  ({100-pf:.0f}%)', ha='center',
                       va='center', color='white', fontsize=10.5, fontweight='bold')
    axS.axvline(base, ls='--', lw=1.6, color=INK)
    axS.text(base+.8, 1.62, f'expected {base:.0f}% female\nif there were no pattern',
             fontsize=8.6, color=INK, va='top')
    axS.set_yticks([0,1]); axS.set_yticklabels(['DOWN-regulated','UP-regulated'],
                                               fontsize=10.5, fontweight='bold')
    axS.set_xlim(0,100); axS.set_ylim(-.5, 2.0); axS.set_xticks([])
    axS.tick_params(left=False)
    axS.spines[['top','right','left','bottom']].set_visible(False)
    axS.set_title('Of the sex-biased genes, what share is female-biased?',
                  fontsize=10.5, fontweight='bold', loc='left', pad=6)

    pad = S['z'].abs().max() * 2.15
    for j, (lab, B, col) in enumerate(blocks):
        ax = fig.add_subplot(gsp[1, j]); B = B.reset_index(drop=True)
        if len(B):
            y = np.arange(len(B))
            ax.barh(y, B['z'], .78, color=col, edgecolor='white', linewidth=.45)
            for i, r in enumerate(B.itertuples()):
                ax.text(r.z + (pad*.03 if r.z > 0 else -pad*.03), i,
                        f"z {r.z:+.2f}  eff {r.eff:+.3f}  #{r.rank:,}",
                        va='center', ha='left' if r.z > 0 else 'right',
                        fontsize=6.0, color=col)
            ax.set_yticks(y); ax.set_yticklabels(B['gene'], fontsize=6.6)
            for tl, le_ in zip(ax.get_yticklabels(), B['LE']):
                tl.set_color(col)
                if le_: tl.set_fontweight('bold')
        ax.axvline(0, color=INK, lw=1.1)
        ax.set_xlim(-pad, pad); ax.set_ylim(rows-.4, -.6)
        ax.set_title(f'{lab}\n{len(B)} genes', fontsize=9.5,
                     fontweight='bold', color=col, loc='left', pad=8)
        ax.set_xlabel('z-score', fontsize=9, fontweight='bold')
        ax.grid(axis='x', color=GRID, lw=.7, ls='--'); ax.set_axisbelow(True)
        ax.spines[['top','right','left']].set_visible(False)
        ax.tick_params(left=False, labelsize=7.5)

    g  = lambda k,i: stat[k][i] if k in stat else float('nan')
    sg = lambda q: 'SIGNIFICANT' if (q == q and q < 0.05) else 'not significant'
    home = set(cmap_to_gtex.get(lineage, []))
    fig.suptitle(
        f'{drug}  ×  {lineage}   →   {tissue}\n'
        f'Category {"A  HOME" if tissue in home else "B  CROSS"}  ·  {src}  ·  '
        f'{"full" if mode=="full" else "|z| > 2"}, {len(V1):,} genes  ·  '
        f'ALL {len(S)} sex-biased genes shown  (\u2640 {nF} · \u2642 {nM})\n'
        f'FEMALE  NES {g("F",0):+.3f}  FDR {g("F",2):.4f}  {sg(g("F",2))}      |      '
        f'MALE  NES {g("M",0):+.3f}  FDR {g("M",2):.4f}  {sg(g("M",2))}      |      '
        f'bold = leading edge', fontsize=11, fontweight='bold', y=.998)
    plt.tight_layout(rect=[0, 0, 1, .965])
    tag = (sig_id or cell_line or 'lineage').replace(':', '_')
    if save:
        fn = f'{BASE}/figures/panel_{drug}_{tissue[:16]}_{tag}_{mode}.png'
        plt.savefig(fn, dpi=165, bbox_inches='tight'); print("saved:", fn)
    plt.show()
    print(f"\n{src}   [{mode}]   ranked list {len(V1):,}")
    print(f"sex-biased matched : {len(S)}   (\u2640 {nF} / \u2642 {nM})")
    print(f"   UP   \u2640 {upF:>4}   \u2642 {upM:>4}")
    print(f"   DOWN \u2640 {dnF:>4}   \u2642 {dnM:>4}   -> DOWN is "
          f"{100*dnF/max(1,dnF+dnM):.0f}% female vs {base:.0f}% expected")
    S[['gene','ens','sex','z','rank','eff','LE']].sort_values('z', ascending=False)\
     .to_csv(f'{BASE}/results/sexgenes_{drug}_{tissue[:16]}_{tag}_{mode}.csv', index=False)
    return S


def sex_enrichment(drug, lineage, tissue, mode='filtered', sig_id=None,
                   cell_line=None, n_perm=1000, max_labels=30, save=True):
    V1, src = build_V1(drug, lineage, mode, sig_id, cell_line, quiet=True)
    if V1 is None: return None
    halves = split_sets(tissue)
    Fk, Mk = f'{tissue}__FEMALE', f'{tissue}__MALE'
    idx = set(V1.index); mF, mM = halves[Fk] & idx, halves[Mk] & idx
    libN = len(halves[Fk]) + len(halves[Mk])
    print(f"library {libN:,}  ->  matched {len(mF)+len(mM):,}  (F {len(mF)} / M {len(mM)})")

    G = _gsea(V1, tissue, n_perm, quiet=True)
    if G is None: return None
    stat = G['stat']

    def walk(S_):
        cum, hit = _running_es(V1, S_)
        if hit.sum() == 0: return np.zeros(len(V1)), hit, 0, np.array([], int)
        i, es, le = _leading_edge(cum, hit); return cum, hit, i, le
    cumF, hitF, iF, leF = walk(halves[Fk])
    cumM, hitM, iM, leM = walk(halves[Mk])

    z, x = V1.values, np.arange(1, len(V1)+1)
    zero = int((z > 0).sum())
    home = set(cmap_to_gtex.get(lineage, []))
    sg = lambda q: 'SIGNIFICANT' if (q == q and q < 0.05) else 'not significant'

    fig = plt.figure(figsize=(19, 11))
    gsp = fig.add_gridspec(4, 1, height_ratios=[3.1, 1.15, 1.15, 1.5], hspace=.18)
    ax = fig.add_subplot(gsp[0])
    ax.plot(x, cumF, color=F_COL, lw=1.5, label='female-biased set')
    ax.plot(x, cumM, color=M_COL, lw=1.5, label='male-biased set')
    ax.axhline(0, color=INK, lw=1.0)
    if hitF.sum(): ax.axvline(iF+1, color=F_COL, ls='--', lw=1.1, alpha=.7)
    if hitM.sum(): ax.axvline(iM+1, color=M_COL, ls='--', lw=1.1, alpha=.7)
    ax.set_ylabel('Running enrichment score', fontsize=11, fontweight='bold')
    ax.set_xlim(0, len(V1)); ax.set_xticklabels([])
    ax.grid(color=GRID, lw=.7, ls='--'); ax.set_axisbelow(True)
    ax.spines[['top','right']].set_visible(False)
    ax.legend(loc='lower left', fontsize=9.5, framealpha=.95)
    ax.text(.012, .045,
        f"SOURCE (V1) = the drug\n   {drug} x {lineage[:26]}\n   {src}\n"
        f"   ranked list  {len(V1):,} genes  "
        f"({'no filter' if mode=='full' else '|z| > 2'})\n{'-'*54}\n"
        f"TARGET (V2) = GTEx sex-biased genes\n   {tissue}\n"
        f"   library {libN:,}  (F {len(halves[Fk]):,} / M {len(halves[Mk]):,})\n"
        f"   matched {len(mF)+len(mM):,}  (F {len(mF)} / M {len(mM)})\n{'-'*54}\n"
        f"FEMALE  NES {stat['F'][0]:+.3f}  FDR {stat['F'][2]:.4f}  [{sg(stat['F'][2])}]\n"
        f"        leading edge {len(leF)} of {int(hitF.sum())}\n"
        f"MALE    NES {stat['M'][0]:+.3f}  FDR {stat['M'][2]:.4f}  [{sg(stat['M'][2])}]\n"
        f"        leading edge {len(leM)} of {int(hitM.sum())}",
        transform=ax.transAxes, fontsize=8.0, family='monospace',
        va='bottom', ha='left',
        bbox=dict(boxstyle='round,pad=0.55', facecolor='#FAFAFA',
                  edgecolor=INK, linewidth=1.0))

    for k, (hit, le, col, lab, tot) in enumerate([
            (hitF, leF, F_COL, 'FEMALE-biased', len(mF)),
            (hitM, leM, M_COL, 'MALE-biased',   len(mM))]):
        axb = fig.add_subplot(gsp[1+k])
        axb.vlines(np.where(hit)[0]+1, 0, 1, color=col, lw=.9, alpha=.55)
        if len(le):
            axb.vlines(le+1, 0, 1, color=col, lw=2.0)
            sel = le if len(le) <= max_labels else le[np.argsort(-np.abs(z[le]))][:max_labels]
            for i in sorted(sel):
                axb.text(i+1, 1.06, f"{ENS2SYM.get(V1.index[i], V1.index[i])} "
                         f"{'\u2640' if col==F_COL else '\u2642'}",
                         rotation=90, fontsize=5.8, color=col, ha='center', va='bottom')
        axb.set_xlim(0, len(V1)); axb.set_ylim(0,1)
        axb.set_yticks([]); axb.set_xticklabels([])
        axb.set_ylabel(f'{lab}\n{tot} genes\nleading edge {len(le)}',
                       fontsize=8.6, fontweight='bold', color=col, rotation=0,
                       ha='right', va='center', labelpad=14)
        axb.spines[['top','right','left','bottom']].set_visible(False)
        axb.set_facecolor('#F4F6F7')

    axz = fig.add_subplot(gsp[3])
    axz.fill_between(x, z, 0, color=GREY, lw=0)
    axz.axhline(0, color=INK, lw=1.0); axz.axvline(zero, color=INK, ls='--', lw=1.1)
    axz.text(zero + len(V1)*.006, z.min()*.55, f'zero cross at rank {zero:,}',
             fontsize=8.6, color=INK)
    axz.text(len(V1)*.012, z.max()*.72, 'left  =  up-regulated by the drug',
             fontsize=9.5, fontweight='bold', color=INK)
    axz.text(len(V1)*.988, z.min()*.72, 'right  =  down-regulated by the drug',
             fontsize=9.5, fontweight='bold', color=INK, ha='right')
    axz.set_xlim(0, len(V1))
    axz.set_ylabel('z-score', fontsize=10, fontweight='bold')
    axz.set_xlabel('Rank in the drug signature  (SOURCE / V1)',
                   fontsize=11, fontweight='bold')
    axz.grid(axis='x', color=GRID, lw=.7, ls='--'); axz.set_axisbelow(True)
    axz.spines[['top','right']].set_visible(False)

    fig.suptitle(
        f'Are the genes this drug modulates enriched for sex-biased genes of {tissue}?\n'
        f'{drug}  ·  {lineage}  ·  '
        f'Category {"A  HOME" if tissue in home else "B  CROSS"}  ·  {src}  ·  '
        f'{"no filter" if mode=="full" else "|z| > 2"}, {len(V1):,} genes  ·  '
        f'{len(mF)+len(mM)} GTEx sex-biased genes matched',
        fontsize=12.5, fontweight='bold', y=.985)
    plt.tight_layout(rect=[0, 0, 1, .955])
    tag = (sig_id or cell_line or 'lineage').replace(':', '_')
    if save:
        fn = f'{BASE}/figures/enrich_{drug}_{tissue[:18]}_{tag}_{mode}.png'
        plt.savefig(fn, dpi=165, bbox_inches='tight'); print("saved:", fn)
    plt.show()
    print(f"\n{src}   [{mode}]")
    print(f"FEMALE  NES {stat['F'][0]:+.3f}  FDR {stat['F'][2]:.4f}  {sg(stat['F'][2])}")
    print(f"MALE    NES {stat['M'][0]:+.3f}  FDR {stat['M'][2]:.4f}  {sg(stat['M'][2])}")
    return stat

print("CELL B ready")
