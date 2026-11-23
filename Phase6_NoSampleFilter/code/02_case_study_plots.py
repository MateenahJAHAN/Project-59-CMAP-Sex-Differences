"""Render the four recorded Phase 6 split-library cases from the public expression matrix."""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO/'fellowship'))
from sources import prepare
from run import outside_checkout

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
        return None, None
    idx = [screen.positions[item] for item in ids]
    arr = np.asarray(screen.Z[sorted(idx)])[np.argsort(np.argsort(idx))]
    mz = pd.DataFrame(arr.T, index=screen.gene_ids)
    mz = mz.iloc[:, 0] if len(ids) == 1 else mz.mean(axis=1)
    n0 = len(mz)
    if mode != 'full':
        mz = mz[mz.abs() > 2]
        if not quiet: print(f"   |z| > 2 : {n0:,} -> {len(mz):,} genes")
        if len(mz) == 0:
            return None, None
    cl = sorted(s['cell_iname'].dropna().unique())
    sx = sorted(s['donor_sex'].dropna().unique())
    lines = f'cell line {cl[0]}' if len(cl)==1 else f'{len(cl)} cell lines'
    src = f"{len(ids)} signature{'s' if len(ids)>1 else ''}  ·  {lines}  ·  donor labels {'/'.join(sx)}"
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
        return None
    res = gseapy.prerank(rnk=V1, gene_sets=use, permutation_num=n_perm,
                         min_size=MIN_SET, max_size=100000, seed=42,
                         no_plot=True, outdir=None, threads=1, verbose=False)
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
    print(f"Reference tissue: {'lineage-matched' if tissue in home else 'other tissue'}")
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
    if len(S) == 0: return None

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
    axS.text(base+.8, 1.62, f'{base:.0f}% female-biased\namong all matched genes',
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
        f'{drug}  ×  {lineage.replace("_"," ")}   →   {tissue.replace("_"," ").replace("EBV.transformed","EBV-transformed")}\n'
        f'{"Lineage-matched GTEx tissue" if tissue in home else "Other GTEx tissue"}  ·  {src}  ·  '
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
    plt.close(fig)
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
        f"LINCS2020 drug signature\n   {drug} x {lineage[:26]}\n   {src}\n"
        f"   ranked list  {len(V1):,} genes  "
        f"({'no filter' if mode=='full' else '|z| > 2'})\n{'-'*54}\n"
        f"GTEx sex-biased gene sets\n   {tissue}\n"
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
    axz.set_xlabel('Rank in the LINCS2020 drug signature',
                   fontsize=11, fontweight='bold')
    axz.grid(axis='x', color=GRID, lw=.7, ls='--'); axz.set_axisbelow(True)
    axz.spines[['top','right']].set_visible(False)

    fig.suptitle(
        f'Enrichment of female- and male-biased GTEx gene sets in {tissue.replace("_"," ").replace("EBV.transformed","EBV-transformed")}\n'
        f'{drug}  ·  {lineage.replace("_"," ")}  ·  '
        f'{"Lineage-matched GTEx tissue" if tissue in home else "Other GTEx tissue"}  ·  {src}  ·  '
        f'{"no filter" if mode=="full" else "|z| > 2"}, {len(V1):,} genes  ·  '
        f'{len(mF)+len(mM)} GTEx sex-biased genes matched',
        fontsize=12.5, fontweight='bold', y=.985)
    plt.tight_layout(rect=[0, 0, 1, .955])
    tag = (sig_id or cell_line or 'lineage').replace(':', '_')
    if save:
        fn = f'{BASE}/figures/enrich_{drug}_{tissue[:18]}_{tag}_{mode}.png'
        plt.savefig(fn, dpi=165, bbox_inches='tight'); print("saved:", fn)
    plt.close(fig)
    print(f"\n{src}   [{mode}]")
    print(f"FEMALE  NES {stat['F'][0]:+.3f}  FDR {stat['F'][2]:.4f}  {sg(stat['F'][2])}")
    print(f"MALE    NES {stat['M'][0]:+.3f}  FDR {stat['M'][2]:.4f}  {sg(stat['M'][2])}")
    return stat

def main():
    global BASE, screen, gseapy, sig_known, gctx_sigs, cmap_to_gtex, gid_to_ensg, sbg, ENS2SYM
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', required=True, type=outside_checkout)
    parser.add_argument('--out', required=True, type=outside_checkout)
    args = parser.parse_args()
    data = args.data.expanduser().resolve()
    BASE = args.out.expanduser().resolve()
    prepare(data)
    for folder in ('figures','results','prepared'):
        (BASE/folder).mkdir(parents=True,exist_ok=True)
    os.environ.update(P59_PUBLIC=str(REPO), P59_META=str(data), P59_MATRIX=str(data),
                      P59_GTEX=str(data), P59_PREPARED=str(BASE/'prepared'))
    alternate = Path(sys.prefix)/'phase6'
    if alternate.exists():
        sys.path.insert(0,str(alternate))
    sys.path.insert(0,str(REPO/'verification/source_replay/phases2_5'))
    import replay as screen
    screen.prepare()
    screen.init()
    gseapy = screen.gseapy
    if gseapy.__version__ != '1.3.2':
        raise RuntimeError('Install the fellowship environment to replay the four cases with GSEApy 1.3.2.')
    sig_known, gctx_sigs = screen.known, set(screen.positions)
    cmap_to_gtex, gid_to_ensg = screen.map5, screen.gene_map
    sbg = pd.read_csv(data/'signif.sbgenes.txt',sep='\t')
    gi = screen.gene_info.dropna(subset=['ensembl_id']).drop_duplicates('ensembl_id')
    ENS2SYM = dict(zip(gi.ensembl_id,gi.gene_symbol))
    cases = pd.read_csv(REPO/'Phase6_NoSampleFilter/results/case_studies_split_library.csv')
    results = []
    for case in cases.itertuples(index=False):
        cell = 'BJAB' if case.scope=='single cell line BJAB' else None
        mode = 'full' if case.z_filter.startswith('none') else 'filtered'
        opts = dict(mode=mode, cell_line=cell)
        args = (case.drug,case.lineage,case.gtex_tissue)
        rank, _ = build_V1(*args[:2], **opts, quiet=True)
        fit = _gsea(rank, case.gtex_tissue, quiet=True)
        for sex, label in [('F','female'),('M','male')]:
            nes, _, q = fit['stat'][sex]
            results.append(dict(drug=case.drug,lineage=case.lineage,tissue=case.gtex_tissue,scope=case.scope,
                 library=label,ranked_genes=len(rank),NES=nes,FDR=q,archived_NES=getattr(case,label+'_NES'),
                 archived_FDR=getattr(case,label+'_FDR')))
        sex_panel(*args, **opts)
        sex_enrichment(*args, **opts)
    pd.DataFrame(results).to_csv(BASE/'results/split_library_comparisons.csv',index=False)


if __name__ == '__main__':
    main()
