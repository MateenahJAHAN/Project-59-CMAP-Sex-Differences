"""Figures (matplotlib; Okabe-Ito colours validated for CVD separation; PNG 300 dpi and PDF)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import gaussian_kde
from common import OUT, Y8, load_cells

FIG = OUT + '/figures'
os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({'font.size': 15, 'axes.titlesize': 16, 'axes.labelsize': 15, 'xtick.labelsize': 13,
                     'ytick.labelsize': 13, 'legend.fontsize': 13, 'axes.spines.top': False,
                     'axes.spines.right': False, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'axes.edgecolor': '#555555', 'axes.labelcolor': '#222222', 'text.color': '#222222'})
COL = {'Male': '#0072B2', 'Female': '#D55E00', 'Unknown': '#009E73'}
MRK = {'Male': 'o', 'Female': 's', 'Unknown': 'D'}
INK2 = '#555555'


def save(fig, name):
    fig.savefig(f'{FIG}/{name}.png', dpi=300, bbox_inches='tight')
    fig.savefig(f'{FIG}/{name}.pdf', bbox_inches='tight')
    plt.close(fig)


S = pd.read_csv(OUT + '/partA_local_states.tsv', sep='\t')
LP = pd.read_csv(OUT + '/partA_local_lineplate_scores.tsv', sep='\t')
cut = pd.read_csv(OUT + '/partA_local_cutoffs.tsv', sep='\t').set_index('score')['cutoff']
cutl = pd.read_csv(OUT + '/partA_local_cutoffs_linelevel.tsv', sep='\t').set_index('score')['cutoff']

# ---------------- Figure 1: XIST against Y score
fig = plt.figure(figsize=(17, 7.6))
gs = fig.add_gridspec(2, 3, width_ratios=[4, 0.8, 5], height_ratios=[0.8, 4], wspace=0.08, hspace=0.06)
ax = fig.add_subplot(gs[1, 0]); axt = fig.add_subplot(gs[0, 0], sharex=ax); axr = fig.add_subplot(gs[1, 1], sharey=ax)
for a_, nm in [('Male', 'Male'), ('Female', 'Female'), ('Unknown', 'Unknown')]:
    d = LP[LP.annotation == a_]
    ax.scatter(d.y_score, d.XIST, s=14, c=COL[a_], marker=MRK[a_], alpha=0.45, linewidths=0, label=f'{nm} ({len(d)})')
ax.axvline(cut['y_score'], color='#222222', lw=1.6); ax.axhline(cut['XIST'], color='#222222', lw=1.6)
ax.set_xlabel('Y score, mean log2(CPM+1) of 8 Y genes'); ax.set_ylabel('XIST, log2(CPM+1)')
ax.legend(title='line x plate, annotation', loc='center right', bbox_to_anchor=(1.0, 0.6), frameon=False, markerscale=1.8)
for v, axm, orient in [(LP.y_score, axt, 'x'), (LP.XIST, axr, 'y')]:
    k = gaussian_kde(v, bw_method='scott'); g = np.linspace(v.min(), v.max(), 1024)
    if orient == 'x':
        axm.fill_between(g, k(g), color='#999999', alpha=0.5, lw=0); axm.axvline(cut['y_score'], color='#222222', lw=1.6)
    else:
        axm.fill_betweenx(g, k(g), color='#999999', alpha=0.5, lw=0); axm.axhline(cut['XIST'], color='#222222', lw=1.6)
    axm.axis('off')
axt.set_title(f'A. {len(LP)} line x plate pseudobulks\ncut-offs at the KDE minimum: Y {cut["y_score"]:.2f}, XIST {cut["XIST"]:.2f}',
              loc='left', fontsize=14)
ax2 = fig.add_subplot(gs[:, 2])
for a_ in ['Male', 'Female', 'Unknown']:
    d = S[S.annotation == a_]
    ax2.scatter(d.y_score, d.XIST, s=90, c=COL[a_], marker=MRK[a_], edgecolors='white', linewidths=1.2,
                label=f'{a_} ({len(d)})', zorder=3)
ax2.axvline(cut['y_score'], color='#222222', lw=1.6, label='primary cut-offs')
ax2.axhline(cut['XIST'], color='#222222', lw=1.6)
ax2.axvline(cutl['y_score'], color='#222222', lw=1.4, ls='--', label='line-level cut-offs')
ax2.axhline(cutl['XIST'], color='#222222', lw=1.4, ls='--')
lab = ['A498', 'RKO', 'SW48', 'HT-29', 'NCI-H2347', 'NCI-H2122', 'COLO 205', 'CFPAC-1', 'NCI-H2030', 'NCI-H596']
off = {'A498': (-9, 0, 'right'), 'RKO': (9, 2, 'left'), 'SW48': (9, -12, 'left'), 'HT-29': (9, 4, 'left'),
       'NCI-H2347': (-9, -8, 'right'), 'NCI-H2122': (-9, 4, 'right'), 'COLO 205': (6, 8, 'left'),
       'CFPAC-1': (0, 10, 'center'), 'NCI-H2030': (0, -20, 'center'), 'NCI-H596': (9, 2, 'left')}
for r in S[S.cell_name.isin(lab)].itertuples():
    dx, dy, ha = off.get(r.cell_name, (6, 6, 'left'))
    ax2.annotate(r.cell_name, (r.y_score, r.XIST), xytext=(dx, dy), textcoords='offset points', ha=ha,
                 fontsize=12, color='#222222')
ax2.set_xlim(0.3, 7.2)
ax2.set_xlabel('Y score (line level)'); ax2.set_ylabel('XIST (line level)')
n_st = S.state.value_counts()
ax2.set_title(f'B. {len(S)} lines after the foreign-cell filter\n' + ', '.join(f'{k} {v}' for k, v in n_st.items()),
              loc='left', fontsize=14)
ax2.legend(loc='center right', bbox_to_anchor=(1.0, 0.66), frameon=False)
ax2.text(5.9, 2.35, f"{int((S.state == 'XY').sum())} XY lines", ha='center', fontsize=12, color=INK2)
save(fig, 'fig1_xist_vs_y_score')

# ---------------- Figure 2: arm scores
fig, ax = plt.subplots(figsize=(10.5, 8.5))
for a_ in ['Male', 'Female', 'Unknown']:
    d = S[S.annotation == a_]
    ax.scatter(d.yp_score, d.yq_score, s=90, c=COL[a_], marker=MRK[a_], edgecolors='white', linewidths=1.2,
               label=f'{a_} ({len(d)})', zorder=3)
ax.axvline(cut['yp_score'], color='#222222', lw=1.6); ax.axhline(cut['yq_score'], color='#222222', lw=1.6)
x0, x1 = 0.8, S.yp_score.max() + 0.4
y0, y1 = 0.6, S.yq_score.max() + 0.7
ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
for txt, x, y, ha, va in [('Y intact', x1 - 0.1, y1 - 0.05, 'right', 'top'), ('Yp lost', x0 + 0.1, y1 - 0.05, 'left', 'top'),
                          ('Yq lost', x1 - 0.1, cut['yq_score'] - 0.1, 'right', 'top'),
                          ('Y lost', cut['yp_score'] - 0.1, y0 + 0.1, 'right', 'bottom')]:
    ax.text(x, y, txt, ha=ha, va=va, fontsize=15, color=INK2, style='italic')
lab_off = {'NCI-H2030': (-10, 10, 'right'), 'CFPAC-1': (10, 4, 'left'), 'COLO 205': (10, 4, 'left'),
           'HCT15': (10, 4, 'left'), 'NCI-H2122': (-10, 8, 'right')}
for r in S[((S.annotation == 'Male') & ~S.arm_state.isin(['Y intact', 'Y lost'])) |
           ((S.annotation != 'Male') & (S.arm_state != 'Y lost'))].itertuples():
    dx, dy, ha = lab_off.get(r.cell_name, (10, 4, 'left'))
    txt = f'{r.cell_name} ({r.arm_state})' if r.annotation == 'Male' else f'{r.cell_name} ({r.annotation}, {r.arm_state})'
    ax.annotate(txt, (r.yp_score, r.yq_score), xytext=(dx, dy), textcoords='offset points', ha=ha, fontsize=12)
ax.set_xlabel('Yp score, mean log2(CPM+1) of RPS4Y1, ZFY')
ax.set_ylabel('Yq score, mean log2(CPM+1) of 6 Yq genes')
ax.set_title(f'Arm scores per line; cut-offs Yp {cut["yp_score"]:.2f}, Yq {cut["yq_score"]:.2f} (line x plate KDE)', loc='left', fontsize=14)
ax.legend(loc='upper left', bbox_to_anchor=(0.0, 0.93), frameon=False, title='annotation')
save(fig, 'fig2_arm_scores')

# ---------------- Figure 3: per-cell Y UMIs against depth
cells = load_cells(['drug', 'cell_line_id', 'cell_name', 'total_umi', 'y8', 'XIST'])
curve = pd.read_csv(OUT + '/partA_local_depth_curve.tsv', sep='\t')
mos = pd.read_csv(OUT + '/partA_local_mosaic.tsv', sep='\t').set_index(['cell_name', 'marker'])
three = ['A549', 'NCI-H460', 'LoVo']
fig, axes = plt.subplots(2, 3, figsize=(18, 10.5), gridspec_kw=dict(height_ratios=[3, 2.4], hspace=0.32, wspace=0.22))
for j, nm in enumerate(three):
    lid = S.loc[S.cell_name == nm, 'cellosaurus_id'].iloc[0]
    d = cells[cells.cell_line_id == lid]
    d = d[d.XIST < 3]  # the foreign-cell filter of these XIST-negative lines
    assert len(d) == int(S.loc[S.cell_name == nm, 'n_vehicle_cells'].iloc[0])
    ax = axes[0, j]
    x = np.log10(d.total_umi.to_numpy()); y = d.y8.to_numpy()
    hb = ax.hexbin(x, np.log2(y + 1), gridsize=(46, 22), norm=LogNorm(), cmap='Blues', mincnt=1, linewidths=0.2)
    yt = [n for n in [0, 1, 3, 7, 15, 31, 63] if np.log2(n + 1) <= np.log2(y.max() + 1) + 0.5]
    ax.set_yticks([np.log2(n + 1) for n in yt]); ax.set_yticklabels([str(n) for n in yt])
    ax.set_xlabel('depth, log10 UMIs per cell'); ax.set_ylabel('Y UMIs per cell (8 genes)')
    ax.set_title(f'{nm}\n{len(d):,} vehicle cells; {np.mean(y == 0):.1%} with 0 Y UMIs', loc='left', fontsize=14)
    for t in [3000, 8000]:
        ax.axvline(np.log10(t), color=INK2, lw=1.2, ls=':')
        ax.text(np.log10(t), ax.get_ylim()[1], f' {t:,}', ha='left', va='top', fontsize=11, color=INK2)
    cb = fig.colorbar(hb, ax=ax, pad=0.01); cb.set_label('cells')
    ax = axes[1, j]
    c = curve[(curve.cell_name == nm) & (curve.marker == 'Y') & (curve.n_cells >= 30)]
    ax.plot(c.median_depth, c.zero_frac, 'o-', color='#0072B2', lw=2, ms=8, label='observed')
    ax.plot(c.median_depth, c.homogeneous_nb_zero_frac, 's--', color='#D55E00', lw=2, ms=7,
            label='every cell positive (NB part of (b))')
    ax.plot(c.median_depth, c.homogeneous_poisson_zero_frac, '^:', color='#009E73', lw=2, ms=7,
            label='every cell positive (Poisson, rate r+)')
    ax.set_xscale('log'); ax.set_ylim(-0.02, 1.0)
    ax.set_xlabel('depth bin, median UMIs per cell'); ax.set_ylabel('share of cells with 0 Y UMIs')
    m = mos.loc[(nm, 'Y')]
    ax.set_title(f'(b) share at 2,000 UMIs {m.b_pi_2000:.3f}\nrange over (a) to (c): {m.share_range_lo:.2f} to {m.share_range_hi:.2f}',
                 loc='left', fontsize=13)
    if j == 0:
        ax.legend(frameon=False, fontsize=12, loc='upper right')
save(fig, 'fig3_y_umis_vs_depth_A549_H460_LoVo')

# ---------------- Figure 4: agreement
ag = pd.read_csv(OUT + '/partA_local_agreement.tsv', sep='\t')
order = ['XY', 'XX-like', 'Y-lost', 'XIST-lost', 'discordant', 'ambiguous']
fig, axes = plt.subplots(1, 3, figsize=(21, 7.2), gridspec_kw=dict(wspace=0.55))
summ = pd.read_csv(OUT + '/partA_local_agreement_summary.tsv', sep='\t').set_index('comparison')
pilot_cells = int(ag.pilot_n_cells.sum())
for ax, col, title, key in [(axes[0], 'C_state', 'Task C (external RNA-seq)', 'Task C, state (primary cut-offs), ambiguous counted as a state [A1]'),
                            (axes[1], 'A_state', 'Task A (CMap lines)', 'Task A, state'),
                            (axes[2], 'pilot_state', f'pilot table ({pilot_cells:,} cells)', 'pilot table, state')]:
    d = ag[ag[col].notna() & ~ag[col].isin(['not_measured', 'no_data'])]
    cols = [s for s in order if s in set(d[col])]
    rows = [s for s in order if s in set(d.state)]
    M = pd.crosstab(d.state, d[col]).reindex(index=rows, columns=cols, fill_value=0)
    ax.imshow(M.to_numpy(), cmap='Blues', vmin=0, vmax=max(1, M.to_numpy().max()) * 1.15, aspect='auto')
    for i in range(M.shape[0]):
        for k in range(M.shape[1]):
            v = M.iat[i, k]
            if v:
                ax.text(k, i, str(v), ha='center', va='center', fontsize=15,
                        color='white' if v > M.to_numpy().max() * 0.55 else '#222222',
                        fontweight='bold' if rows[i] == cols[k] else 'normal')
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=35, ha='right')
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows)
    ax.set_xlabel(title); ax.set_ylabel('this caller (Tahoe vehicle cells)')
    s = summ.loc[key]
    ax.set_title(f'{int(s.n_agree)} of {int(s.n_lines_with_both)} lines agree ({s.share:.1%})', loc='left', fontsize=14)
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
save(fig, 'fig4_agreement')

# ---------------- Figure 5: mosaic shares by estimator
mo = pd.read_csv(OUT + '/partA_local_mosaic.tsv', sep='\t')
mo['lab'] = mo.cell_name + ' (' + mo.marker + ')'
mo = mo.sort_values('b_pi_2000')
fig, ax = plt.subplots(figsize=(11, 10.5))
yy = np.arange(len(mo))
ax.hlines(yy, mo.share_range_lo, mo.share_range_hi, color='#bbbbbb', lw=6, zorder=1)
for col, mk, cl, nm in [('a3000_zero_frac', 'v', '#009E73', '(a) zero share, cells >= 3,000 UMIs'),
                        ('a8000_zero_frac', '^', '#56B4E9', '(a) zero share, cells >= 8,000 UMIs'),
                        ('b_pi_2000', 'o', '#0072B2', '(b) NB mixture share at 2,000 UMIs'),
                        ('c_pi', 's', '#D55E00', '(c) Poisson mixture share')]:
    ax.scatter(mo[col], yy, marker=mk, s=70, c=cl, edgecolors='white', linewidths=0.8, zorder=3, label=nm)
ax.axvspan(0.15, 0.85, color='#eeeeee', zorder=0)
ax.text(0.5, len(mo) - 0.35, 'rule (b) band: 0.15 to 0.85', ha='center', va='bottom', fontsize=12, color=INK2)
ax.set_ylim(-0.7, len(mo) + 0.2)
ax.set_yticks(yy); ax.set_yticklabels(mo.lab, fontsize=12)
for t, r in zip(ax.get_yticklabels(), mo.itertuples()):
    if r.mosaic_verdict.startswith('mosaic'):
        t.set_fontweight('bold')
ax.set_xlabel('share of marker-negative cells'); ax.set_xlim(-0.02, 1.0)
ax.set_title('Marker-negative share by estimator (bold: passes the mosaic rule)', loc='left', fontsize=14)
ax.legend(frameon=False, loc='lower right', fontsize=11)
save(fig, 'fig5_mosaic_shares_by_estimator')
print('figures written to', FIG)
