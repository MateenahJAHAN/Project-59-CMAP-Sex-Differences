"""Draw manuscript Figure 4 using the local model, independently refitted values and real panel doses."""
from pathlib import Path
import json
import sys

import pandas as pd

from reproduce import HERE, REPOSITORY, link

sys.path.insert(0, str(REPOSITORY / 'manuscript/plotting'))
import fig4_partB
import p59style

work = Path(sys.argv[1]).resolve()
model = work / 'results/G'
out = work / 'results/figure4'
inputs = out / 'plot_inputs/part_b'
inputs.mkdir(parents=True, exist_ok=True)


def read(path):
    return pd.read_csv(path, sep='\t', float_precision='round_trip')


summary = json.loads((out / 'figure4_nulls_check.json').read_text())
h1 = read(model / 'G_H1.tsv')
controls = read(model / 'G_controls.tsv')
sets = read(model / 'G_sets_counts.tsv')
groups = sets[(sets.panel == 'P1') & (sets.family == 'factorial')]
points = []
for term, letter, *_ in fig4_partB.NULLS:
    if term == 'TP53':
        select = controls[(controls.panel == 'P1') & (controls.exposure == fig4_partB.TP53)]
        n0 = int(select.n_g0.item())
    else:
        select = h1[(h1.panel == 'P1') & (h1.term == term) & (h1.analysis == 'primary')]
        n0 = int(groups.n_g3.sum())  # the common reference group has neither marker
    row = select.squeeze()
    result = summary[term]
    assert row.panel_stat == result['observed'] and row.panel_p == result['panel_p']
    values = {'observed':row.panel_stat, 'p':row.panel_p, 'n_at_or_above':result['at_or_above'],
              'n_g1':int(row.n_g1), 'n_g0':n0, 'n_tested':int(row.n_tested), 'n_q_lt_01':int(row.n_q_lt_01)}
    points.extend(dict(panel=letter, element=k, value=v) for k, v in values.items())

doses = read(work / 'results/D/D_panel_doses.tsv')
panels = read(work / 'results/D/D_panel_summary.tsv').set_index('panel')
for panel, frame in doses.groupby('panel'):
    modal_dose = frame.loc[frame.n_signatures.idxmax(), 'pert_dose']
    dose = f'{modal_dose:g}' if modal_dose >= 10 else f'{modal_dose:.1f}'
    points.append(dict(panel='D', element='row_label', series=panel, name='head',
                       text=f'{dose} µM, {panels.at[panel, "pert_itime"]}'))
config = json.loads((HERE / 'config/G_defaults.json').read_text())
points.append(dict(panel='D', element='margin', value=config['plan_constants']['E3_margin_pct']))
pd.DataFrame(points).to_csv(inputs / 'fig4_partB_plotted.tsv', sep='\t', index=False)
for name in ('figure4_nulls.tsv.gz', 'figure4_intervals.tsv'):
    link(out / name, inputs / name)
figure, _ = fig4_partB.build(inputs.parent)
p59style.save(figure, str(out / 'fig4_partB'))
