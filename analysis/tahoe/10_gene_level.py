"""Descriptive: per-gene log2(CPM+1) of the Y genes, TXLNGY, PRKY and XIST per line, after the primary
foreign-cell filter (to read arm states, A5). No state is called here."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from common import OUT, Y8, log, load_cells

LOG = open(OUT + '/logs/10_gene_level.log', 'a')
genes = ['XIST'] + Y8 + ['TXLNGY', 'PRKY']
cells = load_cells(['drug', 'cell_line_id', 'total_umi', 'y8'] + genes)
S = pd.read_csv(OUT + '/partA_local_states.tsv', sep='\t').set_index('cellosaurus_id')
keep = np.ones(len(cells), bool)
yn = cells.cell_line_id.map(~S.y_high).to_numpy(bool)
xn = cells.cell_line_id.map(~S.xist_high).to_numpy(bool)
keep &= ~(yn & (cells.y8 >= 3).to_numpy()) & ~(xn & (cells.XIST >= 3).to_numpy())
c = cells[keep]
g = c.groupby('cell_line_id', observed=True)[['total_umi'] + genes].sum()
assert (g['total_umi'] == S.loc[g.index, 'pseudobulk_umi']).all()
out = pd.DataFrame(index=g.index)
for x in genes:
    out[x] = np.log2(1e6 * g[x] / g['total_umi'] + 1)
out.insert(0, 'cell_name', S.loc[out.index, 'cell_name'])
out.insert(1, 'annotation', S.loc[out.index, 'annotation'])
out.insert(2, 'state', S.loc[out.index, 'state'])
out.insert(3, 'arm_state', S.loc[out.index, 'arm_state'])
assert np.allclose(out[Y8].mean(axis=1), S.loc[out.index, 'y_score'], atol=1e-4)
out.reset_index().rename(columns={'cell_line_id': 'cellosaurus_id'}).sort_values('cell_name').to_csv(
    OUT + '/partA_local_gene_level_scores.tsv', sep='\t', index=False, float_format='%.4f')
log('gene-level log2(CPM+1) written', LOG)
print(out[out.annotation == 'Male'].sort_values('arm_state').round(2).to_string())
