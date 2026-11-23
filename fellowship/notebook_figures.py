"""Render the Phase 1, 2 and 6 notebook plots using newly computed results."""
import argparse
from contextlib import chdir, redirect_stdout
import io
import json
from pathlib import Path

from run import outside_checkout

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]


def notebook(phase):
    path = next(REPO.glob(f'Phase{phase}*/**/*.ipynb'))
    return [''.join(cell['source']) for cell in json.loads(path.read_text())['cells']]


def plot_cell(source, variables):
    # Archived plot code writes relative to the fresh per-phase output directory.
    with redirect_stdout(io.StringIO()):
        exec(source.replace('plt.show()', ''), variables)
    plt.close('all')


def pilot(out, data, cells):
    tests = pd.read_csv(out/'phase1/gene_tests.csv.gz', dtype={'gene_id':str})
    annot = pd.read_csv(data/'geneinfo_beta.txt', sep='\t', dtype={'gene_id':str})
    tests = tests.merge(annot[['gene_id','gene_symbol']], on='gene_id', validate='one_to_one')
    strong = tests[(tests.p_adj<.05) & (tests.difference.abs()>1)]
    variables = dict(res_named=tests, strong_named=strong)
    for n in (28,29):
        plot_cell(cells[n], variables)


def fisher(out, data, cells):
    tests = pd.read_csv(out/'comparison/phase2_all_tissue_tests.csv')
    compounds = pd.read_csv(out/'comparison/phase2_compounds.csv').rename(columns={'q_original_call':'category'})
    vorinostat = tests.loc[tests.drug=='vorinostat'].copy()
    # The notebook rounded fold enrichment before drawing this chart.
    vorinostat['fold'] = vorinostat['fold'].round(2)
    variables = dict(results_df=vorinostat, df=compounds)
    for n in (49,73):
        plot_cell(cells[n], variables)


def unfiltered(out, data, cells):
    base = Path.cwd()
    for folder in ['results','figures']:
        (base/folder).mkdir(exist_ok=True)
    signatures = pd.read_parquet(out/'prepared/signatures.parquet')
    tissues = pd.read_csv(out/'phase6/tissue_results.csv.gz')
    best = tissues.groupby(['drug','lineage']).category.min()
    links = signatures.join(best, on=['cmap_name','cell_lineage'])
    counts = links.groupby(['cmap_name',links.category.fillna('D')]).size().unstack(fill_value=0)
    counts['total'] = counts.sum(axis=1)
    counts.to_csv(base/'results/signature_ABCD_per_drug.csv')
    variables = dict(BASE=str(base), sig_known=signatures, gctx_sigs=set(signatures.sig_id),
                     top50_names=signatures.cmap_name.unique(), np=np, pd=pd, plt=plt)
    start = cells[43].index('S = sig_known[')
    for source in [cells[43][start:], cells[44], cells[48], cells[50]]:
        plot_cell(source, variables)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=outside_checkout, required=True)
    parser.add_argument('--out', type=outside_checkout, required=True, help='Existing fellowship run directory')
    parser.add_argument('--phases', type=int, nargs='+', choices=(1,2,6), default=[1,2,6])
    args=parser.parse_args()
    out=args.out.resolve()
    data=args.data.resolve()
    for n in args.phases:
        dest = out/'notebook_figures'/f'phase{n}'
        dest.mkdir(parents=True,exist_ok=True)
        with chdir(dest):
            {1:pilot,2:fisher,6:unfiltered}[n](out,data,notebook(n))
        print(f'Phase {n}: {sum(1 for _ in dest.rglob("*.png"))} notebook figures in {dest}')


if __name__=='__main__':
    main()
