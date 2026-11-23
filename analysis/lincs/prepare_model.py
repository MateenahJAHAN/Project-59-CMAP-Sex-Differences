"""Assemble the registered design and verify all reconstructed LINCS response tensors."""
import csv
import gzip
import hashlib
from pathlib import Path
import sys

import numpy as np
import pandas as pd

from reproduce import HERE, file_hash, link

work = Path(sys.argv[1]).resolve()
results = work / 'results'
design = HERE / 'design'


def table(path, **options):
    return pd.read_csv(path, sep='\t', **options)


def write_tp53_lines():
    annotations = table(results / 'E/E_line_annotations.tsv', dtype=str).set_index('cell_iname')
    lines = (work / 'data/E/tensors/P1_lines.txt').read_text().splitlines()
    a = annotations.loc[lines].copy()
    known = a[a.tp53_status.isin(['mutant', 'wild_type'])].copy()
    both = known.groupby('cell_lineage').tp53_status.nunique()
    use = known[known.cell_lineage.isin(both[both == 2].index)].copy()
    use['tp53_mutant'] = (use.tp53_status == 'mutant').astype(int)
    use[['cell_lineage', 'tp53_mutant']].to_csv(results / 'E/E_pc1_lines.tsv', sep='\t')


def check_activity():
    # Design strata use within-panel correlations between cell lines. D fixed the
    # split before the model. Recompute all scores before reusing that fixed split.
    scores = []
    for panel in ('P1', 'P1b', 'P2'):
        tensor = np.load(work / 'data/D/tensors' / (panel + '_resp.npy'))
        drugs = table(results / 'D/tensors' / (panel + '_compounds.tsv')).cmap_name
        for j, drug in enumerate(drugs):
            y = tensor[~np.isnan(tensor[:, j, 0]), j].astype(np.float64)
            score = np.nan
            if len(y) > 1:
                score = np.corrcoef(y)[np.triu_indices(len(y), 1)].mean()
            scores.append((panel, drug, score, len(y)))
    current = pd.DataFrame(scores, columns=['panel', 'cmap_name', 'activity_mean_pairwise_r', 'n_lines'])
    current['panel_median'] = current.groupby('panel').activity_mean_pairwise_r.transform('median')
    current['above_panel_median'] = current.activity_mean_pairwise_r > current.panel_median
    prior = table(design / 'taskD/D_activity_scores.tsv', float_precision='round_trip')
    pd.testing.assert_frame_equal(prior, current, check_exact=False, atol=1e-12, rtol=0)
    current.to_csv(results / 'D/D_activity_scores_rebuilt.tsv', sep='\t', index=False)


def check_tensors():
    frozen = table(design / 'taskE/E_hashes.tsv', dtype=str).set_index('object')
    actual = table(results / 'E/E_hashes.tsv', dtype=str).set_index('object')
    if not actual.index.is_unique:
        raise ValueError('Duplicate well-array identities')
    checks = []
    for key, src in frozen.iterrows():
        path = work / src.local_path
        if src.sha256 != actual.at[key, 'sha256'] or file_hash(path) != src.sha256:
            raise ValueError(f'Well-array checksum differs: {key}')
        checks.append((str(path.relative_to(work)), src.sha256))
    expected = table(design / 'taskD/D_hashes.tsv').set_index('path').sha256
    for row in table(results / 'D/D_tensor_hashes.tsv').itertuples():
        if row.path not in expected.index or expected[row.path] != row.sha256:
            raise ValueError(f'Level-5 tensor manifest differs from the archived fingerprint: {row.path}')
        if file_hash(work / row.path) != row.sha256:
            raise ValueError(f'Level-5 panel checksum differs: {row.path}')
        checks.append((row.path, row.sha256))
    pd.DataFrame(checks, columns=['file', 'sha256']).to_csv(results / 'tensor_checks.tsv', sep='\t', index=False)


with (HERE / 'config/pseudolabel_checks.tsv').open() as bank_digests:
    for bank in csv.DictReader(bank_digests, delimiter='\t'):
        path = Path(bank['file'])
        generated = results / path.parts[0][4:] / Path(*path.parts[1:])
        with gzip.open(generated, 'rb') as stream:
            checksum = hashlib.file_digest(stream, 'sha256').hexdigest()
        if checksum != bank['uncompressed_sha256']:
            raise ValueError(f'Pseudo-label assignments differ: {path}')
write_tp53_lines()
check_activity()
check_tensors()
with (HERE / 'config/G_inputs_pinned.tsv').open() as stream:
    for row in csv.DictReader(stream, delimiter='\t'):
        name = Path(row['input'])
        group, rest = name.parts[0][4:], Path(*name.parts[1:])
        source = design / name
        if not source.exists():
            source = results / group / rest
        if source.is_file():
            link(source, work / 'model_inputs' / name)
# Exploratory scripts use the same historical layout in their work directory.
for src, dst in [
    (work / 'model_inputs', 'data/G/inputs'), (work / 'data/D/tensors', 'data/G/tensors/D'),
    (work / 'data/E/tensors', 'data/G/tensors/E'), (work / 'sources/siginfo_beta.txt', 'data/G/meta/siginfo_beta.txt')
]:
    link(src, work / dst)
(work / 'provenance').mkdir(exist_ok=True)
