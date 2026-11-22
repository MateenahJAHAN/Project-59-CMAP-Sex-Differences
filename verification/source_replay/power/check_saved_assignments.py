"""Compare saved power replicates with the archive before the full run has finished."""
import argparse
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd


def compare_rows(records, archive_path, ids, keys):
    fresh = pd.DataFrame(records)
    # The original driver writes missing numeric and text values with the same TSV conventions.
    fresh = pd.read_csv(io.StringIO(fresh.to_csv(sep='\t', index=False)), sep='\t',
                        float_precision='round_trip')
    old = pd.read_csv(archive_path, sep='\t', float_precision='round_trip')
    old = old[old.rep.isin(ids)]
    fields = [name for name in old.columns if name != 'time_s']
    if set(fresh.columns) != set(old.columns):
        raise ValueError(f'Unexpected columns in {archive_path.name}')
    a = fresh[fields].set_index(keys).sort_index()
    b = old[fields].set_index(keys).sort_index()
    if not a.index.is_unique or not b.index.is_unique or not a.index.equals(b.index):
        raise ValueError(f'Mismatched scientific record keys in {archive_path.name}')
    diff = a.compare(b)
    return {'rows': len(a), 'scientific_fields': len(fields), 'changed_rows': len(diff),
            'changed_fields': sorted(set(diff.columns.get_level_values(0)))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('checkpoints', type=Path)
    parser.add_argument('archive', type=Path)
    parser.add_argument('report', type=Path)
    args = parser.parse_args()
    files = sorted(args.checkpoints.glob('rep_*.npz'))
    if not files:
        raise ValueError('No saved assignments')
    rows, plants, ids = [], [], []
    bad_intervals = []
    entries = 0
    with np.load(args.archive/'power_e2.npz', allow_pickle=False) as data:
        intervals = data['E2']
    for file in files:
        index = int(file.stem.split('_')[-1])
        with np.load(file, allow_pickle=False) as data:
            calls, planted = json.loads(data['records'].item())
            assert len(calls) == 26 and len(planted) == 120
            assert all(row['rep'] == index for row in calls+planted)
            values = data['intervals']
            assert values.shape == intervals[index].shape
            if not np.array_equal(values, intervals[index], equal_nan=True):
                bad_intervals.append(index)
            entries += values.size
        ids.append(index)
        rows.extend(calls)
        plants.extend(planted)
    evidence = {
        'checkpoint_group': args.checkpoints.name,
        'planned_replicates': int(intervals.shape[0]),
        'completed_replicates_checked': len(ids),
        'replicate_ids': ids,
        'scenario_comparison': compare_rows(rows, args.archive/'power_runs.tsv.gz', ids,
                                            ['rep', 'term', 'scenario']),
        'planted_comparison': compare_rows(plants, args.archive/'power_planted.tsv.gz', ids,
                                          ['rep', 'term', 'scenario', 'compound']),
        'interval_entries_compared': int(entries),
        'replicates_with_interval_differences': bad_intervals,
    }
    args.report.write_text(json.dumps(evidence, indent=2)+'\n')
    print({k: v for k, v in evidence.items() if k != 'replicate_ids'})


if __name__ == '__main__':
    main()
