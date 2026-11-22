"""Compare released and source-rebuilt figures, including all plotted fields and PDF pages."""
from collections import Counter
import argparse
import csv
import hashlib
import json
from pathlib import Path
import fitz


FIGURES = ('fig1_design', 'fig2_label_vs_pattern', 'fig3_agreement', 'fig4_partB',
           'figS1_ambient', 'figS2_a549', 'figS3_power')
PROVENANCE = {'source', 'source_file', 'row_source', 'col_source', 'source_line', 'source_column'}


def table(path):
    with path.open(newline='') as handle:
        rows = csv.DictReader(handle, delimiter='\t')
        return rows.fieldnames, list(rows)


def compare_one(name, released, old_tables, rerun):
    base = released / (name+'.pdf')
    new = rerun / (name+'.pdf')
    with fitz.open(base) as a, fitz.open(new) as b:
        assert len(a) == len(b) == 1, name
        pa, pb = a[0], b[0]
        same_size = pa.rect == pb.rect
        same_text = pa.get_text().split() == pb.get_text().split()
        same_pixels = same_size and pa.get_pixmap(dpi=144).samples == pb.get_pixmap(dpi=144).samples
    original_pdf_sha = hashlib.sha256(base.read_bytes()).hexdigest()
    rebuilt_pdf_sha = hashlib.sha256(new.read_bytes()).hexdigest()
    old_cols, old_rows = table(old_tables / (name+'_plotted.tsv'))
    new_cols, new_rows = table(rerun / (name+'_plotted.tsv'))
    if new_cols != old_cols:
        raise ValueError(f'{name}: plotting table columns differ')
    scientific_diffs, provenance_diffs = [], set()
    for i, (old, current) in enumerate(zip(old_rows,new_rows), start=2):
        for col in old_cols:
            if old[col] == current[col]:
                continue
            if col in PROVENANCE:
                provenance_diffs.add(i)
            else:
                scientific_diffs.append({'table_line': i, 'field': col,
                                         'archived': old[col], 'recomputed': current[col]})
    scientific_fields = [col for col in old_cols if col not in PROVENANCE]
    def counts(records):
        return Counter(tuple(record[col] for col in scientific_fields) for record in records)
    old_counts, new_counts = counts(old_rows), counts(new_rows)
    def changed_records(counter):
        return [{'count': count, 'record': dict(zip(scientific_fields, key))}
                for key, count in counter.items()]
    removed = changed_records(old_counts - new_counts)
    added = changed_records(new_counts - old_counts)
    return {'figure': name, 'released_pdf_sha256': original_pdf_sha, 'rebuilt_pdf_sha256': rebuilt_pdf_sha,
            'identical_pdf': original_pdf_sha == rebuilt_pdf_sha, 'identical_page_size': same_size,
            'identical_pdf_text': same_text, 'identical_pixels_144dpi': same_pixels,
            'archived_plotted_rows': len(old_rows), 'rebuilt_plotted_rows': len(new_rows),
            'rows_with_changed_source_fields': len(provenance_diffs),
            'positional_field_differences': len(scientific_diffs),
            'positional_row_differences': len({change['table_line'] for change in scientific_diffs}),
            'scientific_records_removed': removed, 'scientific_records_added': added}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--released-figures', type=Path, required=True)
    parser.add_argument('--reference-tables', type=Path, required=True)
    parser.add_argument('--rebuilt', type=Path, required=True)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--figures', choices=FIGURES, nargs='+', default=FIGURES)
    args = parser.parse_args()
    records = [compare_one(name, args.released_figures, args.reference_tables, args.rebuilt)
               for name in args.figures]
    args.report.write_text(json.dumps({'figures': records}, indent=2)+'\n')
    for row in records:
        print(row['figure'], row['rebuilt_plotted_rows'], 'rows,', len(row['scientific_records_removed']), 'removed,',
              len(row['scientific_records_added']), 'added scientific records, PDF byte match:', row['identical_pdf'])


if __name__ == '__main__':
    main()
