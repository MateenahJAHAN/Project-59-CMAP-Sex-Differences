"""Recreate the four released Part B tables from the fresh model outputs."""
from pathlib import Path
from csv import DictReader, writer
from hashlib import sha256
import gzip, io, json
import argparse
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--work', type=Path, required=True)
HERE = parser.parse_args().work.resolve()
REPO = Path(__file__).resolve().parents[2]
OUT=HERE/'results/public_tables'
OUT.mkdir(exist_ok=True)
PUBLIC=REPO/'verification/results'
def records(file, delimiter='\t'):
    with file.open(encoding='utf-8', newline='') as stream:
        return list(DictReader(stream, delimiter=delimiter))

def table(file, header, rows):
    stream = io.StringIO(newline='')
    out = writer(stream, delimiter='\t', lineterminator='\n')
    out.writerow(header)
    out.writerows(rows)
    data = stream.getvalue().encode()
    if file.suffix == '.gz':
        data = gzip.compress(data, mtime=0)
    file.write_bytes(data)

def part_b():
    source=HERE/'results/G'
    expected=('G_H1.tsv','G_H3.tsv','G_H2_per_compound.tsv','G_gates.tsv','G_run_manifest.tsv')
    assert all((source/path).is_file() for path in expected), 'Fresh full Part B run is still writing results'
    terms={'factY':'Y','factX':'XIST'}
    def primary(rows):
        return [row for row in rows if row['panel']=='P1' and row['term'] in terms and row['analysis']=='primary']
    h1=primary(records(source/'G_H1.tsv'))
    assert len(h1)==2 and all(row['contrast_ok']=='True' for row in h1)
    h1_original=(
       'n_lines','n_g1','n_g0','n_lineages','B','n_tested','panel_stat','panel_p',
       'n_p_lt_005','n_q_lt_01','dmso_p','dmso_rank_if_inserted',
    )
    h1_public=(
       'n_lines','n_with_indicator','n_without_indicator','n_lineages','n_permutations',
       'n_tested_compounds','panel_statistic','panel_p','n_nominal_p_lt_0_05','n_bh_q_lt_0_10',
       'vehicle_p','vehicle_rank_if_inserted',
    )
    table(OUT/'part_b_h1.tsv',('term',*h1_public),
        [(terms[row['term']],*(row[field] for field in h1_original)) for row in h1])
    h3=[row for row in records(source/'G_H3.tsv') if row['panel']=='P1' and row['term'] in terms
          and row['role']=='primary' and row['scale']=='level3']
    assert len(h3)==2 and all(row['unit']=='percent of mean S' for row in h3)
    h3_original=('n_lines','df','estimate','lo90','hi90','lo95','hi95','labelfree_P_inside_margin')
    h3_public=('n_lines','residual_df','estimate_pct','lo90_pct','hi90_pct','lo95_pct','hi95_pct','probability_interval_inside_margin')
    table(OUT/'part_b_h3.tsv',('term',*h3_public),
        [(terms[row['term']],*(row[f] for f in h3_original)) for row in h3])
    h2=primary(records(source/'G_H2_per_compound.tsv'))
    p_fields=('compound','tested','in_ranking','p','q','finding')
    table(OUT/'part_b_compounds.tsv',('term','compound','tested','in_ranking','p','bh_q','finding'),
        sorted((terms[row['term']],*(row[f] for f in p_fields)) for row in h2))
    gates=[row for row in records(source/'G_gates.tsv') if row['gate'] in ('PC1','K4') and row['panel']=='P1']
    assert len(gates)==2 and all(row['passes']=='True' and row['compound']=='AMG-232' for row in gates)
    table(OUT/'part_b_controls.tsv',('control','response_scale','compound','n_lines','p','rank','n_ranked','passes'),
        [(row['gate'],row['scale'],row['compound'],row['n_lines'],row['p'],row['rank'],row['n_ranked'],row['passes']) for row in gates])
part_b()
summary=[]
for p in sorted(OUT.glob('*.tsv')):
    ref=PUBLIC/p.name
    a,b=sha256(p.read_bytes()).hexdigest(),sha256(ref.read_bytes()).hexdigest()
    summary.append(dict(file=p.name,rows=len(records(p)),candidate_sha256=a,published_sha256=b,bytes_equal=a==b))
(OUT/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
for row in summary:
    print(row['file'], row['rows'], 'rows;', 'exact' if row['bytes_equal'] else 'different')
if not all(row['bytes_equal'] for row in summary):
    raise SystemExit("Part B result differs from the released table. Inspect comparison.json.")

# Compare the complete analysis outputs to the independent earlier raw replay.
# Guard text, absolute paths and package provenance are checked separately.
evidence = REPO / 'verification/source_replay/part_b'
rows = []
for manifest, folder in [('model', HERE / 'results/G'), ('capture', HERE / 'capture'),
                         ('exploratory', HERE / 'results/G_exploratory')]:
    for row in records(evidence / (manifest + '_files.tsv')):
        name = row['file']
        if name in ('G_hash_check.tsv', 'G_run_manifest.tsv', 'G_decide_options.tsv'):
            continue
        artifact = folder / 'figures' / name if name.endswith('.png') else folder / name
        if name == 'G_H2_gene_level.tsv.gz':
            # The compressor includes the creation time. Compare every decompressed byte.
            observed = sha256(gzip.decompress(artifact.read_bytes())).hexdigest()
            expected = 'e16df5e0d3bf3ebd58a642a543977099315767d7aec6373dcd14afbb80cc5c42'
        else:
            observed = sha256(artifact.read_bytes()).hexdigest()
            expected = row['sha256_rerun']
        rows.append(dict(file=name, group=manifest, sha256=observed, expected=expected, matches=observed == expected))

manifest = {row['item']: row['value'] for row in records(HERE / 'results/G/G_run_manifest.tsv')}
assert manifest['label_status'] == 'real' and manifest['reduced_settings'] == ''
assert manifest['guard_execution'] == 'post-publication reproduction of registered analysis'
hashes = records(HERE / 'results/G/G_hash_check.tsv')
expected_hashes = {r['object']:r['expected_sha256'] for r in records(REPO / 'analysis/lincs/config/model_hash_checks.tsv')}
assert {r['object']:r['expected_sha256'] for r in hashes} == expected_hashes
assert len(hashes) == len(expected_hashes)
assert all(r['status'] in ('equal', 'equal (alternative)') for r in hashes)
decisions = {row['decide_id']: row['option'] for row in records(HERE / 'results/G/G_decide_options.tsv')}
cfg = json.loads((REPO / 'analysis/lincs/config/G_defaults.json').read_text())['decide']
assert decisions == {k: v if isinstance(v, str) else f"unit={v['unit']};interval={v['interval']}" for k, v in cfg.items()}
for name in ('figure4_intervals.tsv', 'figure4_nulls.tsv.gz'):
    new, ref = HERE / 'results/figure4' / name, evidence / name
    payload = (lambda p: gzip.decompress(p.read_bytes())) if name.endswith('.gz') else (lambda p: p.read_bytes())
    a, b = sha256(payload(new)).hexdigest(), sha256(payload(ref)).hexdigest()
    rows.append(dict(file=name, group='figure4', sha256=a, expected=b, matches=a == b))
for suffix in ('pdf', 'png'):
    name = 'fig4_partB.' + suffix
    observed = sha256((HERE / 'results/figure4' / name).read_bytes()).hexdigest()
    expected = sha256((REPO / 'manuscript/figures' / name).read_bytes()).hexdigest()
    rows.append(dict(file=name, group='figure4', sha256=observed, expected=expected, matches=observed == expected))
(OUT / 'full_comparison.json').write_text(json.dumps(rows, indent=2) + '\n')
failed = [r['file'] for r in rows if not r['matches']]
if failed:
    raise SystemExit(f"{len(failed)} files differ from the earlier source replay. Inspect full_comparison.json: {', '.join(failed)}")
print(f'{len(rows)} scientific result, capture and figure-source files agree with the independent replay.')
