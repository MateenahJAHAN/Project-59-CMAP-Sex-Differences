"""Verify the entire aggregate package using the Python standard library."""
from pathlib import Path
import csv, hashlib, json
base=Path(__file__).resolve().parents[1]
def read(path):
    with (base/path).open(newline='') as file:return list(csv.DictReader(file,delimiter='\t'))
def receipt(name):return json.loads((base/'receipts'/name).read_text())
for row in read('SHA256SUMS.tsv'):
    path=base/row['file']
    assert hashlib.sha256(path.read_bytes()).hexdigest()==row['sha256'],path
sources=receipt('source_verification.json')
assert len(sources['sciplex_GEO_files'])==16 and sources['emerald_exact_file_matches']==116
src=receipt('sciplex_source_coverage.json')
assert sum(x['matrix_rows'] for x in src.values())==1171291549
assert sum(x['n_vehicle'] for x in src.values())==23570
obs=receipt('tahoe_observations.json')
assert (obs['rows'],obs['full_vehicle_rows'],obs['missing_or_different'])==(100648790,2205786,0)
em=read('receipts/emerald_vehicle_coverage.tsv')
assert len(em)==116 and sum(int(r['vehicle_rows']) for r in em)==130498
lines=read('tables/F1_emerald_states.tsv')
assert len(lines)==52 and sum(int(r['n_vehicle_cells']) for r in lines)==100253
assert len(read('tables/F1_emerald_units.tsv'))==207
for fn, pos, tot in [('str_fresh_EmeraldBay.tsv',47,52),('str_fresh_sci-Plex.tsv',3,3)]:
    rr=read('tables/'+fn)
    assert sum(r['match']=='True' for r in rr)==pos and sum(r['comparable']=='True' for r in rr)==tot
assert len(read('receipts/cellosaurus_record_checks.tsv'))==1102
z=read('tables/sciplex_vs_fresh_bulk.tsv')
assert 'bulk_y' not in z[0] and 'bulk_xist' not in z[0]
assert len(z)==6 and [r['experiment']+'/'+r['line'] for r in z if r['agree']=='False']==['sciPlex3/A549']
assert all(v['new_sha256']==v['old_sha256'] for v in receipt('compare_amelogenin.json').values())
for name,expected in [('CCLE',(773,734,714)),('CMP',(1082,1031,1007))]:
    rr=[r for r in read('tables/str_fresh_bulk.tsv') if r['source']==name]
    assert (len(rr),sum(r['comparable']=='True' for r in rr),sum(r['agree']=='True' for r in rr))==expected
    assert all(r['call_changed']=='False' for r in rr)
primary=receipt('primary_caller_comparison.json')
assert len(primary)==33 and all(v['new_sha256']==v['old_sha256'] for v in primary.values())
agg=receipt('tahoe_aggregate_comparison.json')
assert agg['partial_tables']==7640 and not agg['different_partial_tables'] and not agg['per_shard_differences']
assert len(agg['aggregate_outputs'])==21 and all(z['byte_identical'] for z in agg['aggregate_outputs'].values())
raw=receipt('tahoe_aggregate_checks.json')['checks']
assert raw['cells']==2205786 and raw['line_plate_pairs_equal_obs']==raw['line_plate_pairs']==700
labels=read('tables/partA_local_label_audit_summary.tsv')
assert any(r['n_labelled']=='49' and r['n_without_usual_pattern']=='27' for r in labels)
bulk=read('tables/patterns_tahoe_bulk.tsv')
assert len(bulk)==49 and sum(r['agree']=='True' for r in bulk)==47
cross=read('tables/patterns_cross_atlas.tsv')
assert len(cross)==47 and sum(r['match']=='True' for r in cross)==42
assert sum(r['eligible']=='True' for r in cross)==41 and sum(r['eligible']==r['match']=='True' for r in cross)==40
d=receipt('str_figure_inputs.json')
assert sum(v['fields_compared'] for v in d['tables'].values())==14351
assert not any(v['changed_fields'] for v in d['tables'].values())
for r in read('tables/figure3_agreement_by_source.tsv'):
    dataset=r['source'];rows=read('tables/figure3_supp_table_amelogenin_by_line.tsv')
    counts={k:sum(x[dataset+'_vs_STR']==k for x in rows) for k in ['agree','type A','type B']}
    assert int(r['n_compared'])==sum(counts.values()) and int(r['n_agree'])==counts['agree']
status=receipt('run_status.json')
if 'all_EmeraldBay_expression_shards_and_identity_controls' in status['completed']:
    coverage=read('receipts/emerald_all_condition_coverage.tsv')
    assert len(coverage)==116 and sum(int(r['rows']) for r in coverage)==1831648
    pairs=read('tables/emerald_line_pair_correlations.tsv')
    assert len(pairs)==1326 and len({(r['cellosaurus_id_a'],r['cellosaurus_id_b']) for r in pairs})==1326
    differences=receipt('external_difference_details.json')
    assert differences['per_identifier_nearest_lines_changed']==0 and differences['per_identifier_counts_changed']==0
if 'two_caller_consensus' in status['completed']:
    cc=read('tables/tahoe_partA_states_consensus.tsv')
    assert len(cc)==50 and len({r['cellosaurus_id'] for r in cc})==50
    k9=receipt('k9_summary.json')
    assert sum(r['state_loc']==r['state_f2'] for r in cc)==int(k9['k9_primary_same_state'])
    assert sum(r['state_loc']==r['state_f2'] and r['confidence_loc']==r['confidence_f2'] for r in cc)==int(k9['k9_strict_same_state_and_confidence'])
    assert len(read('tables/A3_pairs.tsv'))==28 and len(read('tables/A3_replicates.tsv'))==70
    assert len(read('receipts/F2_filelist.tsv'))==3395 and len(read('receipts/F2_shard_summary.tsv'))==3388
    ix=read('receipts/F2_dmso_index.tsv')
    assert len(ix)==4048 and sum(int(r['vehicle_rows']) for r in ix)==2205786
    assert len(read('receipts/F2_index_vs_obs.tsv'))==700
if 'all_22_fits_with_manifest_pinned_numpy_pandas_and_scipy' in status['completed']:
    versions=receipt('mosaic_version_comparison.json')
    assert versions['table_s1_rows']==1694
    assert versions['source_vs_pinned_different_fields']==versions['source_vs_main_different_fields']==2
    assert versions['pinned_environment']['numpy']=='2.5.3' and versions['main_environment']['numpy']=='2.5.4'
    assert versions['pinned_environment']['pandas']=='3.0.5' and versions['main_environment']['pandas']=='3.1.1'
    assert versions['pinned_environment']['scipy']==versions['main_environment']['scipy']=='1.18.1'
    bounds=read('tables/table_S1_mosaic_sensitivity.tsv')
    assert [(r['accession'],r['archived'],r['main'],r['pinned']) for r in bounds]==[('CVCL_0399','0.008423','0.008421','0.008420'),('CVCL_1724','0.259381','0.259388','0.259378')]
    primary={(r['cell_name'],r['marker']):r['mosaic'] for r in read('tables/F2_mosaic.tsv')}
    pinned={(r['cell_name'],r['marker']):r['mosaic'] for r in read('tables/F2_pinned_mosaic.tsv')}
    assert primary==pinned and len(pinned)==22
    assert receipt('mosaic_numeric_summary.json')['consensus_verdicts_unchanged']==22
    digests=receipt('mosaic_scientific_sha256.json')['scientific_output_sha256']
    assert digests['supplementary_table_S1.tsv']==versions['S1_sha256']
    for ref,runtime,name in [('F2_mosaic.tsv','main','F2_mosaic.tsv'),('F2_mosaic.tsv','pinned','F2_pinned_mosaic.tsv'),('mosaic_resolution.tsv','main','mosaic_resolution.tsv'),('mosaic_resolution.tsv','pinned','mosaic_resolution_pinned.tsv')]:
        assert hashlib.sha256((base/'tables'/name).read_bytes()).hexdigest()==digests[ref][runtime]
    assert all(z['pixel_differences_vs_fresh']['original_tables_current_renderer']==0 for nm in ['F1','F2'] for z in receipt(nm+'_rendering_comparison.json').values())
print('All file checks and cohort, call, and STR counts passed.')
