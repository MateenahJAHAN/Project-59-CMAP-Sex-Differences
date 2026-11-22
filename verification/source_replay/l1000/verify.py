"""Check a local rerun against the published aggregate results and fingerprint receipts."""
import argparse
import hashlib
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'results')
    parser.add_argument('--reference',type=Path,default=HERE/'reference')
    parser.add_argument('--states',type=Path,help='Optional local Task A table to match to this run receipt')
    args=parser.parse_args()
    out=json.loads((args.output/'receipt.json').read_text())
    published=json.loads((HERE/'results/receipt.json').read_text())
    original=json.loads((args.reference/'reference_fingerprints.json').read_text())
    assert out['verification']['rounded_marker_medians_sha256']==original['rounded_marker_medians_sha256']
    assert out['sampling']['cohort_sha256']==original['cohort_sha256']
    assert out['source']['shape']==[188708,12328]
    assert out['sampling']['cohort_rows']==2800 and out['sampling']['cohort_lines']==127
    assert out['sampling']['qc_pass_DMSO_in_matrix']==75705
    assert out['verification']['selected_raw_bytes']==276147200
    assert out['verification']['selected_well_rows_sha256']==published['verification']['selected_well_rows_sha256']
    assert out['sampling']['archived_cohort_matches_every_field']
    assert out['verification']['all_marker_medians_match_archive_to_four_decimals']
    if args.states:
        assert sha(args.states)==out['verification']['input_sha256'][args.states.name]
    assert len(out['markers'])==len({r['matrix_index'] for r in out['markers'].values()})==9
    metrics=out['results']
    assert metrics['spearman_n']==124
    assert math.isclose(metrics['spearman_rho'],0.21317656242837743)
    for name,pos,neg,expected in [('original_AUC_XY_vs_Y_lost',36,19,0.533625730994152),
                                  ('alternate_AUC_all_Y_high_vs_all_Y_low',37,90,0.6510510510510511)]:
        record=metrics[name]
        assert (record['positive_lines'],record['negative_lines'])==(pos,neg)
        assert math.isclose(record['auroc'],expected)
    print('Verified full-cohort fingerprints, source dimensions, marker map and aggregate statistics.')


if __name__=='__main__':
    main()
