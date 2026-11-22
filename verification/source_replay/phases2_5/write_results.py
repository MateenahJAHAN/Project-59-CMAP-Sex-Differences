"""Write the completed 50-compound results from the per-compound comparisons."""
import csv
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

def read_json(file):
    return json.loads((HERE / file).read_text())

def read_csv(file):
    with (HERE / file).open() as f:
        return list(csv.DictReader(f))

three = read_json('comparison/phase3.json')
four = read_json('comparison/phase4.json')
five = read_json('sensitivity_gseapy_1.3.1/comparison/phase5.json')
compare_pairs = read_csv('comparison/phase4_pair_counts.csv')
drugs = {}
for row in compare_pairs:
    group = drugs.setdefault(row['drug'], [set(),set()])
    for i, key in enumerate(('', 'replay_')):
        for label in ['A','B']:
            if float(row[key+label]) > 0:
                group[i].add(label)
with (HERE / 'comparison/phase4_compounds.csv').open('w',newline='') as f:
    writer = csv.writer(f,lineterminator='\n')
    writer.writerow(['compound','archive_best_category','replay_best_category'])
    for name,(a,b) in sorted(drugs.items()):
        writer.writerow([name,min(a,default='C'),min(b,default='C')])
changed4 = [r['compound'] for r in read_csv('comparison/phase4_compounds.csv') if r['archive_best_category'] != r['replay_best_category']]
dp = read_csv('phase3_dl_pdmp/target_tissue_scores.csv')
ft = read_csv('phase3_literal/home_tissue_scores.csv')
ft = [r for r in ft if r['home_tissue'] == 'Cells_EBV.transformed_lymphocytes' and r['lineage'] == 'haematopoietic_and_lymphoid_tissue']
assert len(ft) == len(dp) == 4
foster = ', '.join(f"{float(r['FDR']):.4f}" for r in ft if r['scheme'].startswith('tied_order_'))
dp_ties = ', '.join(f"{float(r['FDR q-val']):.4f}" for r in dp if r['scheme'].startswith('tied_order_'))
a4,b4 = four['original_categories'], four['replay_categories']
retained4 = a4['A']+a4['B']-four['archive_only']
text = f"""# CMap and GTEx fellowship screens

The source-data replay covers all 50 compounds in each of the four phases. The analyses used 36,523 LINCS2020 Level 5 signatures and GTEx v8 tissue gene sets. All 450,255,544 expression values matched the provider matrix byte for byte. The GTEx gene-set source also matched. [LINCS comparison](input_verification/comparison.json) · [GTEx comparison](input_verification/gtex_comparison.json)

| Phase | What reproduced | What changed |
| --- | --- | --- |
| 2. Gene-list overlap | All 50 compound categories, numbers of differentially expressed genes and signature counts across 616,400 gene tests and 2,200 tissue comparisons. | AS-605240 had two genes overlapping liver, DUSP4 and STX1A. Its adjusted *p* was 0.02471 across the 20 tissues with an overlap and 0.05437 across all 44 eligible tissues. Including all 44 left zero compounds below 0.05. [Tissue results](comparison/phase2_all_tissue_tests.csv) |
| 3. Within-lineage rankings | {three['pairs']} eligible compound–lineage pairs produced {three['tissue_results']:,} tests. The pair counts matched for all 50 compounds. | The archived own-tissue, other-tissue-only and no-hit counts were 3, 23 and 24. The replay gave 4, 23 and 23. Of the 23 no-hit compounds, {three['tested_compounds_without_a_hit']} underwent testing and {three['compounds_with_no_eligible_lineage']} had no eligible lineage. The number of enrichments changed for {three['compounds_with_different_evidence_counts']} compounds. [Compound results](comparison/phase3_compounds.csv) |
| 4. Split reference sets | {four['pairs']} eligible pairs produced {four['all_tested_sets']:,} tests. Both archived CSVs contained the same 9,260 rows. All 603 counts in the saved notebook output agree with both exports. | The archive contained {a4['A']:,} own-tissue hits and {a4['B']:,} other-tissue hits. The replay gave {b4['A']:,} and {b4['B']:,}. Of the 9,179 archived hits, {retained4:,} reproduced and {four['archive_only']} crossed above or reached FDR 0.05. Another {four['replay_only']} crossed below 0.05. [Changed calls](comparison/phase4_changed_threshold_calls.csv) · [Export checks](comparison/phase4_historical_exports.json) |
| 5. Pooled rankings after filtering | All {five['rows']:,} categories, skip reasons, donor-group counts and home-tissue assignments matched. All 620 tested normalized enrichment scores (NES) and false discovery rates (FDR) matched at saved precision with GSEApy 1.3.1. | GSEApy 1.2.1 changed 14 FDR values by at most 0.001 and preserved every classification. The notebook did not record its installed version. [Software comparison](comparison/phase5_software_versions.csv) · [Workbook checks](comparison/auxiliary_exports.json) |

## The two compounds that changed Phase 3 category

- **Fostamatinib** gained an own-tissue hit. The blood-derived cell lines enriched the GTEx EBV-transformed lymphocyte set (NES 1.18534, FDR 0.0410). Of 12,276 ranked genes, 12,002 shared their score with another gene. Reordering only ties gave FDR values of {foster}. The archived file stores two other-tissue scores for this drug. [Full notebook replay and comparison](phase3_literal/archived_and_literal_scores.csv) · [Tie order results](phase3_literal/home_tissue_scores.csv)
- **DL-PDMP** moved from no hit to other-tissue only. Large-intestine cell lines enriched the same EBV-transformed lymphocyte set (NES 1.32482, FDR 0.0355). Reordering only ties gave FDR values of {dp_ties}. The archive contains the no-hit classification and no individual tissue scores for this drug. [Full notebook replay and comparison](phase3_dl_pdmp/archived_and_literal_scores.csv) · [Tie order results](phase3_dl_pdmp/target_tissue_scores.csv)

The original notebook functions and original cmapPy reader reproduced these two compounds’ fresh results. Each tie-order test preserved the gene scores, seed 42 and 1,000 GSEA permutations. The archive retains {three['archived_individually_scored_tissue_results']} individually scored records out of {three['original_total_evidence']} archived hits. Absent scores for unreported tissues cannot be numerically reconstructed from that export. The full replay yielded {three['replay_total_evidence']} hits. [Archived scores paired with the replay](comparison/phase3_archived_leading_evidence.csv) · [All replayed hits](comparison/phase3_replayed_significant_tissues.csv)

## The Phase 4 permutation check

For NVP-BEZ235, the placenta comparison used 12 signatures from one female-labelled line. The female-biased GTEx minor salivary gland set had archived NES −1.747 and FDR 0.005. The notebook replay gave NES −1.000 and FDR 0.503. The observed enrichment score stayed at −0.280708 under three other permutation seeds, each of which gave FDR 0.005. Increasing the permutation count from 200 to 1,000 and 2,000 with the original seed gave FDRs of 0.001 and 0.0005. This demonstrates sensitivity for this set. It cannot identify the cause of every changed call. The GTEx direction label describes the reference tissue donors. These pooled CMap profiles do not estimate a treatment difference by donor sex. [Test results](phase4_literal/target_seed_tests.csv) · [All 88 sets](phase4_literal/archived_and_replay.csv)

All 603 pair-level signature counts agreed with their archived counts. The archive and replay each recorded 81 pairs without a hit. Eighty of those pairs were the same. Fostamatinib in placenta gained at least one hit, and thalidomide in pancreatic lines lost its last hit. The pair-level A/B/C counts changed for {four['pairs_with_changed_counts']} of {four['pairs']} pairs. Collapsing to the best category across all tested pairs changed {len(changed4)} of 50 compounds{(' (' + ', '.join(changed4) + ')') if changed4 else ''}. [Pair counts](comparison/phase4_pair_counts.csv) · [Compound counts](comparison/phase4_compounds.csv)

## Figures and manuscript

The [Phase 3 figure](figures/phase3_within_lineage_counts.png) puts a compound with zero hits at zero and marks the five untested compounds separately. The [Phase 4 figure](figures/phase4_tested_reference_sets.png) uses every tested set in each compound’s denominator. The archived chart mixed counts of enriched sets and single placeholders for pairs with zero hits. The [Phase 5 figure](figures/phase5_testable_tissues.png) distinguishes 620 completed tests from 39,376 untested tissue comparisons. One own-tissue and 79 other-tissue hits reproduced. Phase 5 pooled the expression profiles, so the sign of its NES cannot identify which donor group responded more.

The altered historical calls do not enter the manuscript's primary tests of expression state and drug response. Supplement Note S19's compound overlaps reproduced exactly (45/115 on P1 and P1b, 25/186 on P2). [Panel overlap](comparison/manuscript_note_s19.json)
"""
(HERE / 'RESULTS.md').write_text(text)
