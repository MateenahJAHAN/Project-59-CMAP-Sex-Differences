# Part A raw-data replay

Independent extraction of the EmeraldBay expression counts recovered **100,253 endpoint-vehicle cells** from 52 lines and **30,245 time-zero cells**. Every checked field in all 130,498 per-cell rows matched the archived table. All 207 line-by-sample rows were byte-identical. No line changed its state, Y-arm class or confidence. The largest numerical difference in the line table was 3.47 × 10⁻¹⁸.

All three sci-Plex matrices were read in full, covering **1,171,291,549 sparse matrix entries**. The table of 23,570 vehicle cells matched the archived compressed file byte-for-byte. Five of six line-by-experiment marker patterns matched the freshly derived bulk RNA pattern. A549 expressed both Y markers and XIST in sciPlex3 and only Y markers in bulk.

The DNA check re-parsed 1,102 Cellosaurus records. All 19 STR output tables matched their archived checksums. Joining the newly calculated RNA scores directly to the re-parsed STR calls confirmed 47/52 in EmeraldBay and 3/3 distinct sci-Plex lines. New bulk RNA measurements retain all 773 CCLE and 1,082 Cell Model Passports line calls. Agreement with DNA typing was 714/734 for CCLE and 1,007/1,031 for Cell Model Passports. [Calls by accession](tables/str_fresh_bulk.tsv)

Tahoe RNA and reference-stock DNA typing agreed in 48/50 lines. The [five-source figure input](tables/figure3_agreement_by_source.tsv) and its [per-line decisions](tables/figure3_supp_table_amelogenin_by_line.tsv) were rebuilt from the new RNA measurements and re-parsed Cellosaurus records.

The pinned Tahoe-100M observation file passed its SHA-256 check. All 100,648,790 rows matched the local metadata projection across the eleven source fields used by the analysis. There are 2,205,786 full-filter endpoint-vehicle cells. Across 1,910 expression shards, 4,212 eligible row groups contained 2,205,786 qualifying cells. All 7,640 per-shard scientific files matched the archived files byte-for-byte. The [21 aggregated files](receipts/tahoe_aggregate_comparison.json) also matched, including the vehicle-cell marker counts and pseudobulk expression. All 700 line–plate counts matched the provider metadata. The [source receipt](receipts/tahoe_raw_shards.json) records the read counts and diagnostic differences.

The primary Tahoe caller reproduced the 50-line state table, its sensitivity and identity checks, all per-gene scores and its mosaic tests. All [33 freshly calculated tables](receipts/primary_caller_comparison.json) matched their archived bytes. Of 49 lines with a donor-sex annotation, 27 lacked its usual Y/XIST pattern. 685 of 700 plate-level calls agreed with their line. The Y/XIST pattern matched bulk RNA for 47/49 lines. Across Tahoe and EmeraldBay the pattern matched in 42/47 shared lines and 40/41 with high or medium confidence in both. [Tahoe state calls](tables/partA_local_states.tsv) · [Bulk comparisons](tables/patterns_tahoe_bulk.tsv) · [Shared cell lines](tables/patterns_cross_atlas.tsv)

Across all 1,831,648 EmeraldBay cells, every source file matched its deposited SHA-256. For 1,326 pairs of line profiles over 11,120 genes, the median correlation was 0.790565 (maximum 0.953835). [Pairwise correlations](tables/emerald_line_pair_correlations.tsv) · [116 source checks](receipts/emerald_all_condition_coverage.tsv)

The original positional comparison paired different cell lines in three near-tied pairs. Checking by accession reproduced every deep-cell count and nearest-line identity. The largest absolute numerical difference across the 17 scientific tables was 1.03e-11. [Differences by field and accession](receipts/external_difference_details.json)

The second Tahoe caller independently reduced 2,205,786 cells and retained 2,175,916 after its marker filter. The two callers agreed on the state in 50/50 lines, state and confidence in 50/50, and Y-arm state in 50/50. The alternate line-level cutoffs agreed with the primary state in 48/50. [Consensus calls](tables/tahoe_partA_states_consensus.tsv) · [Line-level sensitivity](tables/linelevel_sensitivity_differences.tsv) · [Mosaic tests](tables/mosaic_resolution.tsv)

The secondary caller's mixture fit gave different last digits for two lower bounds in the 1,694-row Supplementary Table S1. For LoVo Y, the archived, current and dependency-pinned values were 0.008423, 0.008421 and 0.008420. For SW48 XIST they were 0.259381, 0.259388 and 0.259378. All 22 binary mosaic decisions reproduced in both fits. A549 Y and SW48 XIST continue to meet the mosaic rule in only one of the two callers. All eight secondary diagnostic plots were pixel-identical when fresh and archived values were rendered using the same renderer. [Two endpoint changes](tables/table_S1_mosaic_sensitivity.tsv) · [22-fit software sensitivity](receipts/mosaic_version_comparison.json)

## Files and commands

The tables directory contains the complete line-level and EmeraldBay sample-level results. The receipts include the 116 EmeraldBay shard records, all 16 GEO file URLs and hashes, STR record identifiers and checksums, and aggregate and per-cell comparisons.

Run `python3 scripts/verify.py` to verify the package and recompute the reported counts. The original sci-Plex sparse-matrix parser and scientific constants are in `scripts/`. Place the 16 source files listed in `receipts/source_verification.json` under `RUN/data/F1/sciplex/`. From this directory run `python3 scripts/process_sciplex.py RUN sciPlex2` and repeat with `sciPlex3` and `sciPlex4`. This recreates the per-cell vehicle records, sample labels and identity accumulators for each full matrix. The parser requires NumPy, pandas, SciPy and PyArrow.

To reproduce the EmeraldBay deep-cell identity analysis, use a scratch directory for the downloaded expression shards. Each shard is checked against its deposited SHA-256 and released from scratch space after reduction. The three pinned metadata files and the 116 source hashes are in the receipts. The commands use Python with NumPy, pandas, SciPy, PyArrow and requests.

```sh
python3 scripts/prepare_emerald.py RUN
python3 scripts/process_emerald.py RUN 2
python3 scripts/explore_similarity.py RUN
```
