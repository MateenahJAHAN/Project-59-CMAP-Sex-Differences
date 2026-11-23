# What counts as an agreement?

In Tahoe-100M and external bulk RNA-seq, **47 of 49** comparable lines had the same measured Y/XIST pattern. COLO 205 was below the Tahoe Y cutoff and above the bulk cutoff. NCI-H596 showed the reverse. When we compare **the names assigned to the states**, 46 of 49 agree. A498 is XIST-only in both cultures. The Tahoe rule calls XIST-only expression in a male-labelled line discordant. The bulk rule calls any XIST-only pattern XX-like. Cellosaurus records A498 as male and Cell Model Passports records it as female. [Three lines underlying the two comparisons](results/tahoe_bulk_disagreements.tsv)

## Two denominators

Of 596 male-labelled and 302 female-labelled GDSC2/PRISM models eligible by lineage and annotation, 17 and two lacked RNA measurements. Both markers were low for 256 male-labelled and 122 female-labelled lines. Counting low Y alone adds 16 XIST-only male-labelled lines, giving 272/579 with measured RNA. Counting low XIST alone adds three Y-only female-labelled lines, giving 125/300.

For LINCS2020, 129 donor-labelled names had RNA measurements. HCT116 had opposing Y calls in CCLE and Sanger, leaving 128 resolved state names. Of these, 61 were Y-lost, XIST-lost or discordant. U266B1 received the name XY in both sources because Y was high. XIST was high in CCLE and low in Sanger. Consequently, 127 lines had both combined marker flags, of which 61 had an atypical pattern. [HCT116 and U266B1 source calls](results/lincs_exceptions.tsv)

The [claim table](results/claims.tsv) preserves the earlier 128-line wording for the LINCS RNA denominator alongside its corrected value. It also gives the two GDSC2/PRISM contrasts separately. [Completed source checks](results/coverage.json) · [Input SHA-256 hashes](results/source_sha256.json)

## Supplementary Table S1

All 1,694 accession–source pairs and 1,123 distinct lines were recovered. Of 40,656 table cells, two lower bounds of estimated marker-negative share changed. Using the current libraries, LoVo moved from 0.008423 to 0.008421 and SW48 from 0.259381 to 0.259388. Rerunning the same 22 fits with the study’s pinned libraries gave 0.008420 and 0.259378. [Dependency comparison](../part_a/receipts/mosaic_version_comparison.json) The state, confidence, arm call and subpopulation verdict were identical for every record. [Cell-level comparison](results/s1_cell_comparison.json) · [Independent source roster](results/s1_roster.json)

## Reproduce

From the repository root after running Parts A and B, use `python3 reproduce.py claims --work ../p59-run`. This rebuilds the claim comparisons and checks them against the released audit. A single historical discrepancy is retained in both records.

For a workspace with the older raw-data directory layout, run

```sh
python run.py --run-root /path/to/raw-reproduction-workspace --reference-s1 /path/to/published/supplementary_table_S1.tsv --output rerun
```

The script reads the two bulk reference tables, LINCS calls, vehicle-cell statistics, STR calls and available downstream outputs from the reproduction tree. Use `--supplement /path/to/supplement.tex` if the manuscript source for Table S8 lives elsewhere. The [Level 3 L1000 analysis](../l1000/README.md) independently checks the inferred-score correlation and AUROC.
