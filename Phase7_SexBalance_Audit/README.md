# Phase 7: Cell-line sampling audit

This audit counted distinct cell lines and repeated profiles behind the Phase 6 categories. Several own-tissue enrichments came from lineages represented by only one donor-sex group.

| Category | Female-donor signatures | Male-donor signatures | Distinct female-donor lines | Distinct male-donor lines |
| --- | ---: | ---: | ---: | ---: |
| Own-tissue enrichment (A) | 2,426 | 230 | 17 | 18 |
| Other-tissue enrichment (B) | 2,958 | 2,205 | 59 | 52 |
| Tested without significant enrichment (C) | 1,710 | 2,974 | 60 | 53 |
| Untested (D) | 15,719 | 8,301 | 76 | 57 |

A line can appear in more than one category. Female-donor lines accounted for 91.3% of Category A signatures and contributed about 11 times as many profiles per line as male-donor lines. The bortezomib breast-line pair alone supplied 1,875 of the 2,656 Category A signatures, all from female-donor lines. [Signature counts](results/signatures_by_category_and_sex.csv) · [Cell-line counts](results/cell_lines_by_category_and_sex.csv)

Eleven drug–lineage pairs contributed Category A signatures. Bortezomib in blood lines and YM-155 in skin each had at least two lines of each sex. One skin line contributed 21 of the 28 YM-155 signatures (75.0%). For bortezomib in blood lines, the most profiled line contributed 61 of 255 (23.9%). Its 15 lines had a Kish effective sample size of 6.60. The audit specified no numerical cutoff for concentration. The audit describes sampling support. The pooled enrichment does not estimate a sex-specific bortezomib response. [Donor-sex composition by pair](results/category_A_by_drug_lineage_sex.csv) · [Concentration and effective sample size](results/pseudoreplication_category_A.csv)

The replay changed the best category for eight of 36,523 signatures, all from idarubicin-treated stomach lines. Five female-donor and three male-donor profiles moved from C to B. The 2,656 category A assignments remained the same. The [archived signature-label checksum](results/archived_signature_labels.json) covers the full October export projected to `sig_id` and `sig_cat`, sorted by signature ID as tab-delimited text with a header and LF line endings.

The [archived Colab notebook](notebook/phase7_sex_balance_audit.ipynb) reads the original Phase 6 signature and pair outputs. The result tables above are available directly. A line-weighted reanalysis would require rebuilding the rankings with each cell line as the unit.

[Phase 6 screen](../Phase6_NoSampleFilter/README.md) · [Next: FAERS adverse-event reporting](../Phase8_SideEffects_FAERS/README.md) · [Project and manuscript](../README.md)
