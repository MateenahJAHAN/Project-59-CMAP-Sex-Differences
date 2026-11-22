# Phase 5: Categories for perturbation signatures

This exploratory notebook builds the CMap-to-GTEx lineage map, prepares the 50-compound panel, and assigns each signature the best category recorded for its drug–lineage pair across tissues. The four categories distinguish own-tissue enrichment, other-tissue enrichment, a tested pair with no significant enrichment, and an untestable comparison.

The file named [per_signature_categories.csv](results/per_signature_categories.csv) has one row per drug, lineage and tissue. It records 620 completed comparisons across 19 drug–lineage pairs (1 A, 79 B, 540 C), with 39,376 comparisons untested. The sole own-tissue enrichment was bortezomib in blood lines against the GTEx EBV-transformed lymphocyte set. In the [pair summary](results/per_drug_lineage.csv), the `analysed` flag is true only for four pairs with no recorded skip reason. To find every pair with at least one test result, count those with A + B + C greater than zero.

Start with the [notebook](project59_pipeline.ipynb) and [reference inputs](reference/). The output folder includes a [multi-sheet export with the pair-level table](results/FULL_50drugs_table.xlsx), [per-drug counts](results/all50_drug_counts.xlsx), and [working exports at different levels of aggregation](results/). The additional workbooks cover different rows or values and are kept with their original filenames. Use the Phase 6 tables for counts after the per-pair cell-line minimum was removed.

[Next: full screen](../Phase6_NoSampleFilter/README.md) · [Fellowship overview](../README.md#how-the-analysis-developed)
