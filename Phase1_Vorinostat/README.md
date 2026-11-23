# Phase 1: Vorinostat pilot

This exploratory screen compared level-5 CMap LINCS2020 gene expression after vorinostat treatment across signatures grouped by annotated donor sex. The test used **2,641 signatures** from **116 distinct cell lines**: 1,180 signatures from 53 male-donor lines and 1,461 from 63 female-donor lines. Cell lines could contribute repeated signatures at different doses and times.

The notebook ran a two-sided Mann–Whitney test for each of 12,328 genes and applied Benjamini–Hochberg correction across those tests. In the saved results, 6,877 genes passed adjusted p < 0.05 and 20 also had an absolute difference in mean z-score above 1. The original test used signatures as observations and pooled across lineages, doses and times. These results describe that pilot comparison and cannot isolate a sex effect.

Phase 2 mapped the significant gene records to **6,854 distinct Ensembl IDs**. This is a different count from the 6,877 CMap features that passed the Phase 1 threshold.

[Notebook](Project59_CMAP_SexDifferences.ipynb) · [Signature list](vorinostat_signatures.csv) · [All gene results](vorinostat_all_genes_results.csv) · [Top 20 genes](vorinostat_top20_sex_genes.csv)

[Next: fifty-compound screen](../Phase2_Fishers_50Drugs/README.md) · [Fellowship phases and reproduction](../fellowship/README.md)
