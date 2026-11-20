# Phase 4: Split GTEx gene sets

This exploratory phase averaged the CMap level-5 z-scores within each drug–lineage pair, retained genes with absolute mean z-score above 2, and tested separate female- and male-biased GTEx sets. The identity of the enriched reference set and the sign of the normalized enrichment score describe different parts of the result. These pooled responses cannot estimate a difference between donor-sex groups.

In GTEx v8, positive `effsize` denotes higher expression in females. XIST (ENSG00000229807) has positive values in all 44 reported tissues, including +10.07 in lung. The corrected code places XIST in the female-biased lung set. One working notebook shows an initial reversed split and the subsequent rebuild before enrichment.

- [Notebook with the corrected split](Phase3_PerLineage_GSEA%20new.ipynb)
- [Second saved copy of the corrected code](Phase3_PerLineage_GSEA%20%283%29.ipynb)
- [Working notebook showing the initial split and correction](%3D%3D%3D%3D%3D.ipynb)
- [Classification output](phase4_classifications_final.csv) · [Other saved export](phase4_classifications_final%20%281%29.csv) · [Per-compound category figure](phase4_stacked_per_drug.png)

Both classification exports contain the same 9,260 drug–lineage–reference-set rows in different order (237 A, 8,942 B and 81 C). The Phase 3 table counts compounds, so its counts use a different denominator. The lapatinib lung example enriched seven female-biased sets among eight significant sets in that comparison.

[Next: signature categories](../Phase5_PerSignature_Categories/README.md) · [Fellowship overview](../README.md#how-the-analysis-developed)
