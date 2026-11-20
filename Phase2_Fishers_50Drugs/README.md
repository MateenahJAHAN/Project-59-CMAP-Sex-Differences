# Phase 2: Fifty-compound screen

This exploratory phase expanded the vorinostat pilot to 50 compounds. It mapped CMap gene IDs to Ensembl IDs and used one-sided Fisher’s exact tests to examine overlap between drug gene lists and GTEx v8 tissue sets. The screen pooled signatures across lineages and used annotated donor sex.

The vorinostat example mapped 12,277 of 12,328 records to 12,276 distinct Ensembl IDs. Among the features with adjusted p < 0.05, 6,854 distinct Ensembl IDs mapped to the library. For the lung test, 26 of the 270 genes that passed the additional absolute mean-difference threshold of 0.5 overlapped the lung set in the test universe. Fisher’s exact p was 0.574.

The 50-compound table records one own-tissue category assignment and 49 no-enrichment assignments under this phase’s rules. The single assignment was based on two genes, which limits its interpretation. Phase 3 moved the comparison within lineages.

[Notebook](Project59_GTEx_DrugClassification.ipynb) · [Fifty-compound classifications](top50_drug_classification.csv) · [Vorinostat tissue results](vorinostat_tissue_enrichment.csv)

[Next: within-lineage enrichment](../Phase3_PerLineage_GSEA/README.md) · [Fellowship overview](../README.md#how-the-analysis-developed)
