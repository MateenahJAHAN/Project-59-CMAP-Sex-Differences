# Phase 3: Within-lineage enrichment

This exploratory phase ranked genes within CMap drug–lineage groups and used preranked gene-set enrichment against the combined GTEx v8 tissue sets. It collapsed the outputs to one category per compound. The archived table classifies 3 compounds as A (own-tissue), 23 as B (other tissue) and 24 as C (no significant enrichment under these rules).

The notebook originally read positive GTEx effect sizes as male-biased and used the sign of the normalized enrichment score to assign sex direction. Those labels are unreliable. In GTEx v8, positive `effsize` represents higher expression in females. In GSEA, the sign indicates where a set falls in the ranked gene list. Phase 4 separated the reference sets by the GTEx sign. The stored Phase 3 enrichment results remain available with their original labels.

[Notebook](Phase3_PerLineage_GSEA.ipynb) · [Compound classifications](phase3_classification_final.csv) · [Category figure](phase3_classification_chart.png)

[Next: split tissue gene sets](../Phase4_SplitLibraries_DirectionFix/README.md) · [Fellowship overview](../README.md#how-the-analysis-developed)
