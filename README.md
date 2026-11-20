# Project 59 · Sex-chromosome expression and drug response

**Fatima Fellowship, 2026 · Mateenah Jahan**

Mentor: Dr. Marouen Ben Guebila, Viswanathan Lab at Dana-Farber Cancer Institute and Harvard Medical School

## Manuscript and figures

[**Manuscript (PDF)**](manuscript/manuscript.pdf) · [**Supplement (PDF)**](manuscript/supplement.pdf) · [**All figure PDFs**](manuscript/figures/)

Mateenah Jahan and Marouen Ben Guebila. Both authors approved the November 16, 2026 version 2 manuscript and supplement for release as public drafts. The figure PDFs are from the same version.

[![Figure 1. Study design for the expression-state measurement and drug-response test.](manuscript/figures/fig1_design.png)](manuscript/figures/fig1_design.pdf)

[![Figure 4. Prespecified drug-response test, including the TP53 benchmark and response-size estimates.](manuscript/figures/fig4_partB.png)](manuscript/figures/fig4_partB.pdf)

[Figure 2](manuscript/figures/fig2_label_vs_pattern.pdf) · [Figure 3](manuscript/figures/fig3_agreement.pdf) · [Supplementary figure S1](manuscript/figures/figS1_ambient.pdf) · [S2](manuscript/figures/figS2_a549.pdf) · [S3](manuscript/figures/figS3_power.pdf)

## Earlier CMap and GTEx analyses

The notebooks and results below document the exploratory fellowship phases that preceded the manuscript. The manuscript supersedes the earlier single-candidate interpretation. [Fellowship overview](https://mateenahjahan.github.io/#project-59)

## How the analysis developed

The project started with one compound, vorinostat. It then asked how drug-induced gene expression in CMap LINCS2020 cell lines overlaps with gene sets from GTEx v8 tissues. The unit of analysis changed as the screen grew. By Phase 6 the enrichment compared pooled responses with tissue gene sets. Phase 7 counted the individual lines behind those responses.

| Phase | Work in the archive |
| --- | --- |
| 1 · Vorinostat | Pilot analysis of one drug. [Phase 1](Phase1_Vorinostat/README.md) |
| 2 · Fifty compounds | Broader screen using Fisher’s exact test and tissue reference gene sets. [Phase 2](Phase2_Fishers_50Drugs/README.md) |
| 3 · Lineages | Gene-set enrichment within cell lineages. [Phase 3](Phase3_PerLineage_GSEA/README.md) |
| 4 · Gene sets | Separate female- and male-biased libraries and a check of the ranking direction. [Phase 4](Phase4_SplitLibraries_DirectionFix/README.md) |
| 5 · Signatures | Assigned each perturbation signature the category inherited from its drug–lineage comparisons. [Phase 5 notebook](Phase5_PerSignature_Categories/project59_pipeline.ipynb) |
| 6 · Full screen | Dropped the five-lines-per-sex rule within each drug–lineage pair. The rule still applied during compound selection across lineages, and the screen kept the gene-level threshold of `|z| > 2`. [Phase 6](Phase6_NoSampleFilter/README.md) |
| 7 · Sampling audit | Counted distinct cell lines behind the Phase 6 enrichments and measured how heavily the screen relied on individual lines. [Phase 7](Phase7_SexBalance_Audit/README.md) |

## What the counts show

Phase 6 included **50 compounds, 36,523 signatures, 135 cell lines, 20 lineages and 44 GTEx tissues**. There were 39,996 candidate drug–lineage–tissue comparisons. Of these, 6,122 had a recorded test outcome and 33,874 were Category D. Seventeen compounds had at least one significant comparison under the historical criteria. At the signature level, 2,656 inherited a significant result in their own tissue, 5,163 in a different tissue, 4,684 inherited no significant result and 24,020 were untestable. [Category totals](Phase6_NoSampleFilter/results/category_totals.csv) · [Signature category figure](Phase6_NoSampleFilter/figures/04_signature_ABCD_2panel.png)

![Historical Phase 6 signature counts for 50 compounds. Counts range from 241 to 7,300 and are shown on linear and logarithmic axes.](Phase6_NoSampleFilter/figures/01_signatures_per_drug_linear_log.png)

Female-derived lines accounted for **91.3% of signatures** in the own-tissue category. Those signatures came from **17 female-derived and 18 male-derived lines**. The female-derived lines had about 11 times as many profiles per line. One historical bortezomib case met the audit’s representation rule, which counted available lines and the concentration of measurements. [Sampling audit](Phase7_SexBalance_Audit/README.md) · [Audit tables](Phase7_SexBalance_Audit/results/)

The screen used one combined sex-biased gene set per tissue and lacked a negative control or a multiple-testing correction across the full screen. The bortezomib case study used separate female- and male-biased gene sets. Establishing a sex effect would require a direct comparison between the cell-line groups and appropriate controls. [Historical bortezomib figure](Phase6_NoSampleFilter/figures/05_bortezomib_CellsEBV_enrichment_catA.png)

## Code and data

| Phase | Notebooks, notes and outputs |
| --- | --- |
| 1 · Vorinostat | [Notes](Phase1_Vorinostat/README.md) · [Notebook](Phase1_Vorinostat/Project59_CMAP_SexDifferences.ipynb) · [Gene-level results](Phase1_Vorinostat/vorinostat_all_genes_results.csv) |
| 2 · Fifty compounds | [Notes](Phase2_Fishers_50Drugs/README.md) · [Notebook](Phase2_Fishers_50Drugs/Project59_GTEx_DrugClassification.ipynb) · [Compound classifications](Phase2_Fishers_50Drugs/top50_drug_classification.csv) |
| 3 · Lineages | [Notes](Phase3_PerLineage_GSEA/README.md) · [Notebook](Phase3_PerLineage_GSEA/Phase3_PerLineage_GSEA.ipynb) · [Classifications](Phase3_PerLineage_GSEA/phase3_classification_final.csv) |
| 4 · Gene sets | [Notes](Phase4_SplitLibraries_DirectionFix/README.md) · [Archived notebooks and figures](Phase4_SplitLibraries_DirectionFix/) · [Classifications](Phase4_SplitLibraries_DirectionFix/phase4_classifications_final.csv) |
| 5 · Signatures | [Notebook](Phase5_PerSignature_Categories/project59_pipeline.ipynb) · [Results tables](Phase5_PerSignature_Categories/results/) |
| 6 · Full screen | [Notes](Phase6_NoSampleFilter/README.md) · [Notebook](Phase6_NoSampleFilter/notebook/phase6_no_sample_filter.ipynb) · [Code](Phase6_NoSampleFilter/code/) · [Results tables](Phase6_NoSampleFilter/results/) · [Figures](Phase6_NoSampleFilter/figures/) |
| 7 · Sampling audit | [Notes](Phase7_SexBalance_Audit/README.md) · [Notebook](Phase7_SexBalance_Audit/notebook/phase7_sex_balance_audit.ipynb) · [Results tables](Phase7_SexBalance_Audit/results/) |
| Data and setup | [CMap](cmap/README.txt) · [GTEx](gtex/README.txt) · [Phase 6 source notes](Phase6_NoSampleFilter/sources/README.md) · [Dependencies](requirements.txt) |
