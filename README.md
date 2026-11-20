# Project 59 · Historical CMap and GTEx analyses

**Fatima Fellowship, 2026 · Mateenah Jahan**

Mentor: Dr. Marouen Ben Guebila, Viswanathan Lab at Dana-Farber Cancer Institute and Harvard Medical School

## Current status

This is the public archive of the fellowship’s exploratory CMap and GTEx work. The Phase 6 enrichments and Phase 7 sampling audit are historical. The earlier single-candidate conclusion has been superseded. The current manuscript and supplement are in coauthor review. [Fellowship project and status](https://mateenahjahan.github.io/#project-59)

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

| Material | Link |
| --- | --- |
| Phase 6 notebook and code | [Notebook](Phase6_NoSampleFilter/notebook/phase6_no_sample_filter.ipynb) · [Code](Phase6_NoSampleFilter/code/) |
| Phase 6 outputs | [Results tables](Phase6_NoSampleFilter/results/) · [Figure archive](Phase6_NoSampleFilter/figures/) |
| Phase 7 audit | [Notebook](Phase7_SexBalance_Audit/notebook/phase7_sex_balance_audit.ipynb) · [Audit tables](Phase7_SexBalance_Audit/results/) |
| Data sources and setup | [CMap and GTEx source notes](Phase6_NoSampleFilter/sources/README.md) |
