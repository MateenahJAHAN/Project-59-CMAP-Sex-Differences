# Project 59 — Phase 6: No Sample Filter

**Fatima Fellowship 2026 · Mateenah Jahan**  
Mentor: Dr. Marouen Ben Guebila, Viswanathan Lab, Dana-Farber Cancer Institute / Harvard

Does a compound change the expression of genes that already differ by sex in human tissues? Phase 6 screened 50 compounds by comparing drug-induced expression signatures from Connectivity Map (CMap LINCS2020) with tissue-specific sex-biased gene sets from GTEx v8.

## What “No Sample Filter” means

Phase 6 removed the requirement for at least five male-derived and five female-derived cell lines **within each drug–lineage pair**. The same minimum still applied when selecting compounds, pooling all lineages for each compound. The analysis also kept the gene-level threshold of `|z| > 2`.

## What the manuscript reports

- The screen analyzed **36,523 signatures** from **135 cell lines** across **20 lineages** and **44 tissues**, for **39,996 drug–lineage–tissue comparisons**. Seventeen of 50 compounds yielded at least one significant comparison under the screen’s criteria.
- At the **signature level**, 2,656 (7.3%) inherited a significant result in the cell line’s own tissue, 5,163 (14.1%) in a different tissue, 4,684 (12.8%) had no significant result, and 24,020 (65.8%) were not testable.
- **Repeated profiling drove the apparent female skew.** In the own-tissue category, 91.3% of signatures came from female-derived lines. Those signatures came from 17 female-derived and 18 male-derived lines. The female-derived lines were profiled about 11 times as often per line, so the signature counts overstated the evidence for a sex difference.
- **One candidate survived the cell-line audit.** In the bortezomib case study, the female-biased gene set for EBV-transformed lymphocytes was enriched among suppressed genes (NES −1.606, FDR < 0.0001). The male-biased set did not reach significance (NES +0.999, FDR 0.487). This comparison used 255 signatures across 15 cell lines, with a Kish effective sample size of 6.6.

The screen used one combined sex-biased gene set per tissue. The bortezomib case study used separate female- and male-biased sets. The manuscript describes the findings as hypothesis-generating. The screen lacked a negative control and a study-wide multiple-testing correction, and the audit measured repeated sampling without correcting it in the primary analysis.

## Phase 6 files

| Material | Link |
| --- | --- |
| Analysis overview | [Phase 6 overview](Phase6_NoSampleFilter/README.md) |
| Notebook and scripts | [Notebook](Phase6_NoSampleFilter/notebook/phase6_no_sample_filter.ipynb) · [Code](Phase6_NoSampleFilter/code) |
| Results and figures | [Results tables](Phase6_NoSampleFilter/results) · [Figures](Phase6_NoSampleFilter/figures) |
| Data sources and setup | [CMap and GTEx source notes](Phase6_NoSampleFilter/sources/README.md) |
| Cell-line audit | [Audit tables](Phase7_SexBalance_Audit/results) |
