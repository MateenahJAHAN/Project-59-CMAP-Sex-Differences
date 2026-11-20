# Phase 6: Full CMap and GTEx screen

This exploratory screen pooled CMap perturbation signatures by drug and cell lineage, then compared the ranked genes with GTEx tissue gene sets. Removing the five-lines-per-sex requirement within a drug–lineage pair made more pairs eligible. The requirement still applied when selecting the 50 compounds across all lineages.

The batch kept at least three signatures per pair and genes with absolute mean z-score above 2. GSEA required at least 15 matching genes per set. The batch used one combined sex-biased gene set per GTEx tissue.

## Screen results

The panel covered 36,523 signatures from 135 cell lines across 20 lineages. Comparing the 909 drug–lineage pairs with 44 tissues produced 39,996 candidate comparisons.

| Category | Meaning | Comparisons | Signatures |
| --- | --- | ---: | ---: |
| A | Enrichment in a GTEx tissue mapped to the cell lineage | 29 | 2,656 |
| B | Enrichment in another GTEx tissue | 717 | 5,163 |
| C | Tested, with no significant enrichment | 5,376 | 4,684 |
| D | No test outcome | 33,874 | 24,020 |
| Total | | 39,996 | 36,523 |

The first count uses drug–lineage–tissue comparisons. Each signature inherits the best category of its drug–lineage pair across tissues, in order A, B, C, D. Seventeen compounds reached A or B, including six that reached A. Category D includes lineages with no mapped GTEx tissue, pairs with too few signatures or ranked genes, and gene sets that did not yield GSEA results. [Category totals](results/category_totals.csv) · [Per-compound counts](results/signature_ABCD_per_drug.csv) · [Untested comparison reasons](results/category_D_skip_reasons.csv)

The table above gives the recorded run. A replay using GSEApy 1.3.2 returned 718 B and 5,375 C comparisons, with the A and D totals unchanged. The counts of compounds reaching A or B and reaching A were the same.

[![Phase 6 signature counts by compound on linear and logarithmic axes.](figures/01_signatures_per_drug_linear_log.png)](figures/01_signatures_per_drug_linear_log.png)

## Follow-up comparisons

The four saved analyses cover three drug–lineage–tissue combinations. They used separate female- and male-biased GTEx sets. The bortezomib blood-line comparison against the EBV-transformed lymphocyte tissue set used the full mapped ranking (12,276 distinct Ensembl IDs from 255 signatures and 15 cell lines). The two hypothalamus analyses, one pooled and one restricted to BJAB, and the emetine lung-line analysis used the absolute mean z-score cutoff of 2. The emetine split-set comparison did not reach q < 0.05 for either set. The category in the CSV refers to the earlier pooled batch comparison. [Case tables](results/case_studies_split_library.csv) · [Enrichment figures](figures/)

The bortezomib EBV follow-up gave NES −1.606 for the female-biased set and +0.999 for the male-biased set. GSEA stored the corresponding empirical q values as 0.0 and 0.4874 using 1,000 permutations. A zero in the output reflects the estimator’s finite resolution.

## Limits and inputs

These comparisons measure overlap with GTEx reference sets in pooled ranks. The screen has no direct donor-sex contrast, negative control or multiple-testing correction across the whole screen. Phase 7 found that two own-tissue pairs had at least two cell lines of each sex. Only the bortezomib blood-line pair also passed the concentration check. [Phase 7 audit](../Phase7_SexBalance_Audit/README.md)

The [Colab notebook](notebook/phase6_no_sample_filter.ipynb) preserves the original run. The [batch export](code/01_batch_abcd.py) uses the original Drive layout and derived inputs. The [case-study export](code/02_case_study_plots.py) includes notebook shell commands and Colab setup. The [source notes](sources/README.md) list the inputs and explain how to prepare them. The saved batch file name `per_signature_v4_ABCD.csv` refers to the drug–lineage–tissue table and is distinct from the signature-level output.

[Project and manuscript](../README.md#current-status) · [Earlier phases](../README.md#how-the-analysis-developed)
