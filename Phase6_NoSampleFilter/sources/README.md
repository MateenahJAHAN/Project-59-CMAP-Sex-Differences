# Phase 6 data and parameters

## CMap LINCS2020

The level-5 compound treatment file contains 720,216 signatures and 12,328 L1000 features, including measured and inferred genes. Use the metadata to select signature IDs before reading the expression matrix. A CMap signature represents a perturbation profile, often derived from multiple profiles in the earlier processing levels. [CMap data portal](https://clue.io/data/CMap2020#LINCS2020) · [L1000 paper](https://doi.org/10.1016/j.cell.2017.10.049)

| File | Contents |
| --- | --- |
| `level5_beta_trt_cp_n720216x12328.gctx` | External level-5 expression matrix. |
| `siginfo_beta.txt` | External signature metadata, including cell line, dose, time and quality flags. |
| [Cell metadata](../../cmap/cellinfo_beta%20%282%29.txt) | Archived donor-sex annotation and lineage. |
| [Gene metadata](../../cmap/geneinfo_beta%20%281%29.txt) | CMap feature IDs, Ensembl IDs and symbols. |
| [Compound metadata](../../cmap/compoundinfo_beta%20%281%29.txt) | Compound identifiers and annotations. |

The [fellowship downloader](../../fellowship/sources.py) reads the study signatures from Broad by byte range and verifies each against its recorded hash. Earlier notebooks accessed the GCTX through [cmapPy](https://github.com/cmap/cmapPy).

## GTEx v8

Use `signif.sbgenes.txt` from the GTEx v8 sex-biased expression archive. It contains 104,864 gene–tissue rows across 44 tissues. Positive `effsize` means higher expression in females. Negative `effsize` means higher expression in males. The GTEx table uses a local false sign rate threshold of 0.05. These tissue gene sets are references for comparison with cancer cell lines. [GTEx input](../../gtex/README.txt) · [GTEx study](https://doi.org/10.1126/science.aba3066)

## Reference inputs and settings

- [Lineage map](../../Phase5_PerSignature_Categories/reference/cmap_to_gtex.json) links CMap lineages to GTEx tissues.
- [GTEx gene-set JSON](../../Phase5_PerSignature_Categories/reference/sb_library_sets.json) has one combined set per tissue. Convert its lists to Python sets before intersections.
- [Compound list](../../Phase5_PerSignature_Categories/reference/top50_names.json) defines the Phase 6 panel.
- The notebook builds `sig_known.csv` from LINCS signature metadata and cell-line metadata. It builds `gid_to_ensg.json` by dropping unmapped rows from the gene table and mapping CMap feature IDs to Ensembl IDs. The gene table has 12,277 mapped feature records and 12,276 distinct Ensembl IDs.

The batch and case-study exports use 1,000 GSEA permutations, minimum set size 15 and seed 42. The batch and the individual case studies use different rankings and reference sets, as described in the [Phase 6 report](../README.md). The original notes record GSEApy 1.3.1. The subsequent replay used 1.3.2. [Install and run the pinned versions](../../fellowship/README.md#run-the-fellowship-analyses) · [GSEA method](https://doi.org/10.1073/pnas.0506580102)
