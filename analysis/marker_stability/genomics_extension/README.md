# Autosomal drug responses on repeat Tahoe-100M plates

The mean same-line correlation rose from 0.201 for different drugs to 0.385 for the same drug. Among distinct lines of the same organ, sharing baseline Y/XIST marker status added 0.0109 correlation units (exact two-sided P = 0.204). [Methods, results and figure legend](section.md) · [Figure](landmark_responses.png) · [Vector figure](landmark_responses.svg) · [Numerical results and sensitivity analyses](results.json)

## Reproduce the results

Run these commands from this directory. The included integer count matrix is sufficient to regenerate every test, table and panel. The analysis was run with Python 3.12.11, NumPy 2.5.4, pandas 3.1.1 and Matplotlib 3.11.2. Re-extraction from the full source also uses PyArrow 26.0.1.

```sh
OPENBLAS_NUM_THREADS=1 python analyse.py
OPENBLAS_NUM_THREADS=1 python check.py
python figure.py
```

The independent check reconstructs normalized expression from the compact counts and reruns both crossed Pearson correlations for every eligible condition of every pair. It then enumerates the 249,480 within-organ marker reassignments and compares the effect and P value with the main analysis. [Independent numerical check](independent_check.json)

## Inputs and coverage

- [Protocol](protocol.md), [checksum](protocol.sha256) and [feature eligibility amendment](feature_amendment.md). Chromosome mapping excluded the initially listed TSC22D3 from the glucocorticoid control before treatment outcomes were extracted.
- [Source hashes](sources.json), [all 978 Broad landmark IDs and chromosome assignments](broad_landmark_mapping.tsv) and [retained features](features.tsv). Gene coordinates come from Ensembl GRCh38 release 110. Each retained feature has a unique matching Ensembl ID and Tahoe token. The compact file contains the 948 autosomal landmarks plus three additional control transcripts.
- [Integer counts](counts.npz), [well metadata](wells.tsv) and [line metadata](lines.tsv). The axes of `gene_counts` are well, line and gene in the row orders of the accompanying TSVs. `umi` holds all-gene library totals and `cells` holds captured cell counts. The 192 wells contained 13,486,933 captured cells across all 50 input lines. Cell counts do not give the independent sample size for the marker test.
- [Coverage by drug–dose](coverage_by_condition.tsv), [97 eligible line pairs](distinct_line_pairs.tsv), [same-line correlations](same_line_reproducibility.tsv), [individual glucocorticoid measurements](triamcinolone_by_line.tsv), and [permutation reference quantiles](permutation_quantiles.tsv). NCI-H661 and SW1271 fell below the vehicle minimum. L-thyroxine at 5 µM did not meet the treatment minimum in any paired line. The two plate layouts repeated the inferred drug positions.

To rebuild the count matrix from the original pseudobulks, pass the Tahoe archive and the two reference tables to `extract.py --tahoe <archive> --geneinfo <geneinfo_beta.txt> --chrom <gene_chrom.tsv>`. The expected input hashes are in [sources.json](sources.json). The included integer count matrix supports the expression analyses without those upstream downloads. [Tahoe-100M](https://doi.org/10.1016/j.cell.2026.08.035) · [Broad landmark definition](https://doi.org/10.1016/j.cell.2017.10.049)
