# Inferred Y-gene expression in LINCS2020 controls

Can an inferred L1000 Y score distinguish lines classified by external RNA sequencing? Across 36 XY and 19 Y-lost lines, its AUROC was **0.534** (95% bootstrap CI 0.348 to 0.709 in the original diagnostic). The correlation with the CCLE Y score was **Spearman ρ = 0.213** across the 124 lines with both measurements (two-sided p = 0.0174). These estimates describe the Level 3 DMSO controls.

A second comparison includes every sampled line with a resolved Y-high or Y-low call. It gives **AUROC 0.651** across 37 Y-high and 90 Y-low lines. The original and expanded comparisons use different groups. Only the 36-versus-19 comparison underlies the manuscript’s 0.53.

## Data and sampling

The matrix has 188,708 control profiles and 12,328 genes. From 3,026,460 instrument records we selected the 75,705 passing quality control and annotated as DMSO vehicle. We reproduced the archived sampling rule (at most 30 wells per line, seed 20261114), recovered the same 2,800 wells in 127 lines and read their complete rows from the Broad matrix. The score is the mean of eight gene-specific medians within each line. The marker manifest, matrix identifiers, byte counts, and eligibility counts are in the [run receipt](results/receipt.json). [Broad Level 3 control matrix](https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level3/level3_beta_ctl_n188708x12328.gctx)

[Well identifiers and source-row SHA-256 hashes](results/sampled_well_hashes.tsv.gz) · [Inferred scores for 127 lines](results/inferred_scores.tsv) · [Agreement of the 12,328-gene exploratory table after joining on gene ID](results/all_gene_order_check.json)

## Reproduce

Run these commands from this directory using Python 3.12. The analysis downloads only the matrix metadata and selected rows. The two metadata files are stored locally.

```sh
python -m pip install -r requirements.txt
mkdir -p data
curl -fL https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt -o data/geneinfo_beta.txt
curl -fL https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/instinfo_beta.txt -o data/instinfo_beta.txt
python run.py --data data --states reference/line_state_inputs.tsv --reference reference --output rerun
python verify.py --output rerun
```

The reference directory contains the Task A fields used for sampling and comparison, the original well selections, and the archived marker medians rounded to four decimal places. A full Task A state table with the same column names can also be supplied. To verify the included evidence offline, run `python verify.py`.
