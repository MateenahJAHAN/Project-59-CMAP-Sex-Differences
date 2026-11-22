# Inferred Y-gene expression in LINCS2020 controls

Can an inferred L1000 Y score distinguish lines classified by external RNA sequencing? Across 36 XY and 19 Y-lost lines, its AUROC was **0.534** (95% bootstrap CI 0.348 to 0.709 in the original diagnostic). The correlation with the CCLE Y score was **Spearman ρ = 0.213** across the 124 lines with both measurements (two-sided p = 0.0174). These estimates describe Level 3 DMSO controls. [Analysis receipt](results/receipt.json) · [Bootstrap and gene-level replay checks](results/all_gene_order_check.json)

A second comparison includes every sampled line with a resolved Y-high or Y-low call. It gives **AUROC 0.651** across 37 Y-high and 90 Y-low lines. Only the 36-versus-19 comparison underlies the manuscript’s 0.53.

## Data and sampling

The published matrix contains 188,708 control profiles and 12,328 genes. Of 3,026,460 instrument records, 75,705 met the DMSO vehicle and quality criteria. Sampling at most 30 wells per line with seed 20261114 recovered all 2,800 original wells across 127 lines. We read their complete rows from the Broad matrix. The inferred score is the mean of eight gene-specific medians within each line. For each complete 127-line table, the four-decimal marker medians matched the original results. [Broad Level 3 control matrix](https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level3/level3_beta_ctl_n188708x12328.gctx) · [Reference table fingerprints](reference/reference_fingerprints.json)

## Reproduce

Use Python 3.12 and a local Task A state table generated from CCLE and Cell Model Passports. The script reads the matrix metadata and selected profile rows over HTTPS. It saves the aggregate results and SHA-256 checks for the complete selected cohort, its source rows, and its marker medians.

```sh
python -m pip install -r requirements.txt
mkdir -p data
curl -fL https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/geneinfo_beta.txt -o data/geneinfo_beta.txt
curl -fL https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/instinfo_beta.txt -o data/instinfo_beta.txt
python run.py --data data --states /path/to/A_cellline_sex_states.tsv --reference reference --output rerun
python verify.py --output rerun --states /path/to/A_cellline_sex_states.tsv
```

The Task A table supplies `cell_iname`, `state`, `Yhi` and `y_score_ccle` in the original Task A row order. The rank correlation and both AUROCs are computed locally across every eligible line. To check the saved aggregate receipt offline, run `python verify.py`.
