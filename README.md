# Project 59 — CMAP Sex Differences

This project looks at how the same drug can affect male and female cells differently, using public data from the CMAP / LINCS 2020 dataset. The drug studied here is vorinostat (a histone deacetylase inhibitor).

The whole analysis runs inside one Google Colab notebook: `Project59_CMAP_SexDifferences.ipynb`.

---

## Goal of the project

When scientists test a drug, they usually report one average effect. But male and female cells can react differently. The goal here is simple: find which genes respond to vorinostat differently in male vs female cells.

These sex-different genes can help explain why some drugs work better, or cause more side effects, in one sex than the other.

---

## Data used

- `cellinfo_beta.txt` — metadata about each cell line (sex, tissue, cell type, etc.)
- `siginfo_beta.txt` — info about each drug experiment (drug, cell line, dose, time)
- `level5_beta_trt_cp_n720216x12328.gctx` — the big gene-expression matrix (~33 GB) with z-scores for every drug experiment across 12,328 genes
- `geneinfo_beta.txt` — dictionary that maps gene IDs to readable gene symbols

All files come from the official clue.io LINCS 2020 public S3 bucket.

---

## Steps in the notebook

1. Setup — install `cmapPy`, import pandas and numpy.
2. Download metadata — cellinfo, siginfo, and geneinfo files.
3. Clean the cell info — drop rows where sex is Unknown, keep only Male / Female cell lines.
4. Pick the drug — find drugs tested in both sexes; choose vorinostat.
5. Download the big matrix — the 33 GB `.gctx` file with all gene-expression signatures.
6. Filter signatures — keep only vorinostat experiments, split them into a Male matrix and a Female matrix.
7. Statistical test — run a Mann–Whitney U test for every gene to compare male vs female responses.
8. Multiple-testing correction — apply FDR (Benjamini–Hochberg) to control false positives.
9. Effect size filter — keep genes where the median difference between male and female is large enough to matter (|Δ| ≥ 0.5).
10. Add gene names — map gene IDs to readable symbols using `geneinfo_beta.txt`.
11. Save results — export two CSV files: all genes tested + top sex-different genes.
12. Visualize — volcano plot (effect size vs `-log10(FDR)`) and bar chart of the top sex-different genes.

---

## Key result

The notebook produces a final list of genes that respond to vorinostat in a significantly different way between male and female cells (FDR < 0.05 and a meaningful effect size). These genes are the candidates for follow-up biological interpretation.

Output files saved by the notebook:

- `vorinostat_sex_diff_all_genes.csv` — full per-gene table with stats
- `vorinostat_sex_diff_top_genes.csv` — only the significant, large-effect genes
- `volcano_plot.png` — volcano plot of the results
- `top_genes_bar.png` — bar chart of the strongest sex-different genes

---

## A quick note on -log10(p)

P-values can be very small numbers like 0.0000001. If you plot them as-is, all the important genes squish together near zero.

`log10` just counts how many zeros a number has. Taking `-log10(p)` flips small p-values into large positive numbers, so the most significant genes appear high up on the volcano plot.

| p-value | -log10(p) |
|---------|-----------|
| 0.05    | 1.3       |
| 0.001   | 3         |
| 1e-6    | 6         |

Higher on the y-axis = stronger statistical evidence of a sex difference.

---

## How to run

1. Open `Project59_CMAP_SexDifferences.ipynb` in Google Colab.
2. Run the cells from top to bottom.
3. The 33 GB `.gctx` download needs a stable connection and enough Colab disk space — a high-RAM runtime is recommended.
4. Final CSVs and plots are saved in the Colab working directory and can be downloaded at the end.

---

## Tools used

- Python 3
- cmapPy — for reading `.gctx` files
- pandas, numpy — data handling
- scipy.stats — Mann–Whitney U test
- statsmodels — FDR correction
- matplotlib — plots

---

## Author

**Mateenah Jahan** — Project 59

Supervised by [**Marouen Ben Guebila, PhD**](https://labs.dana-farber.org/viswanathanlab/people/marouen-ben-guebila-phd) (Viswanathan Lab, Dana-Farber Cancer Institute).
