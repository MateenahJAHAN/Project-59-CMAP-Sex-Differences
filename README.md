Project 59 — CMAP Sex Differences
This project looks at how the same drug can affect male and female cells differently, using public data from the CMAP / LINCS 2020 dataset. The drug studied here is vorinostat (a histone deacetylase inhibitor).

The whole analysis runs inside one Google Colab notebook: Project59_CMAP_SexDifferences.ipynb.

Goal of the project
When scientists test a drug, they usually report one average effect. But male and female cells can react differently. The goal here is simple: find which genes respond to vorinostat differently in male vs female cell lines.

Because cell lines have no hormones, any difference found is genetic (sex chromosomes and baseline expression), not hormonal. These sex-different genes are candidates for follow-up biological interpretation.

Data used
cellinfo_beta.txt — metadata about each cell line (donor sex, tissue, cell type)
siginfo_beta.txt — info about each drug experiment (drug, cell line, dose, time)
geneinfo_beta.txt — dictionary that maps gene IDs to readable gene symbols
level5_beta_trt_cp_n720216x12328.gctx — the big gene-expression matrix (~33 GB), z-scores for 12,328 genes across 720,216 drug-treatment signatures
All files come from the official clue.io LINCS 2020 public source.

Steps in the notebook
Setup — install cmapPy, import pandas and numpy.
Download metadata — cellinfo, siginfo, and geneinfo files.
Clean the cell info — drop cell lines with Unknown sex, keep only Male / Female.
Pick the drug — find drugs tested in both sexes; choose vorinostat.
Download the big matrix — the ~33 GB .gctx file.
Filter signatures — keep only vorinostat experiments, split them into a Male matrix and a Female matrix.
Statistical test — run a Mann–Whitney U test for every gene, comparing male vs female responses.
Multiple-testing correction — apply Benjamini–Hochberg FDR to control false positives.
Effect-size filter — keep genes where the mean difference between male and female is large (|difference| > 1.0 in z-score units), so the result is not just statistically significant but also meaningful.
Add gene names — map gene IDs to readable symbols using geneinfo_beta.txt.
Save results — export the per-gene results and the standout genes as CSV files.
Visualise — a volcano plot and a bar chart of the standout genes.
Key result
Out of 12,328 genes, 6,877 were statistically significant (FDR < 0.05). After adding the effect-size filter (|difference| > 1.0), this narrows to 20 standout genes — 7 respond more strongly in male cells and 13 more strongly in female cells.

CDKN1A, a known target of HDAC inhibitors, is among them — which supports that the pipeline is working correctly. None of the 20 genes are sex-chromosome genes, so the differences lie in downstream gene regulation.

Output files saved by the notebook:

vorinostat_all_genes_results.csv — full per-gene table with statistics for all 12,328 genes
vorinostat_top20_sex_genes.csv — the 20 significant, large-effect genes
vorinostat_signatures.csv — the list of vorinostat experiments used
volcano_vorinostat.png — volcano plot of the results
top20_genes_vorinostat.png — bar chart of the 20 standout genes
A quick note on -log10(p)
P-values can be very small numbers like 0.0000001. If you plot them as-is, all the important genes squish together near zero.

log10 just counts how many zeros a number has. Taking -log10(p) flips small p-values into large positive numbers, so the most significant genes appear high up on the volcano plot.

p-value	-log10(p)
0.05	1.3
0.001	3
1e-6	6
Higher on the y-axis = stronger statistical evidence of a sex difference.

How to run
Open Project59_CMAP_SexDifferences.ipynb in Google Colab.
Run the cells from top to bottom.
The ~33 GB .gctx download needs a stable connection and enough Colab disk space.
Final CSVs and plots are saved in the Colab working directory and can be downloaded at the end.
Tools used
Python 3
cmapPy — for reading .gctx files
pandas, numpy — data handling
scipy.stats — Mann–Whitney U test
statsmodels — FDR correction
matplotlib — plots
Author
Mateenah Jahan — Fatima Fellowship 2026, Project 59

Mentor: Marouen Ben Guebila, PhD — Dana-Farber Cancer Institute
