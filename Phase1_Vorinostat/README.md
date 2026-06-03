# Phase 1 — Vorinostat Single-Drug Pilot

> First phase. Built and validated the male-vs-female differential expression
> pipeline on one well-studied drug: **vorinostat** (an FDA-approved HDAC
> inhibitor for cutaneous T-cell lymphoma).

---

## 🎯 The goal

Before scaling to many drugs, test the pipeline end-to-end on one drug:

1. Can we pull a drug's signatures from the 33 GB CMAP gctx file?
2. Can we split them by donor sex and run a fair statistical test?
3. Does the result recover known biology (sanity check)?

---

## 📋 What I did — step by step

### Step 1 — Download CMAP metadata

I downloaded three small text files from the Broad Institute's clue.io
servers:

- `cellinfo_beta.txt` — 240 cell lines with donor sex labels
- `siginfo_beta.txt` — metadata for each drug experiment
- `geneinfo_beta.txt` — gene ID translator

These let me figure out which cell lines are male-donor vs female-donor
and which signatures belong to which drug.

### Step 2 — Filter to drug-treatment experiments

I kept only rows where:
- `pert_type == 'trt_cp'` (a real drug treatment, not a control)
- The cell line's donor sex was known (M or F, not Unknown)

This left **596,117 valid drug experiments** in the dataset.

### Step 3 — Pick the drug

I needed a drug tested in both sexes with enough samples. **Vorinostat**
was chosen because:
- It is FDA-approved (so well-studied)
- It was tested in **116 cell lines** (53 male donors + 63 female donors)
- It is a histone deacetylase (HDAC) inhibitor — a class with known
  chromatin-level effects, good for validation

### Step 4 — Download the big z-score matrix

The CMAP `level5_beta_trt_cp_n720216x12328.gctx` file is **33 GB**. It
holds z-scores for 720,216 drug-treatment signatures × 12,328 genes.

Using `cmapPy.pandasGEXpress.parse`, I extracted **only vorinostat's 116
columns** — no need to load the full 33 GB.

### Step 5 — Split into male and female matrices

I matched each signature ID with its donor sex and split the matrix:
- Male submatrix: 12,328 genes × 53 experiments
- Female submatrix: 12,328 genes × 63 experiments

### Step 6 — Mann-Whitney U test per gene

For each gene, I compared its z-score distribution between male and
female cells using the Mann-Whitney U test. This is a non-parametric
test — robust to outliers and does not assume normal distributions.

Output per gene: difference of medians, raw p-value.

### Step 7 — FDR correction

Running 12,328 tests at once would give many false positives just by
luck. I applied the Benjamini-Hochberg FDR correction to get adjusted
p-values that account for multiple testing.

### Step 8 — Apply effect-size filter

Statistical significance alone is not enough — a tiny difference can be
significant if you have many samples. I also required a meaningful
biological effect:

```
strong DE gene = FDR < 0.05 AND |z-score difference| > 1.0
```

### Step 9 — Annotate with gene symbols

Using `geneinfo_beta.txt`, I mapped each gene's Entrez ID to its
readable symbol (e.g., 1026 → CDKN1A).

### Step 10 — Visualise

- **Volcano plot:** effect size on x-axis, −log10(FDR) on y-axis. The
  strong DE genes appear at the corners.
- **Top-20 bar chart:** the 20 standout genes ranked by absolute effect size.

---

## 📊 Results

### Numbers

| Filter | Genes |
|--------|-------|
| Total genes tested | 12,328 |
| Statistically significant (FDR < 0.05) | 6,877 |
| Strong DE (FDR < 0.05 AND \|Δ\| > 1.0) | **20** |
| └ stronger in male cells | 7 |
| └ stronger in female cells | 13 |

### Key finding

**CDKN1A** appears in the top 20. CDKN1A (also called p21) is a well-known
transcriptional target of HDAC inhibitors — recovery of this gene confirms
the pipeline is detecting real biology and not noise.

### What this tells us about vorinostat

- Vorinostat causes **sex-specific gene expression changes** in 20 genes
  that survive both statistical and effect-size filters.
- None of the 20 genes are sex-chromosome-encoded — they are autosomal.
  So the difference is not coming directly from X or Y chromosome activity.
- The sex difference must arise from **downstream regulatory networks** that
  process the HDAC inhibition differently between male-donor and
  female-donor cell backgrounds.

### What this phase did NOT do

- Did not compare findings against tissue-level natural sex biology (GTEx).
- Did not check if the sex effect is tissue-specific or generic.
- Used pooled signatures across all vorinostat cell lineages.

Those questions were addressed in Phase 2 and Phase 3.

---

## 📁 Files in this folder

| File | What it is |
|------|------------|
| `Project59_CMAP_SexDifferences.ipynb` | The Colab notebook for Phase 1 |
| `vorinostat_all_genes_results.csv` | Full per-gene results (12,328 rows) |
| `vorinostat_top20_sex_genes.csv` | The 20 strong DE genes with symbols |
| `vorinostat_signatures.csv` | The 116 signature IDs used |
| `volcano_vorinostat.png` | Volcano plot |
| `top20_genes_vorinostat.png` | Bar chart of top 20 |

---

## 🛠️ How to reproduce

1. Open `Project59_CMAP_SexDifferences.ipynb` in Google Colab.
2. Run all cells top to bottom.
3. The big gctx download takes about 15 minutes.
4. Full analysis takes ~20 minutes total.
5. Output CSVs and PNGs save to the Colab working directory.

---

## ➡️ Why Phase 2 was needed

Phase 1 answered "does sex matter for one drug?" but not "does this sex
effect line up with a specific tissue's natural sex biology?"

To answer that, we needed:
- A reference list of sex-biased genes per tissue → **GTEx sbgenes**
- A way to compare drug effects vs tissue reference → **enrichment test**
- More drugs to find real patterns → **scale to top 50**

That is what Phase 2 does. Open the `Phase2_Fishers_50Drugs/` folder.
