# Phase 2 — 50 Drugs with Fisher's Exact Test

> Second phase. Scaled the Phase 1 pipeline from one drug to 50. Added the
> GTEx sex-biased gene library and Marouen's A/B/C classification framework.
> Used **Fisher's exact test** for enrichment — a simpler method that
> turned out to be too coarse (Phase 3 fixes this).

---

## 🎯 The goal

Phase 1 showed that vorinostat has sex-specific gene effects. Phase 2 asks
two new questions:

1. Do other drugs also show male-vs-female differences?
2. Do those differences match the **natural sex biology of a specific tissue**?

To answer #2 we needed a reference: which genes are naturally sex-biased
in each body part? That comes from GTEx.

---

## 📋 What I did — step by step

### Step 1 — Download GTEx sex-biased genes

I downloaded `GTEx_Analysis_v8_sbgenes.tar.gz` from gtexportal.org. After
unzipping, the main file is `signif.sbgenes.txt` with 104,864 rows.

Each row is one (gene, tissue) pair that is significantly sex-biased
(LFSR ≤ 0.05). Columns: gene, tissue, effect size, standard error,
reliability score.

### Step 2 — Build the sb-gene library

I converted the flat table into a Python dictionary:

```
{ "Lung":   {gene1, gene2, ...},
  "Liver":  {geneA, geneB, ...},
  "Brain_Cortex": {...},
  ... 44 tissues total }
```

Sizes range from 473 genes (lymphocytes) to 4,558 (sun-exposed skin).

### Step 3 — Build the CMAP ↔ GTEx translator

CMAP labels cell lines using broad lowercase names (`lung`, `liver`,
`central_nervous_system`). GTEx uses formal names (`Lung`, `Liver`,
`Brain_Cortex`).

I built a manual mapping:

```
'lung'                              → ['Lung']
'breast'                            → ['Breast_Mammary_Tissue']
'skin'                              → [2 Skin tissues]
'liver'                             → ['Liver']
'large_intestine'                   → [2 Colon tissues]
'haematopoietic_and_lymphoid_tissue'→ [3 blood tissues]
'central_nervous_system'            → [13 Brain regions]
... etc
```

**Coverage:** 70 % of known-lineage CMAP cell lines (111 of 158).
Lineages like ovary, prostate, endometrium, bone have no GTEx
equivalent and remain unmapped.

### Step 4 — Filter to top 50 drugs

From all 32,917 drugs in CMAP, I kept only those with:
- ≥ 5 male-donor cell lines AND
- ≥ 5 female-donor cell lines

That left 2,408 eligible drugs. I picked the **top 50 by total cell line
count** for maximum statistical power.

### Step 5 — Per-drug pipeline (looped over 50 drugs)

For each drug, the function `analyze_drug()` did:

1. Pulled all the drug's signatures from gctx (**pooled across all lineages**).
2. Ran Mann-Whitney M vs F for each of 12,328 genes.
3. Applied FDR correction → strong DE genes (FDR < 0.05, |Δ| > 0.5).
4. Mapped Entrez gene IDs to Ensembl IDs (for matching with GTEx).
5. Ran **Fisher's exact test** against each of the 44 GTEx tissue
   sb-gene sets:
   - Universe = all genes with Ensembl ID (12,276)
   - Tested if the overlap between drug DE genes and tissue sb-genes
     was bigger than chance
6. Applied FDR correction across the 44 tests per drug.
7. Assigned A / B / C based on which tissues passed (if any).

### Step 6 — Visualise

I built two charts:
- A 50-drug horizontal bar chart, coloured by category
- A 44-tissue enrichment chart for vorinostat specifically

---

## 📊 Results

### Category breakdown

| Category | Count |
|----------|-------|
| A (same-tissue match) | **1** |
| B (cross-tissue match) | 0 |
| C (no significant tissue signal) | 49 |

### The one Cat A drug — and why it is suspicious

**AS-605240** was flagged as Category A:
- Best tissue: Liver
- Fold enrichment: 28×
- FDR adjusted p: 0.025

The fold value of 28× sounds impressive. **But:** AS-605240 had only
**2 strong DE genes**. With only 2 genes in the overlap, Fisher's
exact test gives **inflated fold values** that look real but are
actually small-numbers statistical artifacts.

**My honest interpretation:** This is not real biology. The hit is
mathematically real but biologically meaningless. I reported it
honestly to Marouen with this caveat, and the new Phase 3 method
correctly classifies AS-605240 as Category C (no signal) once the
proper GSEA test is used.

### Why 49 of 50 drugs landed in Category C

Two reasons, both fixable:

**Reason 1 — Pooling diluted tissue-specific signal.**

I pooled all of a drug's cell lines into one M vs F test. So vorinostat's
lung cells got mixed with brain cells, skin cells, liver cells, etc.
A real lung-specific sex effect gets washed out by signal from other
tissues. The test sees only an average across tissues, which is usually
weak.

**Reason 2 — Fisher's exact test loses information.**

Fisher's only looks at the strong-DE list (a binary "in or out" per gene).
Many genuinely sex-affected genes have p-values just above the cutoff
and get discarded. GSEA, which uses the full ranked list, would catch
those subtle consistent shifts.

Both problems are fixed in Phase 3 — see `Phase3_PerLineage_GSEA/`.

---

## 📁 Files in this folder

| File | What it is |
|------|------------|
| `Project59_GTEx_DrugClassification.ipynb` | The Colab notebook for Phase 2 |
| `top50_drug_classification.csv` | 50 drugs × {category, n_strong_DE, best_tissue, fold, p} |
| `vorinostat_tissue_enrichment.csv` | Vorinostat × 44 tissue Fisher's results |
| `top50_drugs_strong_de.png` | Bar chart of all 50 drugs |
| `vorinostat_tissue_enrichment.png` | Vorinostat × 44 tissue fold-enrichment chart |

---

## 🛠️ How to reproduce

1. Open `Project59_GTEx_DrugClassification.ipynb` in Google Colab.
2. Run all cells top to bottom.
3. GTEx download is fast (31 MB). gctx download is 33 GB (~15 min).
4. The 50-drug loop takes ~13 minutes after data is loaded.
5. Total runtime: about 30 minutes.

---

## ➡️ Why Phase 3 was needed

Phase 2 gave one honest negative result (49/50 Cat C, 1 fake Cat A) but
the method was clearly too blunt. Two specific changes were obvious:

1. **Stop pooling.** Run DE per cell lineage separately.
2. **Use GSEA, not Fisher's.** Use the full ranked list, not a binary cutoff.

That's exactly what Phase 3 does. Open `Phase3_PerLineage_GSEA/` to see
the proper version and the 3 clean Cat A drugs it found.
