# Phase 3 — Per-Lineage Analysis with Proper GSEA

> Third and final phase. Implemented Marouen's actual recommended method:
> per-lineage differential expression + GSEA (`gseapy.prerank`). This is
> the **definitive analysis** of the project.

---

## 🎯 The goal

Phase 2 left two clear problems:

1. Pooling cell lines across lineages dilutes tissue-specific signal.
2. Fisher's exact test throws away rank information.

Phase 3 fixes both at once:

1. Run M vs F separately **within each cell lineage**.
2. Use **GSEA** on the full ranked gene list, not Fisher's on a binary set.

---

## 📋 What I did — step by step

### Step 1 — Install gseapy

Phase 2 used `scipy.stats.fisher_exact`. Phase 3 uses `gseapy`, the standard
Python implementation of Gene Set Enrichment Analysis.

```bash
!pip install gseapy
```

### Step 2 — Reload all metadata + the GTEx library

Same inputs as Phase 2: cellinfo, geneinfo, siginfo, the GTEx sb-gene
dictionary, the CMAP ↔ GTEx tissue translator. Same top 50 drugs.

### Step 3 — Write the per-lineage DE function

Built a new function `per_lineage_DE(drug_name, lineage)` that:

1. Filters siginfo to one drug × one lineage.
2. Checks that the lineage has at least 5 male AND 5 female cell lines
   (otherwise the test is too noisy — skip).
3. Pulls only those signatures from the gctx file.
4. Splits into male and female matrices.
5. Runs Mann-Whitney U test per gene.
6. Applies FDR correction.
7. **Builds a ranked list** with score = sign(Δ) × −log10(p_value):
   - Top of list = strongly male-biased
   - Bottom of list = strongly female-biased
   - Middle = no signal
8. Maps Entrez IDs to Ensembl IDs for downstream matching with GTEx.

**This is the key new piece.** A ranked list contains far more information
than a binary "strong / not strong" classification, and GSEA can use all
of it.

### Step 4 — Write the GSEA function

Built `run_gsea(de_results, gene_sets_dict)` that:

1. Takes the ranked list from step 3.
2. Takes the dictionary of 44 GTEx tissue sb-gene sets.
3. Calls `gseapy.prerank()`:
   - 1000 permutations (gives stable FDR estimates)
   - min set size 15 (gseapy requirement)
   - seed = 42 (reproducibility)
4. Returns a table with NES (normalised enrichment score) and FDR per tissue.

**Direction convention:**
- Positive NES → sb-genes enrich at the TOP → male-biased pattern
- Negative NES → sb-genes enrich at the BOTTOM → female-biased pattern

### Step 5 — Write the full-drug wrapper

Built `analyze_drug_full(drug_name)` that:

1. Lists all the drug's eligible lineages.
2. Runs `per_lineage_DE` for each.
3. Runs `run_gsea` for each lineage's DE result.
4. Collects all (lineage, tissue, NES, FDR) tuples where FDR < 0.05.
5. Sorts them as Cat A evidence (tissue matches drug's home lineage) or
   Cat B evidence (tissue does not match).
6. Decides final A/B/C:
   - Any Cat A evidence → **A**
   - Else any Cat B evidence → **B**
   - Else → **C**

### Step 6 — Test on vorinostat first

Vorinostat ran successfully across 3 lineages (lung, blood, colon) in
about 80 seconds. It came out as **Category A** — its large_intestine
lineage hit Colon_Transverse (NES = −1.189, FDR = 0.029).

This was the first sign that the new pipeline was finding real signals
that Phase 2's Fisher's test had missed.

### Step 7 — Run on all 50 drugs

The 50-drug loop took **~37 minutes**. Checkpoint saves every 5 drugs
in case of a crash.

### Step 8 — Inspect Cat A and Cat B in detail

Filtered the result table and saved separate CSVs for Cat A (3 drugs)
and Cat B (23 drugs) with the full evidence rows.

### Step 9 — Visualise

Built a horizontal bar chart of all 50 drugs, coloured by category
(green = A, yellow = B, gray = C). One picture summarises the whole
screen.

---

## 📊 Results

### Headline

| Category | Phase 2 (Fisher's) | Phase 3 (GSEA) |
|----------|--------------------|----------------|
| A | 1 (artifact) | **3 clean** |
| B | 0 | **23** |
| C | 49 | 24 |

---

### The 3 Category A drugs

#### 1. Tozasertib — strongest hit of the project

Aurora kinase inhibitor (research compound). In the **haematopoietic /
lymphoid** lineage, three blood tissues all showed significant enrichment
in the same direction:

| GTEx tissue | NES | FDR |
|-------------|-----|-----|
| Spleen | **−1.475** | **0.001** |
| Whole_Blood | −1.342 | 0.007 |
| EBV-transformed lymphocytes | −1.193 | 0.016 |

Negative NES = female-biased direction. Three tissues from the same
lineage giving the same answer is strong internal evidence that this is
real biology, not noise.

#### 2. Lapatinib — clean lung match

HER2 / EGFR inhibitor (FDA-approved for HER2-positive breast cancer).
In the lung lineage, Lung sb-genes enriched at the TOP of the ranked list
(positive NES = male-biased direction):

| GTEx tissue | NES | FDR |
|-------------|-----|-----|
| Lung | **+1.267** | **0.032** |

Consistent with reported HER2-driven sex differences in lung and breast
cancer outcomes.

#### 3. Vorinostat — borderline colon match

FDA-approved HDAC inhibitor. In the large_intestine lineage,
Colon_Transverse sb-genes enriched at the BOTTOM (negative NES):

| GTEx tissue | NES | FDR |
|-------------|-----|-----|
| Colon_Transverse | −1.189 | **0.029** |

Borderline (FDR just under the 0.05 cutoff) but biologically plausible —
HDAC inhibitors cause genome-wide chromatin changes that can interact
with tissue-specific transcriptional programmes including sex-biased ones.

---

### Phase 2's artifact correctly rejected

**AS-605240** was Phase 2's only "Category A" hit. Phase 3 classifies it
as **Category C** — no significant tissue enrichment in any of its
lineages.

This is the desired behaviour. The Phase 2 hit was based on only 2 strong
DE genes, and Fisher's exact test gives unreliable enrichment estimates
at such small numbers. Phase 3 uses the entire ranked list of 12,328
genes and sees no consistent signal. It correctly returns no result.

**This validates that the new method has higher specificity.** It avoids
overclaiming.

---

### 23 Category B drugs — future follow-up candidates

These drugs show clear sex-specific patterns but the significantly enriched
tissue is NOT their home tissue. They are candidates for investigating
off-target sex effects.

Top Cat B drugs by number of significant tissue hits:

| Drug | # significant tissues |
|------|-----------------------|
| AM-580 | 26 |
| vemurafenib | 20 |
| ABT-751 | 16 |
| barasertib-HQPA | 12 |
| idarubicin | 11 |
| calcitriol | 7 |
| nutlin-3 | 7 |

---

## 📁 Files in this folder

| File | What it is |
|------|------------|
| `Phase3_PerLineage_GSEA.ipynb` | The Colab notebook for Phase 3 |
| `phase3_classification_final.csv` | 50 drugs × {category, n_lineages, A_evidence_count, B_evidence_count} |
| `phase3_category_A_drugs.csv` | 3 Cat A drugs with full evidence rows |
| `phase3_category_B_drugs.csv` | 23 Cat B drugs with full evidence rows |
| `phase3_classification_chart.png` | Visual summary, all 50 drugs |

---

## 🛠️ How to reproduce

1. Open `Phase3_PerLineage_GSEA.ipynb` in Google Colab.
2. Run all cells top to bottom.
3. Heavy data is reused from Phase 2 (gctx, siginfo).
4. GSEA loop takes about 40 minutes for all 50 drugs (longer than Phase 2
   because each test runs 1000 permutations).

---

## ➡️ What I would do next (Phase 4 ideas)

1. **Validate tozasertib** on isolated blood cell lines (K562 for CML,
   MOLM13 for AML). The triple-tissue hit in this screen is convincing
   but a per-cell-line follow-up would strengthen the claim.
2. **Investigate top Cat B drugs.** AM-580 (26 hits) and vemurafenib (20)
   show the strongest off-target sex patterns. Possibly publishable.
3. **GSEA leading-edge analysis** for each Cat A hit — identify the
   specific genes driving the enrichment.
4. **Compare findings** against published sex-biased drug literature
   to see if any match known results.

---

## 🔑 Key methodological lessons

1. **Method matters more than data.** Same 50 drugs, completely different
   conclusion just by switching pooled → per-lineage and Fisher's → GSEA.
2. **Honest negatives are valuable.** AS-605240 falling to Cat C
   validates the new method's specificity. Knowing what is NOT real is
   as important as knowing what IS real.
3. **Per-lineage analysis is essential for tissue-specific questions.**
   You cannot ask "is this tissue-specific?" if you have already pooled
   away the tissue dimension.
4. **GSEA is more sensitive than Fisher's exact** for ranked data, but
   it needs ≥ 5 samples per group per lineage to be reliable.
