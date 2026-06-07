# Project 59 — CMAP × GTEx Sex-Biased Drug Response

> Four-phase computational pipeline that screens 50 drugs for tissue-specific
> sex differences in gene expression. Uses public CMAP/LINCS 2020 and GTEx v8 data.

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![Colab](https://img.shields.io/badge/run%20on-Colab-orange.svg)](https://colab.research.google.com/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Fatima Fellowship](https://img.shields.io/badge/Fatima%20Fellowship-2026-purple.svg)](https://www.fatimafellowship.com/)

---

## 🎯 The question

Drug studies usually report a single average effect. But male and female cells often respond differently. This project asks:

> **Does a drug's sex-specific effect on gene expression match the tissue's natural sex biology?**

We screen 50 most-tested drugs from CMAP against 44 GTEx tissues and classify each (drug × cell-lineage × tissue) combination into one of three categories:

| Label | Meaning |
|-------|---------|
| **A** | Drug modulates sex-biased genes in its OWN tissue → real signal |
| **B** | Drug modulates sex-biased genes in a DIFFERENT tissue → off-target |
| **C** | No significant tissue-specific sex effect |

---

## 📁 Four phases, four folders

The project progressed in four stages. Each phase has its own folder with the notebook, the results, the figures, and a step-by-step README.

### 📂 [Phase1_Vorinostat/](Phase1_Vorinostat/)
**Single-drug pilot.** Built and tested the male-vs-female differential expression pipeline on one drug (vorinostat, an HDAC inhibitor).
- Found **20 sex-biased genes** (7 male-biased, 13 female-biased)
- Validated the method by recovering **CDKN1A** (a known HDAC target)

### 📂 [Phase2_Fishers_50Drugs/](Phase2_Fishers_50Drugs/)
**Scaled to 50 drugs with Fisher's exact test.** Added the GTEx sex-biased gene library and the A/B/C classification.
- 49 of 50 drugs → Category C (no signal)
- 1 of 50 → Category A: AS-605240 (but only 2 DE genes → **likely artifact**)
- Honest finding: pooled-across-lineages analysis dilutes tissue-specific signal

### 📂 [Phase3_PerLineage_GSEA/](Phase3_PerLineage_GSEA/)
**Per-lineage analysis with GSEA.** Followed Marouen's recommendation — run DE separately per cell lineage and use GSEA (`gseapy.prerank`) instead of Fisher's exact test.
- **3 clean Category A drugs**: tozasertib, lapatinib, vorinostat
- **23 Category B drugs** as follow-up candidates
- Phase 2's "Cat A" hit (AS-605240) correctly falls to Cat C

### 📂 [Phase4_SplitLibraries_DirectionFix/](Phase4_SplitLibraries_DirectionFix/)
**Marouen's three follow-up corrections.** Split GTEx libraries into male and female sets per tissue (88 sets total), used level-5 z-scores directly (no extra DE), and read sex direction from library identity instead of NES sign.
- **237 Cat A signals** (vs 3 in Phase 3) — split libraries recover signal
- **9,260 total classification rows** across 50 drugs × multiple lineages
- Lapatinib correctly reclassified as **female-biased** (Female_Breast hit)
- Sign convention verified using XIST, TSIX, KDM6A, EIF1AX
- **This is the definitive analysis of the project.**

---

## 🏆 Headline result

![Phase 3 classification](Phase3_PerLineage_GSEA/phase3_classification_chart.png)

| Method | Cat A | Cat B | Cat C |
|--------|-------|-------|-------|
| Phase 2 (Fisher's, pooled) | 1 (artifact) | 0 | 49 |
| Phase 3 (GSEA, per-lineage, mixed libraries) | 3 clean | 23 | 24 |
| **Phase 4 (GSEA, per-lineage, split M/F libraries)** | **237** | **8,942** | **81** |

### Top Category A drugs

| Drug | Lineage | GTEx tissue | NES | FDR | Source |
|------|---------|-------------|-----|-----|--------|
| **Tozasertib** | haematopoietic | Spleen | −1.475 | 0.001 | Phase 3 |
| **Tozasertib** | haematopoietic | Whole_Blood | −1.342 | 0.007 | Phase 3 |
| **Tozasertib** | haematopoietic | EBV lymphocytes | −1.193 | 0.016 | Phase 3 |
| **Lapatinib** | lung | Female_Breast_Mammary_Tissue | −1.531 | 0.005 | **Phase 4** |
| **Vorinostat** | large_intestine | Colon_Transverse | −1.189 | 0.029 | Phase 3 |

---

## 🚀 Quick start

```bash
git clone https://github.com/MateenahJAHAN/Project-59-CMAP-Sex-Differences.git
cd Project-59-CMAP-Sex-Differences
pip install -r requirements.txt
```

Then open any phase folder and run the notebook in Google Colab:

| Folder | Notebook | Runtime |
|--------|----------|---------|
| `Phase1_Vorinostat/` | `Project59_CMAP_SexDifferences.ipynb` | ~20 min |
| `Phase2_Fishers_50Drugs/` | `Project59_GTEx_DrugClassification.ipynb` | ~30 min |
| `Phase3_PerLineage_GSEA/` | `Phase3_PerLineage_GSEA.ipynb` | ~40 min |
| `Phase4_SplitLibraries_DirectionFix/` | `Phase4_SplitLibraries.ipynb` | ~55 min |

For step-by-step explanations and results of each phase, open the **README.md inside that phase's folder**.

---

## 📚 Data sources

### CMAP / LINCS 2020 (clue.io)
- `cellinfo_beta.txt` — cell line metadata (donor sex, lineage)
- `siginfo_beta.txt` — drug experiment metadata
- `geneinfo_beta.txt` — gene ID translator (Entrez ↔ Ensembl)
- `level5_beta_trt_cp_n720216x12328.gctx` — z-score matrix, 720,216 signatures × 12,328 genes (33 GB)

### GTEx v8 (gtexportal.org)
- `GTEx_Analysis_v8_sbgenes.tar.gz` — sex-biased genes per tissue (44 tissues × 13,294 genes, LFSR ≤ 0.05)
- Split into 88 male/female libraries in Phase 4 using the `effsize` column

All inputs are public and free.

---

## 🔬 What changed across phases

| | Phase 2 | Phase 3 | Phase 4 |
|---|---------|---------|---------|
| DE strategy | Pooled all lineages | Per-lineage Mann-Whitney | Use level-5 z-scores directly |
| Library | 44 mixed sb-genes | 44 mixed sb-genes | 88 split (M + F per tissue) |
| Enrichment test | Fisher's exact (binary) | GSEA (ranked) | GSEA (ranked) |
| Direction read | NES sign | NES sign | **Library identity** |
| Result | 1 fake Cat A | 3 clean Cat A | 237 Cat A |

For full method detail, open each phase's README.

---

## 🛣️ Next steps (Phase 5 ideas)

1. **Per-signature analysis** — run GSEA per individual signature instead of averaging z-scores across signatures within a (drug, lineage).
2. **Tighten the |z| > 2 cutoff** at the signature level (less averaging means the filter is more selective).
3. **Increase GSEA permutations** to 1,000 (Phase 4 used 200 for speed).
4. **Top Female-biased drugs as a class** — do they share mechanism, target, or scaffold?
5. **Cross-reference Cat A findings** with published drug-sex literature to look for known matches.

---

## 📑 Citation

```
Jahan, M. (2026). Per-lineage GSEA with split male/female libraries for
sex-biased drug response in CMAP/LINCS 2020 against GTEx v8 sex-biased genes.
Fatima Fellowship 2026, Project 59.
Mentor: Dr. Marouen Ben Guebila (Dana-Farber).
GitHub: https://github.com/MateenahJAHAN/Project-59-CMAP-Sex-Differences
```

---

## 👤 Author

**Mateenah Jahan** — Fatima Fellowship 2026
Mentor: [**Marouen Ben Guebila, PhD**](https://labs.dana-farber.org/viswanathanlab/people/marouen-ben-guebila-phd) (Viswanathan Lab, Dana-Farber Cancer Institute)

---

## 📜 License

MIT License. See [LICENSE](LICENSE).
