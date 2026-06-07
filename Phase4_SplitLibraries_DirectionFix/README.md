# Phase 4 — Split Libraries + Direction Fix (Marouen's Feedback)

After Phase 3, Marouen sent three feedback emails. Phase 4 implements all of his corrections.

## Marouen's feedback summary

1. Don't compute DE per drug — level 5 data is already z-scores of differential expression.
2. Build male and female libraries separately per tissue (split by sign of effsize).
3. GSEA on z-scores directly with |z| > 2 cutoff for strong signals.
4. We are looking at drug-cell pairs — a drug can appear in multiple A/B/C categories simultaneously.
5. Direction interpretation bug: sex bias should come from which library (Male_X or Female_X) is enriched, NOT from NES sign.

## Sign convention verification (canonical female-biased X-linked genes)

| Gene | Effsize | Why female-biased |
|------|---------|-------------------|
| XIST | +10.07 | X-inactivation, canonical female marker |
| TSIX | +7.48 | XIST antisense, X-inactivation |
| KDM6A | +0.64 | Escapes X-inactivation |
| EIF1AX | +0.48 | Escapes X-inactivation |

All four are well-established female-biased X-linked genes. They all show POSITIVE effsize, confirming: positive effsize = FEMALE-biased.

## Results

| Metric | Phase 3 (mixed) | Phase 4 (split) |
|--------|-----------------|------------------|
| Cat A signals | 3 | 237 |
| Cat B signals | 23 | 8,942 |
| Cat C | 24 | 81 |
| Total classification rows | 50 | 9,260 |
| (Drug × lineage) pairs | 161 | 603 |
| Runtime | 38 min | 52 min |

## Direction split (the critical fix)

| Direction | Count | Percentage |
|-----------|-------|------------|
| Female | 5,104 | 55.5% |
| Male | 4,075 | 44.5% |

## Lapatinib verification (Marouen's specific question)

Phase 3 falsely labeled lapatinib × lung as "male-biased" (NES = +1.267 against mixed Lung library).

Phase 4 shows 8 FDR-significant libraries:
- 7 of 8 in Female direction (correct for breast cancer drug)
- Top hit: Female_Breast_Mammary_Tissue (NES = −1.531, FDR = 0.005)
- Other strong Female hits: Salivary Gland, Esophagus_Mucosa, Colon, Small_Intestine, Stomach
- Single Male hit: Pituitary (FDR = 0.038)

Marouen's hypothesis confirmed: lapatinib is female-biased in its response, consistent with its mechanism as a HER2/EGFR inhibitor used for breast cancer.

## Files

- Phase4_SplitLibraries.ipynb — the Colab notebook
- phase4_classifications_final.csv — all 9,260 classification rows
- phase4_stacked_per_drug.png — visual summary stacked bars per drug

## Author

Mateenah Jahan — Fatima Fellowship 2026

Mentor: Marouen Ben Guebila, PhD — Dana-Farber Cancer Institute / Harvard
