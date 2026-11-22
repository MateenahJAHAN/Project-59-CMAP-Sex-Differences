# LINCS drug-response analysis

The rebuild read all 30,056 selected drug signatures and 4,330 vehicle signatures from the LINCS2020 Level 5 matrices. The well-level analysis read 61,652 treated wells and 9,108 vehicle wells. All 749 plates had enough vehicle wells for same-plate subtraction. All 45 reconstructed array, axis and identifier files matched the original model inputs byte for byte. The starting measurements are the provider’s Level 3 normalised expression and Level 5 consensus signatures.

| Independent check | Result |
|---|---|
| Measured-state classification | All 129 LINCS line labels and eligibility decisions reproduced. Among 23 upstream Task A tables, 22 were byte-identical. KDE density and ratio diagnostics in the cutoff table differed by at most 1.11 × 10⁻¹⁵. All cutoffs and state assignments agreed. |
| Covariate join | On reloading the line table, six P2 ploidy values differed by at most 4.45 × 10⁻¹⁶. All 30 model and sensitivity designs agreed exactly. The ploidy values that enter the primary compound adjustment also agreed. |
| Ploidy benchmark | Fresh line tables and responses gave panel p values of 0.50445 (P1), 0.67293 (P1b) and 0.07929 (P2), each from 10,000 within-lineage permutations. Every recorded number matched. |
| Negative-control label assignments | All nine Task D and four Task D2 assignment banks had identical decompressed contents. |
| Transcription-factor sets | The 10,319 CollecTRI edges linking transcription factors to landmark genes and the 43 candidate transcription factors reproduced exactly. |
| Idarubicin, 12 µM, 24 h | The high-dose check recovered 51 signatures from 49 lines. Within the analysis cohort, 20 lines remained after the lineage rule, including only two Y-positive lines. The minimum was three per group, so the high-dose effect could not be estimated. |

The control results, membership checks and SHA-256 comparisons are provided as tables. The source annotation inventory links to the provider files. The manuscript’s primary inference requires the complete permutation, interval and compound-level analyses.

## Reproduction and restart

Run `python3 verify.py` to check the integrity of this package. The 15 [LINCS analysis modules](../../../analysis/lincs/src/G/) match the frozen code. The [public configuration](../../../analysis/lincs/config/) uses the same numerical settings. The complete analysis checkout at commit `87640df15a73ac1f2ad65191baee1943cc938d7e` supplies the input manifest and preregistration tag. Provide that checkout, its specified data, and the rebuilt response and label directories to rerun the fits. The scripts accept paths on the command line.

```sh
python join_impact.py --engine ENGINE --rebuild REBUILD --output CHECKS
python direct_primary.py --engine ENGINE --rebuild REBUILD --output PRIMARY
```

For the full frozen-plan run, the engine uses `src/G/run_g.py` in real-data mode, the `partB-prereg-v1` tag and `analysis_plan/G_plan.md`, with the matching SHA-256 supplied in `P59_UNBLIND_PLAN_SHA256`. Retain the scratch results if interrupted and restart that stage into an empty output directory. The permutation count for the main contrasts is 10,000.

## Primary estimates from the fresh join

The 10,000-permutation P1 analysis against the fresh expression, marker and covariate tables is complete. Compare the measured panel and vehicle p values in `direct_H1.tsv` and pooled effects in `direct_H3.tsv`. The independent comparisons found zero differences in the H1, H3 and 232 per-compound rows.

## Idarubicin across doses

The [three dosing conditions](X2all_idarubicin_conditions.tsv) reproduced cell for cell. For the Y term, the exploratory response-size test gave p = 0.108 at 0.08 µM for 6 h (39 lines) and p = 0.0081 at 0.12 µM for 24 h (41 lines). The latter condition had already been examined in the original study. At 12 µM for 24 h, only two eligible lines expressed the Y marker, below the minimum of three. Dose, time and panel differ across conditions, so these estimates do not establish a dose–response relationship.

## Figure 4

The rebuilt line table and expression tensors reproduced [all 30,000 permutation draws](figure4_nulls.tsv.gz) in panels A–C and [all 12 effect estimates and confidence intervals](figure4_intervals.tsv) in panel D. The figure compares Y and XIST with TP53 mutation status on the same primary drug panel.

## Calibration and precision

The [nine calibration analyses](G_calibration_summary.tsv) reproduced all [1,200 pseudo-label rows](G_calibration.tsv) and all [46 historical-control rows](G_calibration_replay.tsv) exactly. The [two positive controls](G_gates.tsv) and the [precision estimates](G_E3_precision_labelfree.tsv) also matched.

## Vehicle referencing

In the primary panel, the exploratory XIST test gave p = 0.376 when each treated well was compared with vehicle wells on the same plate. Subtracting a single vehicle signature per line yielded p = 0.0090, with 20 nominally significant compounds and none at q < 0.1. The normalization changes the panel result. All [four contrasts](X1_H1.tsv), 460 compound rows, and 1,200 calibration rows matched the original files.

## Registered model and Figure 4 files

All [110 captured outputs](capture_files.tsv) matched byte for byte. Of the [35 model files](model_files.tsv), 32 were byte-identical, the gene table differed only in its gzip header, and two metadata tables recorded the new file paths and run time. The [four published tables](published_table_checks.json), [274 plotted values and annotations](fig4_partB_plotted.tsv), and [Figure 4 PDF](../../../manuscript/figures/fig4_partB.pdf) were byte-identical.
