# LINCS drug response

Across 41 cell lines and 114 testable compounds, neither Y-gene expression (p = 0.723) nor XIST expression (p = 0.113) showed a panel-wide association with drug response in the primary analysis. No compound met the registered criteria. The pooled Y estimate was −11.0% of mean response (90% CI −29.1% to 7.1%). The XIST estimate was 1.6% (−16.1% to 19.3%). The intervals leave room for meaningful differences. [Panel tests](part_b_h1.tsv) · [Pooled estimates](part_b_h3.tsv) · [Compounds](part_b_compounds.tsv)

## Vehicle controls and dose

XIST results depended on how vehicle controls were used. Comparing each treated well with vehicle wells on the same plate gave p = 0.376 in the primary panel. Subtracting a single vehicle signature per line gave p = 0.0090, with 20 nominally significant compounds and none at q < 0.1. At 0.12 µM, the corresponding panel results were p = 0.267 and p = 0.0173. Ponatinib and trametinib each reached q = 0.0749 in the exploratory lower-dose, same-plate analysis. At 10 µM for 6 h, all four Y and XIST panel tests gave p ≥ 0.40. [Primary panel](X1_H1.tsv) · [Other panels](X1all_H1.tsv) · [Vehicle diagnostic](X1_common_vector_summary.tsv)

For idarubicin, the exploratory Y response-size test gave p = 0.108 at 0.08 µM for 6 h (39 lines) and p = 0.0081 at 0.12 µM for 24 h (41 lines). The second condition had been examined previously. At 12 µM for 24 h, only two eligible lines expressed the Y marker, below the minimum of three. Both time and panel vary across the dose comparisons. [Three-dose results](X2all_idarubicin_conditions.tsv)

## From assay data to Figure 4

The rebuild read 30,056 drug and 4,330 vehicle signatures from the provider’s LINCS2020 Level 5 matrices, and 61,652 treated wells and 9,108 vehicle wells from Level 3. All 749 plates supplied enough vehicle wells for same-plate subtraction. All 45 rebuilt input files, including response arrays, axes and identifiers, matched byte for byte. [Source inventory](source_annotations.tsv) · [Input comparisons](rebuild_checks.tsv)

| Check | Result |
|---|---|
| Measured states | All 129 line labels and eligibility calls agreed. In the cutoff table, density estimates and ratios differed by at most 1.11 × 10⁻¹⁵. All cutoffs agreed. |
| Covariates | Six reloaded P2 ploidy values differed by at most 4.45 × 10⁻¹⁶. All 30 designs matched exactly. [Design comparisons](model_design_comparison.tsv) |
| Ploidy comparison | Panel p values were 0.50445, 0.67293 and 0.07929 across P1, P1b and P2, each from 10,000 within-lineage permutations. [Ploidy results](ploidy_control_comparison.json) |
| Controls and calibration | Both positive controls, 1,200 calibration assignments, 46 prior-control comparisons and three precision estimates reproduced. [Positive controls](part_b_controls.tsv) · [Calibration](G_calibration_summary.tsv) · [Precision](G_E3_precision_labelfree.tsv) |
| Transcription-factor sets | All 10,319 CollecTRI edges and 43 candidate transcription factors matched. [Set comparisons](regulon_checks.tsv) |

Every one of the 30,000 permutation draws behind Figure 4A–C reproduced. The tail counts were 7,231 for Y, 1,126 for XIST and 1,423 for the post hoc TP53 comparison. All 12 effects and their 90% and 95% intervals in Figure 4D matched. The regenerated PDF, PNG and 274-row plotting table were byte-identical. [Figure 4](../../../manuscript/figures/fig4_partB.pdf) · [Draws](figure4_nulls.tsv.gz) · [Intervals](figure4_intervals.tsv) · [Rendering check](figure4_check.json)

## Manuscript and supplement

All 336 stored Part B values and their printed forms agreed, including all 241 keys used in the text and tables. All 156 consistency checks passed. Dates and historical checkpoints come from the original plan and records. [Value-by-value comparison](manuscript_values.tsv) · [Manuscript checks](manuscript_values_check.json)

All 110 capture files, 18 exploratory files and four public result tables matched byte for byte. Of 35 registered model files, 32 were byte-identical. The compressed gene table had identical uncompressed contents. Two metadata tables recorded 15 new paths and the elapsed run time. [Model files](model_files.tsv) · [Captured outputs](capture_files.tsv) · [Exploratory files](exploratory_files.tsv) · [Public tables](published_table_checks.json)

## Reproduction

Run `python3 verify.py` in this directory to check the package. The 15 [public model modules](../../../analysis/lincs/src/G/) match the frozen code. The [public configuration](../../../analysis/lincs/config/) uses the same numerical settings. The complete checkout at commit `87640df15a73ac1f2ad65191baee1943cc938d7e` supplies the source manifest and the `partB-prereg-v1` tag.

With that checkout and the rebuilt expression and line tables, the commands below run the independent primary analysis and Figure 4 checks. `CAPTURE` is the original directory of retained permutation arrays.

```sh
python join_impact.py --engine ENGINE --rebuild REBUILD --output CHECKS
python direct_primary.py --engine ENGINE --rebuild REBUILD --output PRIMARY
python figure4_nulls.py --engine ENGINE --rebuild REBUILD --archive-capture CAPTURE --output NULLS
python figure4_intervals.py --engine ENGINE --rebuild REBUILD --output INTERVALS
```

The full analysis runs `src/G/run_g.py` in real-data mode with `analysis_plan/G_plan.md` and its SHA-256 in `P59_UNBLIND_PLAN_SHA256`.
