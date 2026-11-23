# Sensitivity of the LINCS drug-response design

These post hoc analyses ask what effects the primary panel could detect. They use the 41 eligible lines, their observed expression tensors and the original counts of marker groups within each lineage. Figure S3 reports performance under the specified simulated effects. It cannot establish the frequency of such effects in cell lines.

| Analysis | Design | Output |
|---|---|---|
| 1. Calibration | Count rejections in the registered pseudo-label tests and compare observed rates with nominal thresholds. Repeat 200 assignments for each marker under two alternative designs. | [Observed rejection rates](full_run/size_panel_test.tsv) · [Design comparisons](full_run/size_diagnostic.tsv) |
| 2. Planted effects and interval coverage | For each of 1,000 within-lineage assignments, fit the Y and XIST terms with 2,000 permutations. Compare unchanged expression with a 10%, 20%, 30% or 50% increase in response magnitude. In sparse scenarios, alter 1, 6 or 23 eligible compounds, either by a 50% response increase or a one-residual-SD shift across 49 of the 978 landmark genes. In two additional scenarios one or two lines retain their original response under a 20% increase. | [Detection by scenario](full_run/power_summary.tsv) · [Per-compound coverage](full_run/coverage_e2.tsv) · [Pooled coverage](full_run/coverage_e3.tsv) |
| 3. TP53 comparison | Fit the comparison exposure and examine AMG-232 in the positive control and factorial line sets. | [Method and outputs](scripts/tp53_yardstick.py) |
| 4. Marker classification | Refit the factorial analysis after omitting the three flagged lines together and individually. Test alternative states, including the mixed A549 observation. | [Affected lines](sensitivity_lines.tsv) · [Analysis](scripts/misclass_e3.py) |
| 5. Response scale | Estimate the pooled response-size difference and heteroskedasticity-robust intervals across panels and scales. Check the influence of weak-expression templates. | [Panel estimates](scripts/e3_table.py) · [Weak-template sensitivity](scripts/e3_weak_templates.py) |
| 6. Detectable effects | Use observed interval widths to calculate the effect magnitudes corresponding to 80% power under the stated approximation. | [Calculation](scripts/detectable_effects.py) |

The simulation uses seed 20261131 for pseudo-label assignments and 20261132 for planted scenarios. The implementation gives the per-assignment and per-term seed derivations. [Simulation code](scripts/power.py) · [Definitions](scripts/common.py) · [Run and verify](README.md)
