# Figure S3 power and calibration

On the primary drug panel, effects were planted into measured expression profiles and Y/XIST labels reassigned within lineage. Each of 1,000 assignments used 2,000 permutations. These post hoc detection probabilities describe the 41-line design and the particular effect planted. They do not estimate the prevalence of biological effects.

| Term | Primary Level 3 detection at a 20% change | Null panel rejections, p < 0.025 | 95% exact interval |
|---|---:|---:|---:|
| Y marker | 34.0% | 3/1000 | 0.06%–0.87% |
| XIST | 40.9% | 5/1000 | 0.16%–1.16% |

The 20% figures are interpolations between measured simulation points. The null counts measure the empirical rejection rate when the reassigned markers carry no effect. The [paired design diagnostic](size_diagnostic.tsv) uses the first 200 assignments. The [configured calibration](size_panel_test.tsv) summarizes a separate set of 1,200 assignment–contrast records.

[Figure S3](figS3_power_values.tsv) contains 66 of 66 coordinates with all assignments usable, across 18 of 18 fully observed curves. The [denominators](figure_denominators.tsv) identify the usable assignments and eligible compound range for each coordinate. Sparse effects were planted only among eligible compounds. Fifty-four binary detection proportions use Wilson intervals. The twelve orange points average the fraction of planted compounds detected per assignment and use normal intervals from the variation among those fractions. The [coverage tables for pooled effects](coverage_e3.tsv) and [per-compound effects](coverage_e2.tsv) separately evaluate the 90% and 95% study intervals.

The [comparison](comparison.json) checks the raw simulation records, every E2 interval and all scientific fields of the revision's tables against the archived study. The [independent recount](receipt.json) derives figure coordinates, intervals, coverage and calibration from replicate records. Source tables and per-line expression values remain in the local study workspace.
