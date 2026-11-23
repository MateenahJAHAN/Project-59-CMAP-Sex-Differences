# Verification tables

The tables report the November 2026 reruns behind the Project 59 revisions and a replay of the earlier CMap/GTEx screen. They contain study-level counts and test statistics, plus comparison labels for the earlier screen. [Manuscript](../manuscript/manuscript.pdf) · [Supplement](../manuscript/supplement.pdf)

The [source-data replay](source_replay/) includes the complete Phase 2–5 analyses, three revised historical figures and tests of the manuscript drug responses from rebuilt inputs.

## Read the results

| Analysis | What the tables show |
| --- | --- |
| [Part A](results/part_a.tsv) | The Tahoe-100M caller used 2,175,916 of 2,205,786 vehicle cells and called 50 lines. The measured expression patterns agreed with external RNA-seq in **47 of 49** lines. Agreement between the states that use donor annotations was **46 of 49**. EmeraldBay agreement was 40 of 41 among lines with high or medium confidence in both datasets, and 42 of 47 across all shared lines. The separate post hoc Amelogenin check agreed with 48 of 50 expression-based Y calls. |
| [Part B panel](results/part_b_h1.tsv) · [Response sizes](results/part_b_h3.tsv) | The primary comparison covered 41 LINCS2020 lines across 11 lineages and 114 tested compounds at 1 µM and 24 hours. The panel p values were 0.7231 for Y and 0.1127 for XIST. The XIST vehicle-control p value was 0.0323. Response-size intervals were wide, especially for Y. |
| [Part B compounds](results/part_b_compounds.tsv) · [Controls](results/part_b_controls.tsv) | The rankings include 114 tested compounds per term. The table also shows the untested AZD-8055 entry and the separate DMSO vehicle result. No ranked compound passed the full decision rules. AMG-232 ranked first in the TP53 positive controls. |
| [Earlier Phase 6 replay](results/phase6_counts.tsv) · [Comparison labels](results/phase6_comparison_labels.tsv.gz) | The replay assigned categories to 39,996 drug–lineage–tissue comparisons. A is enrichment in a tissue mapped to the lineage, B is enrichment in another tissue, C is a test without significant enrichment and D is no test outcome. Seventeen compounds reached A or B, including six that reached A. |

Part A expression calls represent how the cells expressed the markers in the cultures measured. The donor-annotation-dependent names use the recorded label for some categories. The label-independent comparison collapses these names to the measured patterns. [Methods](../manuscript/manuscript.pdf)

The Part B analysis plan was fixed internally before the expression calls were joined to drug response. The panel test has limited power for modest or compound-specific effects. The XIST vehicle result limits interpretation of that contrast, and the pooled interval checks do not establish equivalence. The counts of lines with and without each model indicator are not counts of mutually exclusive expression groups. [Design and sensitivity checks](../manuscript/supplement.pdf)

The Phase 6 replay used GSEApy 1.3.2. Five comparisons moved from C to B and four from B to C. The other 39,987 labels agreed, including every A and D assignment. [Nine changed comparisons](../Phase6_NoSampleFilter/results/replay_category_changes.csv) These labels describe pooled comparisons with tissue gene sets. They do not estimate drug-response differences between donor-sex groups. [Recorded Phase 6 analysis](../Phase6_NoSampleFilter/README.md)

To rebuild all seven figures and typeset both PDFs in a separate folder, use the [paper build](../manuscript/README.md).

## What was rerun

| Analysis | Starting material and comparison |
| --- | --- |
| Phase 1 | Rereading all 2,641 vorinostat signatures from the CMap Level 5 matrix recovered 6,877 adjusted hits in 12,328 two-sided Mann–Whitney tests. The same 20 genes exceeded the effect threshold. |
| Phase 2 | Repeating 616,400 gene tests and 2,200 tissue tests reproduced all 50 compound classifications under the original correction. Including all 44 eligible tissues in the correction moved the AS-605240 liver result from 0.0247 to 0.0544. [All tissue tests](source_replay/phases2_5/comparison/phase2_all_tissue_tests.csv) |
| Phases 3–5 | Two Phase 3 compounds changed category when tied genes changed order. In Phase 4, 9,125 enrichments appeared in both runs, 54 appeared only in the archive and 35 appeared only in the replay. Phase 5 reproduced all 39,996 categories, including all 620 comparisons that yielded a test result. [Per-phase results and figures](source_replay/phases2_5/RESULTS.md) |
| Phases 6–7 | Replaying CMap level-5 scores for 36,523 signatures across 50 drugs matched all 39,996 released replay categories and skip reasons. Separately, the original October comparison table reproduced all 36,523 archived signature assignments. The underlying signatures reproduced the Phase 7 counts by sex and the per-line concentration and effective sample size for all eleven own-tissue drug–lineage pairs. |
| Phase 8 | The 18-row FAERS summary comes from the executed October notebook. The rank tests were recalculated from the printed per-drug counts. The unversioned live API may return different counts today. [FAERS counts and interpretation](../Phase8_SideEffects_FAERS/) |
| Manuscript | Two callers reproduced all 50 Tahoe expression states from vehicle cells. The registered LINCS analysis reproduced all four published drug-response tables using rebuilt inputs. [Cell-level replay](source_replay/part_a/) · [Drug-response replay](source_replay/part_b/) |
| Drug-treated follow-up | Separate reruns covered the 14-plate Tahoe index, the two-plate expression matrix, external RNA and DNA measurements and raw intensity arrays. [Methods and source checks](../analysis/marker_stability/) |

## Run the checks

From the repository root, run this with Python 3.

```sh
python3 verification/verify_results.py
```

The verifier runs offline using the Python standard library. It checks file hashes, Part A counts, Part B ranges and denominators, the reported Benjamini–Hochberg correction and all 39,996 Phase 6 replay labels. For the historical phases it recalculates the adjusted values for 12,328 pilot genes, checks the Phase 2–5 classifications, reconciles the Phase 6–7 signature counts and reproduces the FAERS rank tests from the 18 [per-drug records](../Phase8_SideEffects_FAERS/per_drug_events.csv). These checks start from the exported tables.

To repeat an analysis from provider data, start with the [source-data replay](source_replay/). It links the inputs, scripts and comparisons for each analysis. The [Phase 6 notebook](../Phase6_NoSampleFilter/notebook/phase6_no_sample_filter.ipynb) and [input notes](../Phase6_NoSampleFilter/sources/README.md) document the earlier pooled screen.
