# Verification tables

The tables report the November 2026 reruns behind the Project 59 revisions and a replay of the earlier CMap/GTEx screen. They contain study-level counts and test statistics, plus comparison labels for the earlier screen. [Manuscript](../manuscript/manuscript.pdf) · [Supplement](../manuscript/supplement.pdf)

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

To typeset the [manuscript](../manuscript/manuscript.tex) and [supplement](../manuscript/supplement.tex), run `tectonic manuscript.tex` and `tectonic supplement.tex` from the `manuscript/` directory.

## What was rerun

| Analysis | Starting material and comparison |
| --- | --- |
| Phases 1–2 | All 2,641 vorinostat signatures were reread from the CMap level-5 matrix. The 12,328 two-sided Mann–Whitney tests recovered the same 6,877 adjusted hits and the same 20 genes above the effect threshold. All 572 signatures for AS-605240 were replayed in Phase 2. The other 49 drug classifications were checked against the archived table. [AS-605240 sensitivity to the number of tissues tested](../Phase2_Fishers_50Drugs/as605240_tissue_tests.csv) |
| Phases 3–5 | The 50 Phase 3 classifications, both 9,260-row Phase 4 exports and the 39,996 Phase 5 comparisons were checked as saved. XIST’s positive sign was checked in all 44 tissues in the GTEx v8 source. Fresh enrichment runs for these phases remain outside this verification. |
| Phases 6–7 | Replaying CMap level-5 scores for 36,523 signatures across 50 drugs matched all 39,996 released replay categories and skip reasons. Separately, the original October comparison table reproduced all 36,523 archived signature assignments. The underlying signatures reproduced the Phase 7 counts by sex and the per-line concentration and effective sample size for all eleven own-tissue drug–lineage pairs. |
| Phase 8 | The 18-row FAERS summary comes from the executed October notebook. The rank tests were recalculated from the printed per-drug counts. The unversioned live API may return different counts today. [FAERS counts and interpretation](../Phase8_SideEffects_FAERS/) |
| Manuscript and drug-treated follow-up | The Part A/B checks below read the released tables. Separate follow-up reruns covered the complete 14-plate Tahoe index, the two-plate expression matrix, archived external RNA and DNA measurements and raw intensity arrays. [Follow-up methods and source checks](../analysis/marker_stability/) |

## Run the checks

From the repository root, run this with Python 3.

```sh
python3 verification/verify_results.py
```

The verifier runs offline using the Python standard library. It checks file hashes, Part A counts, Part B ranges and denominators, the reported Benjamini–Hochberg correction and all 39,996 Phase 6 replay labels. For the historical phases it recalculates the adjusted values for 12,328 pilot genes, checks the Phase 2–5 classifications, reconciles the Phase 6–7 signature counts and reproduces the FAERS rank tests from the 18 [per-drug records](../Phase8_SideEffects_FAERS/per_drug_events.csv). These checks start from the exported tables.

A full analysis replay needs the study inputs and workflow. For Part A, start from the Tahoe-100M and EmeraldBay vehicle cells and the external RNA-seq comparison described in the manuscript. For Part B, use the LINCS2020 response panel and join the measured-state assignments before reproducing the within-lineage permutation test and response-size analysis. The archived [Phase 6 notebook](../Phase6_NoSampleFilter/notebook/phase6_no_sample_filter.ipynb) and [input notes](../Phase6_NoSampleFilter/sources/README.md) document the earlier replay’s matrix, metadata, gene sets and GSEA settings. These files support checks of the released tables. Full analysis requires the source input data and pipeline.
