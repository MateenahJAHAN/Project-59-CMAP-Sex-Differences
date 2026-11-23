# Fatima Fellowship

The fellowship began with a comparison of vorinostat signatures from male- and female-donor cell lines. Repeating the analysis across compounds exposed a sampling problem. The same cell line could contribute hundreds of profiles, and several apparent tissue enrichments came entirely from one donor group. The paper asks how well donor sex describes the sex-chromosome expression of the cells used in drug-response experiments. The paper measured those markers in Tahoe-100M and tested their relationship to response in LINCS2020. [Manuscript](../manuscript/manuscript.pdf) · [Supplement](../manuscript/supplement.pdf) · [Paper analyses](../analysis/)

## The eight phases

| Phase | What changed, and what the analysis found |
| --- | --- |
| [1 · Vorinostat](../Phase1_Vorinostat/) | Among 2,641 profiles from 116 lines, 20 genes passed both the adjusted *p* cutoff and an absolute difference in mean z-score above 1. Repeated profiles from a line contributed separately. |
| [2 · Fifty compounds](../Phase2_Fishers_50Drugs/) | Two AS-605240 genes overlapped the GTEx liver set. Its adjusted *p* rose from 0.0247 to 0.0544 when the correction included all 44 eligible tissues. No drug passed the latter test. |
| [3 · Within cell lineages](../Phase3_PerLineage_GSEA/) | Ranking genes within 86 drug–lineage pairs reduced mixing across tissue origins. Two of the 50 compounds changed category on replay. Reordering tied gene scores moved both marginal calls across the threshold. |
| [4 · Direction of enrichment](../Phase4_SplitLibraries_DirectionFix/) | Separating female- and male-biased GTEx sets clarified which reference genes were enriched in a pooled drug response. Across 53,064 tests, 54 archived hits no longer met the FDR cutoff and 35 other comparisons reached it. |
| [5 · A minimum number of lines](../Phase5_PerSignature_Categories/) | Requiring five distinct lines of each donor sex left 620 testable comparisons out of 39,996. One hit involved a reference tissue matched to the lineage, and 79 involved other tissues. |
| [6 · Removing that minimum](../Phase6_NoSampleFilter/) | This brought 6,122 comparisons into the test. Seventeen of 50 compounds reached enrichment in any tissue and six in a tissue mapped to the lineage. The rankings still pooled donor groups. |
| [7 · Counting independent cell lines](../Phase7_SexBalance_Audit/) | Among profiles assigned a matched-tissue hit, 2,426 of 2,656 were from female-donor lines. The group contained 17 distinct female-donor lines and 18 male-donor lines. Two of 11 pairs contained at least two lines of each donor sex. The most profiled line supplied 61 of 255 profiles for bortezomib in blood lines and 21 of 28 for YM-155 in skin. |
| [8 · Adverse-event reports](../Phase8_SideEffects_FAERS/) | A separate comparison of 18 drugs found median fractions of flagged reaction terms of 11.3% and 14.9% in drugs with and without a Phase 6 enrichment (*U* = 27, *p* = 0.3283). Report volume correlated with the fraction flagged (Spearman *ρ* = 0.853). |

## Run the fellowship analyses

Use Python 3.12 on Linux or macOS and run these commands from the repository root.

```sh
sh fellowship/install.sh "$PWD/../p59-env"
../p59-env/bin/python fellowship/run.py \
  --data "$PWD/../p59-inputs" --out "$PWD/../p59-results" --workers 2
```

The first run retrieves LINCS2020 Level 5 perturbation signatures and GTEx v8 tissue gene sets. It needs network access and several GB of free disk space. The downloader reads the 36,523 selected signatures from the Broad compound GCTX by HTTP byte range and verifies each signature and the finished 12,328-gene matrix. The assembled matrix occupies 1.8 GB. The signature selection, input hashes and gene mappings are available for inspection. [Input preparation](sources.py) · [Pinned input record](../verification/source_replay/phases2_5/source_manifest.json) · [Provider matrix verification](../verification/source_replay/lincs_inputs/)

Keep the input and output directories outside the checkout. Use `--prepare-only` to stage the inputs. Add `--phases 1 2` to choose phases. Phase 7 runs Phase 6 when its prerequisite is missing. The default two workers divide Phases 2–5 by compound. Interrupted jobs reuse completed compound files in the same output directory. Use an empty output directory for an independent repeat. GSEA uses seed 42. The install script puts GSEApy 1.2.1 and the Phase 6 replay’s 1.3.2 in separate import paths.

| Output under `--out` | Contents |
| --- | --- |
| `phase1/` | All 12,328 vorinostat gene tests and their numerical comparison with the archived export |
| `phase2_gene_tests/` | All 616,400 gene tests, with gene IDs, Ensembl mappings, mean differences, *p* and adjusted *p* |
| `outputs/phase2/` through `outputs/phase5/` | Compound-level files for all four screens |
| `comparison/` and `figures/` | Fresh calls beside the archives, and three figures rebuilt when Phases 3–5 are run together |
| `notebook_figures/` | Original pilot, fifty-compound overlap and full-screen plots rendered from the fresh phase outputs |
| `phase6/` | Every drug–lineage–tissue result, its FDR or reason for being untested, and any differences from the released replay |
| `phase7/` | Signature counts, distinct cell lines and concentration by pair, for the regenerated screen and the historical classification |
| `phase8/` | Drug-level comparison recalculated from the 18 archived openFDA counts |

The default run also rebuilds the four split-library Phase 6 follow-ups, including eight plots, under `phase6/cases/`. To rebuild those cases on their own, use the case-study script.

```sh
../p59-env/bin/python Phase6_NoSampleFilter/code/02_case_study_plots.py \
  --data "$PWD/../p59-inputs" --out "$PWD/../p59-cases"
```

The original FAERS analysis queried a changing openFDA index. The October 2026 drug-level table survives, but the full reaction-by-sex API responses are absent from the available archive. Phase 8 starts from those drug totals. It reproduces the drug-level comparison and its correlation with volume. Repeating the original queries today would measure a different snapshot. To run that phase alone, use `--phases 8`. It requires no LINCS or GTEx download. [Recorded FAERS counts](../Phase8_SideEffects_FAERS/per_drug_events.csv) · [FDA adverse-event data](https://open.fda.gov/apis/drug/event/)

The runner compares every resulting row in Phases 2–5 with the released replay. To repeat that comparison, run `../p59-env/bin/python fellowship/compare_release.py --out "$PWD/../p59-results"`. It joins by compound, lineage, tissue and, for Phase 4, reference-set direction. It reports numerical changes beyond a 1e−12 absolute or relative tolerance.

The original Phase 3 compound classifications, the Phase 4 threshold crossings and the Phase 6 differences between GSEA versions are retained beside the replays. Use the tables for the chosen run when quoting a number. [Phases 2–5 comparison](../verification/source_replay/phases2_5/RESULTS.md) · [Phase 6 changed calls](../Phase6_NoSampleFilter/results/replay_category_changes.csv)

An independent run from the downloaded matrix recovered every released row in Phases 2–5 across all 50 compounds (2,200, 3,784, 53,064 and 39,996 rows, respectively). Phase 6 reproduced all 39,996 categories and skip reasons. [Phase-by-phase comparison](../verification/source_replay/phases2_5/fellowship_run_comparison.json)
