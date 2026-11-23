# Project 59 in the drug-treated Tahoe cultures

## COLO 205 across independent cultures

Two external RNA studies and multiple DNA platforms show stronger signal from the measured long arm of chromosome Y in COLO 205. The independently cultured Broad sample gave a long-minus-short-arm difference of +1.395 log₂ units. HCT-15 and HT-29 gave +0.054 and −0.073 on the same platform. An older NCI karyotype recorded no identifiable Y in nine metaphases. [Results, figures and source records](independent_colo205/README.md) · [Broad raw DNA and two 450K arrays](independent_colo205/nci_dna/README.md)


[Integrated manuscript](paper.md) · [Integrated supplement](paper_supplement.md) · [Manuscript and supplement (PDF)](../../manuscript/followups/marker_stability/marker_stability.pdf)

For the marker-transfer results, the [cell-level analysis](marker_transfer/analyze.py) and [independent verification](marker_transfer/independent_check.md) give the treated-well comparisons and the COLO 205 threshold sensitivity. The [autosomal expression study](genomics_extension/README.md) tested drug responses across 948 genes on the two repeat plates. Matched conditions gave a larger correlation than mismatched conditions in every eligible line (24/24). Sharing a baseline marker group across distinct lines added 0.0109 correlation units (exact p = 0.204). [Source matrix and gene list](genomics_extension/) · [Data figure](genomics_extension/landmark_responses.png). The [recovery analysis](paper_supplement.md) measures the change in the relative cell count.

A separate [donor-label comparison](donor_response/protocol.md) tests whether differences between male- and female-labelled lines repeat for the same drug and dose. Eight bowel and eleven lung lines qualified. Across 11,550 whole-line reassignments the observed cross-plate score was near the permutation mean (one-sided p = 0.571). [Full donor result](donor_response/results.json) · [Permutation figure](donor_response/donor_response.png)

Across 379 named compounds, lines expressing Y genes or XIST showed almost the same average change in relative recovery between 0.05 and 5 µM as lines expressing neither marker (−0.0051 log₂ units, 95% model interval −0.0673 to 0.0570). Bowel, lung and skin contributed 26 eligible lines. Most concentrations occupied different plates. Concentration and plate effects therefore remain difficult to distinguish. The follow-up measures relative cell counts, a different outcome from the original paper's expression response to drug treatment.

[Vector figure](figure_1.pdf) · [Analysis plan](analysis_plan.md) · [Plate 13 addendum](plan_addendum.md) · [Independent count verification](results/independent_verification.json)

## Reproduce the experiment

From the repository root, run `python3 reproduce.py follow-up --work ../p59-run`. The command regenerates the recovery, marker agreement, landmark-expression and donor-label results and their figures under `../p59-run/marker_stability/`. It also recounts the recovery experiment from individual cell annotations and reruns the independent 249,480-assignment check of landmark expression.

The 8.6 MB `inputs/` directory contains the original fourteen well–line–marker-stratum count files, completion records, archived sample and driver metadata, the frozen expression-state table, and the relevant entries from the stream's SHA-256 manifest. Source data are pinned to [Tahoe-100M release 349e2a5](https://huggingface.co/datasets/tahoebio/Tahoe-100M/tree/349e2a55f5538a3e73b32a06c6bdd8e9e87011e8). Project 59's baseline marker rules and paper are at the [fellowship repository](https://github.com/MateenahJAHAN/Project-59-CMAP-Sex-Differences/tree/d9a94ddedeb65ced19c16f6cb78497d837f161b0).

From this directory with Python 3.12, NumPy 2.5.4, pandas 3.1.1, PyArrow 26.0.1, SciPy 1.18.1 and Matplotlib 3.11.2, run

```sh
python analyse.py --root inputs
python figure.py
```

The recovery check obtains the [Tahoe-100M observation file](https://huggingface.co/datasets/tahoebio/Tahoe-100M/blob/349e2a55f5538a3e73b32a06c6bdd8e9e87011e8/metadata/obs_metadata.parquet) (2,293,981,573 bytes) and verifies its [SHA-256](../../verification/source_replay/part_a/receipts/tahoe_observations.json). It reuses the file when Part A has already downloaded it under `../p59-run/tahoe/index/meta/`. It reads the filter, cell line and well for all 100,648,790 profiles. Summing the 95,624,334 full-filter profiles reconstructs a 1,344-well by 50-line matrix. Every entry, including zeros, is compared with the fourteen released plate tables. The recovery contrast uses the released vehicle classifications. When the work directory contains fresh Part A calls, the verifier also checks the line accessions, tissues and marker assignments against those calls.

The full count matrix is written to `marker_stability/results/independent_well_line_counts.tsv.gz`. The adjacent `independent_verification.json` gives the number of discordant entries, the estimated contrast and the two-sided permutation probability across all 972,972 within-tissue assignments. To repeat the recovery check by itself, use the atlas environment created by the top-level runner.

```sh
../p59-run/envs/atlas/bin/python analysis/marker_stability/verify.py --work ../p59-run
```

`results/lines.tsv` identifies all 50 lines and why each contributes or falls outside the primary contrast. Full per-well counts by group are in the input tables. The statistics for all 379 compounds, observed plate repeats, controls and individual-line exclusions are in [results](results/). This is an independent computational follow-up to the fellowship paper. Codex designed and ran the analyses and prepared the manuscript for review. The fellowship manuscript’s authorship applies to the preceding study.

## Data sources

Tahoe Bio released this version of [Tahoe-100M](https://huggingface.co/datasets/tahoebio/Tahoe-100M/tree/349e2a55f5538a3e73b32a06c6bdd8e9e87011e8) under CC0 1.0. The original Y/XIST classifications and Cellosaurus donor annotations carry the [fellowship table credits and CC BY 4.0 terms](../data/README.md). This analysis derives well-level counts and compares the archived marker rules. The external cell-line assays are credited by accession and sample in the [source records](independent_colo205/README.md). The raw GEO, NCI and CCLE downloads remain with the respective providers. The [file manifest](independent_colo205/sources.tsv) records URLs, sizes and SHA-256 digests for the external RNA, array CGH and GDSC analysis. The [450K and SNP-array records](independent_colo205/nci_dna/README.md) identify their separate source archives and probe definitions.
