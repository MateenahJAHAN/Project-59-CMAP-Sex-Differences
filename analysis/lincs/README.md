# Drug response in LINCS cell lines

Part B asks whether Y-gene expression or XIST predicts drug-induced changes in transcription across cell lines of the same lineage. The primary panel contains 41 eligible lines and 114 tested compounds. Equation inputs, permutation fits, model tables and figure data can all be rebuilt locally.

## Run from provider data

Use Python 3.13 on Linux or macOS. Allow several hours and about 12 GB of free working space. The workflow reads the required rows from the public LINCS matrices over HTTP. The level-3 and matched level-5 well extraction alone transfers 7.55 GB. The whole matrices are hundreds of gigabytes and stay at the provider.

From the repository root:

```sh
python3.13 -m venv .venv-lincs
.venv-lincs/bin/python -m pip install -r analysis/lincs/requirements.txt
.venv-lincs/bin/python analysis/lincs/reproduce.py --work ../p59-lincs-run --stage all --workers 4
```

The work directory must be outside the repository. Completed calibrations are retained there under a fingerprint of the expression tensor, samples, compounds, design, permutation keys and numerical engine. A later invocation can reuse those completed assignments. Smaller worker counts reduce concurrency without changing the registered settings.

| Stage | What it reconstructs |
| --- | --- |
| `download` | Fourteen checksum-pinned source files from CLUE, DepMap, Cell Model Passports, Cellosaurus, OmniPath, Ensembl and PubMed Central. The `labels-c` step also downloads the GDSC workbook. |
| `labels-a`, `labels-c` | Bulk RNA marker states and independently matched cell-line annotations. |
| `level5` | Three landmark-expression panels, their compound and cell-line axes, and their vehicle controls from LINCS2020 level 5. |
| `wells` | Level-3 wells, level-5 signatures, per-plate DMSO centering, landmark-expression matrices and counts. |
| `design` | Joins the rebuilt state tables, TP53 and cell-line covariates. Reconstructs all 13 pseudo-label banks and checks line membership, activity stratification and tensor checksums. |
| `model` | Fits the registered contrasts after the controls and calibration pass. Writes `results/G`, its five figures, and the fit and null arrays in `capture`. |
| `exploratory` | Repeats the vehicle and idarubicin follow-ups in `results/G_exploratory`. |
| `figure4` | Independently refits 30,000 panel-null draws and 12 effect intervals and renders the manuscript PDF and PNG in `results/figure4`. |
| `check` | Compares the four reported result tables, model outputs, capture arrays, exploratory tables and figure data against the previous full source replay. Writes `results/public_tables/full_comparison.json`. |

Run a single stage with the same `--work` path after completing its dependencies. Known source revisions can also be placed at `WORK/sources/<filename>`. The [source manifest](config/source_files.tsv) identifies each by URL, byte length and SHA-256. A rolling provider release with different bytes stops the replay.

## Inputs fixed before fitting

[Design files](design/) include the original pseudo-label assignments, calibration benchmarks and pre-fit control summaries. The source workflow regenerates every assignment from its seed and verifies the uncompressed bytes before fitting the model. The model recomputes the calibration statistics and compares them with the original benchmarks. The three auxiliary tables in [inputs](inputs/) contain chromosome assignments, historical compound overlap and the 30 highest-weight landmark targets for GTEx lymphoblastoid PANDA transcription factors with at least ten prior motif targets. [PANDA source data](https://zenodo.org/records/838734)

The normal run rebuilds the chromosome and historical overlap tables. The model checks the released PANDA reference against its frozen fingerprint.

The public entry point records a post-publication replay. The original registered plan used tag `partB-prereg-v1`, commit `572b97498d8d555e6f3616c99cb3289351223c77`, and plan SHA-256 `32bdc666e9c85a43368a8f6f70269155b24fb0fde0095326c5d2178588a8bd1f`. Numerical fits, calibration thresholds, seeds and permutation counts follow that analysis. The public runner checks the [released code and configuration](src/G/SHA256SUMS.txt) and every opened model input. Seven ploidy values can differ from their frozen text by less than a trillionth. The alternative byte fingerprints are fixed and the lineage-specific median splits were identical. [Ploidy comparison](../../verification/source_replay/part_b/ploidy_control_comparison.json)

[Manuscript](../../manuscript/) · [Part B reproduction evidence](../../verification/source_replay/part_b/)
