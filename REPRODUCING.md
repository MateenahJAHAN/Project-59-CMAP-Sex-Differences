# Reproducing Project 59

The paper has two analyses. Part A measures sex-chromosome expression and checks those measurements against independent assays. Part B tests the measured markers against drug response in LINCS2020. The numbered fellowship phases precede the paper. Their original notebooks, results and replay differences are preserved in the [phase guide](fellowship/README.md).

## Start from GitHub

On macOS or Linux, use Git and Python 3.11 or later. The source analyses also need [uv](https://docs.astral.sh/uv/getting-started/installation/) and an internet connection to the data providers. Install the [Tectonic command-line program](https://tectonic-typesetting.github.io/en-US/install.html) to build the PDFs. The launcher uses uv to install separate Python environments for the exploratory screens, the atlas measurements, the registered permutations and the figures. Set the work directory on a volume with room for downloaded inputs and intermediate results.

```sh
git clone https://github.com/MateenahJAHAN/Project-59-CMAP-Sex-Differences.git
cd Project-59-CMAP-Sex-Differences
python3 reproduce.py verify
python3 reproduce.py all --work ../p59-run --workers 4
```

`verify` runs offline. It rechecks the phase summaries, Part A and B evidence, and all seven released figures. It recalculates the 66 points in Figure S3 from the archived assignment records. The `all` command starts the numbered fellowship analyses, runs both manuscript parts, repeats the planted-effect experiment, draws the figures from those outputs and compiles the two PDFs. The resulting paper and supplement are in `../p59-run/paper-from-run/`. Keep the work directory between stages.

## Run just the analysis you need

Run each command from the repository root. The same `--work` path lets one stage consume another’s results.

| Command after `python3 reproduce.py` | Starting material | Results under the work directory |
|---|---|---|
| `paper --work ../p59-run` | Released numerical tables and TeX sources | `paper/manuscript.pdf`, `paper/supplement.pdf` and seven newly drawn figure PDFs |
| `fellowship --work ../p59-run` | Broad L1000, GTEx and the archived Phase 8 drug counts | `fellowship-results/`, including the eight phase outputs and comparisons with the archived results |
| `part-b --work ../p59-run` | Broad signature and well matrices, bulk RNA and covariate sources | `results/G/` for the registered panel tests, controls and calibration, `capture/` for the permutation arrays, and `results/figure4/` for its draws and effect intervals |
| `part-a --work ../p59-run` | Tahoe-100M, EmeraldBay, sci-Plex, independent bulk RNA and cell-line annotations | `results/tahoe/latest/`, `results/F1/`, `results/F2/`, `results/str/`, `results/depth_matching/` and the complete `results/supplement/supplementary_table_S1.tsv` |
| `transfer --work ../p59-run` | The primary Part B design and the independent Part A assays | `results/transfer/` for the 41-line stock comparison in Note S16 and its sensitivity exclusions |
| `claims --work ../p59-run` | Completed Parts A and B, including the locally reconstructed Table S1 | `results/claims/` for the 96 claim comparisons, the line roster and coverage |
| `l1000 --work ../p59-run` | Part B bulk RNA states and 2,800 control wells from the provider matrix | `results/l1000/` for the independent L1000 marker-score diagnostic |
| `power --work ../p59-run --workers 4` | The completed Part B tensors, model design, calibration and permutation captures | `power/raw/` for 1,000 assignment checkpoints and `power/results/` for detection rates, coverage and diagnostics |
| `paper-from-run --work ../p59-run` | Completed Parts A, B and planted-effect simulation | `plot-inputs/`, seven redrawn figures and both PDFs in `paper-from-run/` |
| `follow-up --work ../p59-run` | Published treated-well and repeat-plate counts from the later Tahoe study | `marker_stability/`, including regenerated marker, recovery, landmark-expression and donor-label analyses |

The follow-up is a separate study. Its [source guide](analysis/marker_stability/README.md) gives the additional commands for the 100-million-profile check, COLO 205 RNA/DNA assays and raw array analyses.

## Downloads and runtime

Allow several hours for the models and permutations. The work can take longer with fewer cores or on a slow connection. As a planning allowance, use a volume with at least 100 GB free for the source workflow; this is headroom, not a measured minimum. The datasets have different access patterns.

- The fellowship runner fetches selected byte ranges from the Broad level-5 matrix. The assembled matrix contains 36,523 signatures × 12,328 genes and occupies 1.8 GB. Phases 2–5 use GSEApy 1.2.1. The Phase 6 replay uses 1.3.2 in its own library directory. [Fellowship setup and individual phases](fellowship/README.md)
- The Tahoe vehicle reduction read 14.29 GB of expression byte ranges across 1,910 source shards. Its observation table is 2.29 GB. The independent EmeraldBay reconstruction reads 116 shards totalling 57.71 GB, and the sci-Plex source files total 3.70 GB. EmeraldBay shards can be discarded after reduction. [Atlas acquisition and callers](analysis/part_a/) · [Source sizes and hashes](verification/source_replay/part_a/receipts/source_verification.json)
- The LINCS model reads both signature-level responses and same-plate vehicle wells. The [source manifest](analysis/lincs/config/source_files.tsv) identifies the exact versions of its additional expression and metadata files. Some providers reuse a filename for new revisions. A different hash stops the run. Use the recorded version of that file in `<work>/sources/`. [Part B model and inputs](analysis/lincs/)
- The power study fits 2,000 permutations per test for each of 1,000 within-lineage assignments. Completed assignments are checkpointed. To resume after interruption, repeat the power command with the same directory and environment. [Power protocol, results and offline recount](verification/source_replay/power/)

## Reading a difference

Phase 2 reproduces both tissue-correction denominators. In Phase 3, two compounds cross the enrichment threshold when tied genes change order. In Phase 4, 54 threshold calls appear only in the archive and 35 only in the replay. In Phase 6, nine B/C tissue classifications changed. Eight of the 36,523 individual signature assignments moved from C to B, all for idarubicin in stomach lines. All 2,656 own-tissue assignments agreed. The comparison files identify the affected compounds or tissues. [Phases 2–5 comparison](verification/source_replay/phases2_5/RESULTS.md) · [Phase 6 differences](Phase6_NoSampleFilter/results/replay_category_changes.csv)

For Phase 8, only the original per-drug counts survive. They support the drug-level comparison for all 18 drugs. The original reaction-by-sex records and API responses needed for a raw historical replay were not archived. Querying the live FDA service today would create a new dataset. [Archived counts and notebook](Phase8_SideEffects_FAERS/)

The historical reproduction and the present-day build answer different questions. The [source replay record](verification/source_replay/) gives the original comparisons, including discrepancies, with the published results. The [paper source](manuscript/README.md) identifies the data behind each figure.
