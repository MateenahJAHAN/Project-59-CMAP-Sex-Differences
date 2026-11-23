# Reproduce the measured RNA states

Part A measures XIST and Y-linked RNA in vehicle-treated cell lines. Two independent programs read Tahoe-100M. EmeraldBay and sci-Plex provide additional cultures. Cellosaurus amelogenin records provide a separate DNA comparison.

## Run from source

From the checkout root on Linux or macOS, use Python 3.12. Choose a work directory on a volume with room for source downloads, intermediate counts, and checkpoints. A full run transfers more than 80 GB. EmeraldBay alone streams 57.7 GB in 116 shards and deletes its downloaded expression shards after processing. The uncompressed per-line sparse matrices from the second Tahoe reducer require additional local space. Allow several hours and use a machine with enough RAM for the full 95.6 million-row Tahoe index. The published workflow uses CPUs.

```sh
python3.12 -m venv .venv-atlas
.venv-atlas/bin/python -m pip install -r analysis/tahoe/requirements.txt
.venv-atlas/bin/python analysis/part_a/run.py --work /absolute/path/to/project59-work --stage all
```

The runner resolves scripts relative to the checkout and writes data and results under `--work`. Both manuscript parts can share the same work directory. The six bulk-source hashes and sizes are in [the bulk manifest](../../verification/source_replay/bulk_RNA/sources.json). Frozen snapshots already obtained from the provider may be placed under `WORK/sources/` using the manifest basenames. The Sanger filename `rnaseq_latest.csv.gz` and the Cellosaurus flat-file URL can change at the provider. The study used Cellosaurus release 57.0. The [Cellosaurus archive instructions](https://www.cellosaurus.org/faq) describe access to older releases through Yareta. Acquisition checks the file digest before accepting a cached or downloaded snapshot.

| Stage | Inputs and computation | Main output under the work directory |
|---|---|---|
| `labels-c` | CCLE and Sanger expression, model IDs, Cellosaurus and screening membership | `results/C/C_cellline_sex_labels.tsv` |
| `tahoe` | Download five Tahoe metadata tables. Read all 3,388 shard footers, select vehicle row groups and reduce original counts. Fit the first caller, sensitivity checks, ambient background, mixtures and identity tests. | `results/tahoe/latest/partA_local_states.tsv`, diagnostic plots under `figures/` in that directory |
| `second` | Read the Tahoe observation table and all shard footers with a second parser. Independently extract vehicle rows, fit states, ambient estimates, mixtures and identity controls. | `results/F2/F2_states.tsv`, `F2_mosaic.tsv`, `figures/F2/` |
| `external` | Reduce all EmeraldBay shards for identity controls and vehicle calls. Parse all three sci-Plex count matrices and their labels. | `results/F1/F1_emerald_states.tsv`, `F1_sciplex_states.tsv`, `figures/F1/` |
| `consensus` | Join the two Tahoe calls by cell-line accession and compare the consensus with EmeraldBay. | `results/consensus/tahoe_partA_states_consensus.tsv`, `mosaic_resolution.tsv`, `A2_emerald_vs_full_agreement.tsv` |
| `depth-matching` | Match depths across seven line–culture pairs using ten thinning draws per pair, four estimators and the [original analysis rule](secondary/depth_matching.md). | `results/depth_matching/A3_pairs.tsv`, `A3_replicates.tsv`, `A3_unmatched.tsv`, `A3_summary.json` |
| `table-s1` | Join the Tahoe consensus, EmeraldBay, three sci-Plex samples, LINCS, GDSC2 and PRISM. Reapply the Cellosaurus label to the nine-point LINCS sensitivity grids and attach the four cross-caller mosaic ranges. | `results/supplement/supplementary_table_S1.tsv`, 1,694 accession/source pairs from 1,123 cell lines |
| `str` | Parse the pinned Cellosaurus flat file for 1,102 accessions. Compare these DNA records with independently measured Y RNA in Tahoe, EmeraldBay, sci-Plex, CCLE and CMP. | `results/str/figure3_supp_table_amelogenin_by_line.tsv`, `figure3_agreement_by_source.tsv` |
| `transfer-s16` | Check the 41 Part B cell lines against the atlas, bulk RNA and STR calls. Requires both manuscript parts in the same work directory. | `results/transfer/supp_table_partB_transfer.tsv`, `sensitivity_lines.tsv` |

`all` executes the stages in dependency order. It also invokes the [bulk label A pipeline](../lincs/reproduce.py) for LINCS cell lines. A single stage may be rerun with `--stage NAME`. Use `tahoe-reduce`, `second-reduce` or `external-reduce` to stop after extracting the counts. Their corresponding `*-call` stages start from those work-directory reductions. Each completed first-caller run gets its own result directory. `results/tahoe/latest` points to the most recent successful run.

## Paper and supplement

| Output | Measurements and inputs |
|---|---|
| Main Figure 2, marker states | `results/consensus/tahoe_partA_states_consensus.tsv`, `results/tahoe/latest/partA_local_states.tsv` and `results/tahoe/latest/partA_local_gene_level_scores.tsv`. The primary caller also writes the cutoffs and label audit. |
| Main Figure 3, RNA and STR | The two `results/str/figure3_*.tsv` tables contain all accession-level calls, comparison types and per-source counts. |
| Supplementary Figure S1, ambient RNA | `results/tahoe/latest/partA_local_plate_background.tsv` and `partA_local_foreign_cells.tsv`. |
| Supplementary Figure S2 and Table S8, A549 | `results/tahoe/latest/partA_local_depth_curve.tsv`, `results/depth_matching/A3_pairs.tsv` and `A3_replicates.tsv`. |
| Supplementary mosaic sensitivity | `results/consensus/mosaic_resolution.tsv` joins the independent estimates. Two fitted endpoints have shown small numerical differences across software environments. The [software comparison](../../verification/source_replay/part_a/receipts/mosaic_version_comparison.json) records binary decisions and unrounded differences. |
| Supplementary Table S16, cross-stock comparison | After Parts A and B finish, run `python reproduce.py transfer --work /absolute/path/to/project59-work`. The [table comparison](references/s1_s16_comparison.json) records the replay of Tables S1 and S16 and the three lines flagged for the Part B sensitivity analysis. |

The historical pilot comparator is [50 Tahoe lines](references/tahoe_pilot_states.tsv). It enters agreement checks and supplies annotated sex for missing external annotations. The primary state estimates come from full vehicle counts. Cellosaurus STR derives from stocks characterized by the contributing laboratories. The methods cannot establish that these were the exact stocks sequenced by an atlas. The A549 measurements are regenerated from the atlas counts. The supplement cites the ATCC product sheet for its historical description of the cell line. The pilot's construction, historical decision records, and the original software-version sensitivity run have separate provenance in [the source replay](../../verification/source_replay/part_a/README.md).

[Source hashes for metadata and annotations](references/sources.json) · [Tahoe observation receipt](../../verification/source_replay/part_a/receipts/tahoe_observations.json) · [EmeraldBay and sci-Plex source manifest](../../verification/source_replay/part_a/receipts/source_verification.json)
