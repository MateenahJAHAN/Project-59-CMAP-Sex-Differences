# Analysis code and derived data

## Marker persistence under treatment

The [marker stability follow-up](marker_stability/) tests the vehicle classifications in 58,259 treated line–wells. It includes an independent COLO 205 RNA and DNA case study, two-plate landmark-gene responses, the donor-label comparison and a relative-recovery analysis.

## Tahoe-100M expression calls

[tahoe/](tahoe/) contains the local implementation used to call Y-gene and XIST expression and evaluate the marker-negative cell analyses. It includes the scoring, filtering, sensitivity and plotting code. The local analyses used Python 3.12.11 with NumPy 2.5.4, pandas 3.1.1 and SciPy 1.18.1. [Dependencies](tahoe/requirements.txt)

Begin with the [Tahoe-100M dataset](https://huggingface.co/datasets/tahoebio/Tahoe-100M/tree/349e2a55f5538a3e73b32a06c6bdd8e9e87011e8). The caller expects vehicle counts and metadata after reduction to `DMSO_TF` cells with `full` status, including per-cell UMI totals, marker counts and line-by-well pseudobulks. It also expects external line annotations for the agreement checks. Obtain the sources and reference files listed in the paper's supplementary methods. The input locations in [`common.py`](tahoe/common.py) and [`01_gene_chrom.py`](tahoe/01_gene_chrom.py) are read from environment variables. [`run_all.sh`](tahoe/run_all.sh) shows the execution order. The original shell launcher expects the interpreter at `$P59_DATA/venv/bin/python`. Set `P59_PARTA_OUT` to a fresh output directory.

| Variable | Folder it names | Files the scripts read there |
| --- | --- | --- |
| `P59_DATA` | Data root | `tahoe/dmso/cells_markers.parquet`, `pseudobulk_counts.parquet`, `pseudobulk_groups.tsv`, `genes.tsv`, `lines.tsv` and `six_lines/plate*.parquet`; `tahoe/index/meta/sample_metadata.parquet`; `venv/bin/python` for the launcher |
| `P59_ANALYSIS` | Analysis checkout | `inputs/C_cellline_sex_labels.tsv`, `inputs/A_cellline_sex_states.tsv` and `inputs/tahoe_pilot_states.tsv`. These three line tables are not distributed in this repository. |
| `P59_PARTA_OUT` | Output folder | `downloads/` and `provenance_downloads.tsv`, written by [`00_downloads.sh`](tahoe/00_downloads.sh) |

```sh
export P59_DATA=/path/to/data P59_ANALYSIS=/path/to/analysis-checkout P59_PARTA_OUT=/path/to/new-output
zsh analysis/tahoe/run_all.sh
```

The two synthetic tests, `test_calling_synthetic.py` and `test_mosaic_synthetic.py`, need no data and no variables.

[Download expression states for Tahoe-100M and EmeraldBay](data/tahoe-emeraldbay-expression-states.tsv). This file contains 102 atlas rows selected from Supplementary Table S1. [Licence and source credits](data/README.md)

## LINCS2020 drug-response analysis

[lincs/](lincs/) contains the model and permutation engine, analysis entry point and numerical configuration. The original run used Python 3.13.12 on macOS arm64 with NumPy 2.5.3, pandas 3.0.5 and SciPy 1.18.1. It used 10,000 permutations and passed `--workers 12` to the runner. [Dependencies](lincs/requirements.txt) · [Configuration](lincs/config/G_defaults.json) · [Entry point](lincs/src/G/run_g.py)

The runner expects the original LINCS2020 level-5 signature tensors and level-3 tensors referenced to same-plate DMSO wells, plus the prepared line table, line-by-compound masks, calibration data and associated covariates. The original reduction selected 978 L1000 landmark genes and organized signatures by line, compound and panel. P1 uses ASG signatures at the `1.11 uM` source label, predominantly measured at 1.0 µM and 24 hours. The well-level reduction subtracts same-plate DMSO means before forming the level-3 tensors. The configuration specifies the input layout. The original `real` mode also verifies the input and source checksums and the plan held in the internal Git tag `partB-prereg-v1`. The historical tag name is a private version-control label. Re-executing that mode requires the original analysis checkout, plan and reduced inputs under the providers' terms. This snapshot exposes the implemented methods and settings. See the [manuscript and supplement](../manuscript/) for the figures, statistical tables and input sources.

The [repository software licence](../LICENSE) covers the code. Derived data have a [separate licence and attribution](data/README.md).
