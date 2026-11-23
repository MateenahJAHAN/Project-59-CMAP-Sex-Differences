# Part A, check A3: depth-matched shares of marker-negative cells across cultures

Written on 2026-11-15 (US Eastern time). Committed before any depth-matched share was computed. The git history is the proof of order.

## The frozen text this file carries out

`PartA_frozen_rules.md`, section B, row A3: "A549 Y-negative and SW48 XIST-negative shares in EmeraldBay and sci-Plex within 0.15 of Tahoe's, on depth-matched data (binomial thinning of UMIs to Tahoe's depth distribution), under the same estimators".

## What had been seen before this file was written

- Shares without depth matching, from task F1's Poisson mixture: A549 0.385 (EmeraldBay), 0.364 (sciPlex2), 0.452 (sciPlex3) and 0.505 (sciPlex4); SW48 0.417 (EmeraldBay). Source: `03_cloud_results/task-F1/results/F1/F1_agreement.tsv`, the share rows.
- Frozen Tahoe ranges across estimators (a) to (c) and both callers: A549 0.193 to 0.426, SW48 0.259 to 0.569. Source: `results/PartA/mosaic_resolution.tsv`, columns `share_lo_both` and `share_hi_both`.
- Mean UMIs per A549 vehicle cell: Tahoe 1,906, EmeraldBay 3,311, sciPlex2 4,953, sciPlex3 1,598, sciPlex4 1,460. Source: the packet's `numbers.json`, keys `depth_*_a549`, each pointing to its source file.
- No depth-matched share had been computed.

## Deviation from the frozen text, written before running

Thinning can only remove UMIs. sciPlex3 and sciPlex4 are shallower than Tahoe, so their cells cannot be thinned "to Tahoe's depth distribution". In every pair, the deeper dataset is therefore thinned to the shallower dataset's depth distribution.
- EmeraldBay and sciPlex2 are deeper than Tahoe, so for them this is the frozen text.
- For sciPlex3 and sciPlex4, Tahoe is the side that is thinned.

## Data

- **Tahoe.** The local caller's vehicle cells: `WORK/tahoe/dmso/cells_markers.parquet` (DMSO_TF, full pass). These are the same cells as the frozen local table, after its foreign-cell filter (`n_vehicle_cells` in `partA_local_states.tsv`).
  - Per-plate, per-UMI background from `partA_local_plate_background.tsv`. Y uses the line set "Y-negative lines with both arms low (partial-Y lines removed)"; XIST uses "XIST-negative lines".
  - The cloud caller's per-cell data are not on this machine, so A3 uses one caller.
- **EmeraldBay.** `F1_emerald_vehicle_cells.tsv.gz`, rows with drug DMSO_TF only, with states from `F1_emerald_states.tsv`.
- **sci-Plex.** A per-cell vehicle table rebuilt with F1's labelling rules and checked against F1's unit totals (`PROJECT59_DEEP_DIVE_2026-11-13/08_partA_A2_A3/sciplex_cells/`), with states from `F1_sciplex_states.tsv`.
- **Background in EmeraldBay and sci-Plex,** computed like Tahoe's. Per unit (EmeraldBay sample, sci-Plex batch), it is marker UMIs over total UMIs, in vehicle cells with fewer than 3 marker UMIs, from lines that are negative for that marker in that dataset.
  - For Y, these are lines with arm state "Y lost".
  - For XIST, these are lines flagged XIST-negative.
- **Foreign-cell filter, as in the frozen callers.** It is applied before thinning.
  - When Y is tested in a line that is XIST-negative in that dataset, drop cells with 3 or more XIST UMIs.
  - When XIST is tested in a Y-negative line, drop cells with 3 or more Y UMIs.
  - The tested marker is never used to filter.

## Pairs

| Pair | Line, marker | Tahoe against | Role |
|---|---|---|---|
| 0 | A549, Y | EmeraldBay | test |
| 1 | A549, Y | sciPlex2 | test |
| 2 | A549, Y | sciPlex3 | test |
| 3 | A549, Y | sciPlex4 | test |
| 4 | SW48, XIST | EmeraldBay | test; EmeraldBay has fewer than the 1,000 cells the mosaic tests need, and the report says so |
| 5 | NCI-H460, Y | EmeraldBay | context only, no decision |
| 6 | LoVo, Y | EmeraldBay | context only, no decision |

Pairs 5 and 6 are lines without a known Y-negative subpopulation. They show how far thinning alone moves each estimator.

## Depth matching

1. In each pair, the dataset with the higher mean UMIs per cell is thinned.
2. Rank its cells by depth. The cell at rank i of n gets a target depth t: the shallower dataset's empirical depth quantile at (i - 0.5)/n, from numpy.quantile with method "linear".
3. If t is below the cell's depth d, keep every UMI independently with probability t/d. Marker UMIs become y' ~ Binomial(y, t/d) and the other UMIs become Binomial(d - y, t/d); the new depth d' is their sum. If t is not below d, keep the cell as it is.
4. Run ten thinning replicates, with seed 20261115 + 100 x pair + replicate.
5. The shallower dataset is not changed. The per-UMI background is not changed by thinning.

## Estimators

These are the local caller's functions in `partA_local/code/mosaic_models.py`, unchanged:
- (a) the zero-marker fraction among cells with at least 3,000 UMIs, and among cells with at least 8,000 UMIs;
- (b) the negative-binomial-plus-background mixture with a logistic share, taken at 2,000 UMIs;
- (c) the Poisson mixture with background.

When each is used:
- (a) is used at a threshold only if both sides have at least 100 cells above it after matching.
- (b) is used only if its fit converged and is not at a bound on either side.

## Decision rule

1. For each pair and each estimator used on both sides, take the absolute difference between the two shares. On the thinned side, use the median over the 10 replicates. The estimator agrees if the difference is at most 0.15.
2. A pair agrees if every estimator used agrees.
3. sciPlex2, sciPlex3 and sciPlex4 come from one study and may share a stock. They are reported one by one. Together they count as one culture, which agrees if at least 2 of the 3 experiments agree.
4. A549 passes A3 if it agrees in EmeraldBay and in sci-Plex. SW48 passes if it agrees in EmeraldBay.
5. Report k of 3 cultures (Tahoe, EmeraldBay, sci-Plex) for A549.

How the result reads:
- **A3 passes:** "Once depth is matched, A549's Y-negative cells appear at similar shares in [k] of 3 cultures."
- **A3 fails:** "The A549 Y-negative share depends on the stock or the dataset, not only on the line name."

In both cases the mosaic claim stays dropped (`results/PartA/README.md`); A3 does not reopen it.

## Also reported, descriptive

- The same estimators without depth matching, for every dataset.
- The range over the thinning replicates.
- Cells and deep cells on each side.
- Fit convergence for (b) and (c).

## Addendum 1: background batches in sci-Plex

Written before the sci-Plex per-cell table existed, so before any sci-Plex share under this rule could be computed. "sci-Plex batch" above means F1's background batch, the `bg_batch` column of `F1_ambient.tsv`:
- sciPlex3: the RT plate (rt1 to rt4), because its culture plates hold one line each.
- sciPlex4: the whole sample.
- sciPlex2: it has only A549, so it has no Y-negative line. Its Y background is the pooled sciPlex3 rate under the definition above, as in F1's plan (`F1_REPORT.md`, Decisions 10).
