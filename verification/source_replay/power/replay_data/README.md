# Replay Figure S3 and check Note S17

Use Python 3.11 or later to reconstruct all 66 coordinates and their Monte Carlo intervals. This command runs offline with the standard library. The optional output directory receives the recalculated coordinates, the number of usable assignments for every point and the individual cases behind Note S17.

```sh
python3 replay_s3.py --out replay
```

The packet contains 22,000 assignment records across the two terms and eleven plotted scenarios. The input for the panel test is the permutation exceedance count b. Both terms use p = (b + 1) / 2,001 and p < 0.025, equivalent to b ≤ 49. The pooled endpoint compares each 95% interval with zero. Its paired shift is the planted estimate minus the null estimate from the same assignment. Detection of a planted compound requires BH q < 0.1.

The orange series averages the proportion of planted compounds detected *within each assignment*. The green series counts assignments with one or more detected planted compounds. The orange intervals use the standard error of those assignment-level fractions. Other intervals use Wilson's formula for event counts. A missing measurement is excluded from that endpoint's denominator. The null has no planted compounds, so a planted-hit share is undefined. Null panel rejections and pooled interval exclusions have their own valid denominators.

## Single-compound simulations

The file `single_plants.tsv.gz` contains all 4,000 single-plant cases: the 50% response increase (C1) and the 49-gene shift (G1) for Y and XIST. Of 34 planted compounds with q < 0.1, 4 had at least one permutation statistic reach the observed statistic. The other 30 had zero exceedances. The BH step-up calculation adjusts p values jointly across the eligible drugs in an assignment. A single *planted* compound can pass with positive exceedances when other compounds also pass. These paired records establish what happened at 2,000 permutations. A rate at 10,000 permutations would require running that comparison.

| Term | Scenario | Assignment | Planted compound | Exceedances | p | BH q | Other rejections |
|---|---|---:|---|---:|---:|---:|---:|
| factX | G1 | 444 | daunorubicin | 2 | 0.001499 | 0.085457 | 1 |
| factY | C1 | 78 | nutlin-3 | 4 | 0.002499 | 0.094953 | 2 |
| factY | C1 | 211 | BRD-K63750851 | 3 | 0.001999 | 0.075962 | 2 |
| factY | C1 | 723 | NVP-AUY922 | 4 | 0.002499 | 0.094120 | 2 |

`assignments.tsv.gz` contains the tested-drug count, panel exceedances, sparse-scenario planted counts and hits, and the pooled estimates and intervals used by the figure. Blank fields indicate an inapplicable or unavailable endpoint. `reference_coordinates.tsv.gz` supplies the figure values used to check the recalculation. `verification.json` records the full-precision comparisons and denominators. `provenance.json` identifies the source of this packet.
