# Independent verification of marker retention after drug exposure

The nine marker counts, library depths, cell counts and RNA classifications agree exactly across two routes through the archived screen. The main analysis summed individual cell marker records. The audit independently extracted gene-token counts from all fourteen sealed whole-transcriptome pseudobulk files, joined their well and marker-stratum totals, then recomputed classifications. The input hashes matched. There were 67,018 observed well–line entries and 182 with no captured cells, from 95,624,334 cells. An empty entry has no expression classification.

## Denominators and concentration of disagreements

The 200-cell criterion retained 58,259 treated entries from 48 lines and 1,263 vehicle entries from 47. NCI-H661 qualified only in treatment. NCI-H2122 and NCI-H596 had no eligible treated entry. All 50 lines had a pooled vehicle call consistent with the archived reference. Every test held its plate out of the vehicle reference aggregation. The marker cutoffs had previously been learned using all vehicle plates.

With every represented line weighted equally, the primary treated agreement was 98.711% over 48 lines. For the 47 lines represented in both conditions, it was 98.683% in treatment and 99.468% in vehicle. The difference was −0.785 percentage points. The line-resampling percentile interval reproduced at −2.302 to −0.012 percentage points. Removing COLO 205 reduced the difference to −0.033 percentage points. This interval describes variation between lines within this screen.

The eight-gene Y score crossed its threshold in 788 treated entries and XIST crossed in five. One entry crossed both thresholds, giving 792 joint mismatches. COLO 205 contributed 773 of 792 and all seven vehicle mismatches. Its vehicle references ranged from 3.962 to 4.001 against a Y cutoff of 4.02168. Its treated median was 4.081 (interquartile range 3.942–4.203). The original manuscript already identified COLO 205 as near its cutoff and discordant with external bulk cultures. The extension quantifies this sensitivity under drug exposure. [Original manuscript and supplement](https://github.com/MateenahJAHAN/Project-59-CMAP-Sex-Differences/tree/d9a94ddedeb65ced19c16f6cb78497d837f161b0/manuscript)

## Sensitivity to the classification rule

Two sets of cutoffs were already archived in the original fellowship analysis. The primary cutoffs came from the distribution across line–plate measurements. The alternate cutoffs came from the distribution of one pooled measurement per line. Both rows below apply to identical unfiltered cells and the same 58,259 eligible treated entries. Each uses its own archived vehicle definition and holds the target plate out of the RNA-count reference.

| Quantity | Line–plate cutoffs | Line-level cutoffs |
|---|---:|---:|
| Y cutoff | 4.02168 | 3.74152 |
| XIST cutoff | 4.68498 | 5.52651 |
| COLO 205 baseline call | Neither | Y-positive |
| COLO 205 treated entries matching baseline | 508/1,281 | 1,227/1,281 |
| Equal-line agreement across treatment | 98.711% | 99.783% |
| Treated minus vehicle, common 47 lines | −0.785 percentage points | +0.158 percentage points |
| Agreement across the two plates | 4,104/4,142 | 4,130/4,142 |
| Pairs with the same departure from baseline | 20 (neither → Y-positive) | 2 (Y-positive → neither) |

All 22 repeated departures counted under the two rules occurred in COLO 205. The directions and denominators depend on the baseline rule. The data support a boundary-sensitive RNA classification. Reproducible changes in chromosome state would require an independent culture, direct chromosome measurement and an established decision rule. [Full alternative-threshold results](independent_check_thresholds.json) · [Threshold comparison code](independent_check_thresholds.py)

## Sensitivity to marker-dependent cell exclusion

The primary analysis kept every cell passing atlas QC. As a retrospective sensitivity, we applied the earlier state-dependent contamination rule. In baseline-Y-low lines it removed cells carrying three or more Y-marker UMIs. In baseline-XIST-low lines it removed cells carrying three or more XIST UMIs. Across the screen this removed 1,384,481 cells. The reference and test wells were recalculated. On the same 58,201 treated entries eligible under both processing choices, 790 disagreed with the vehicle call before filtering and 167 after. Because the rule removes observations that could overturn the baseline, the uncensored comparison is the appropriate estimate of observed call retention.

Same-plate vehicle references gave 99.007% equal-line agreement across 58,188 eligible treated entries. Library thinning gave 98.984% with the specified seed. Two further seeds gave 99.033% and 98.992%. Shared-well RNA, composition and biological differences remain possible explanations for the residual disagreement. [Control results](independent_check_controls.json)

## Replicate plates

Plates 6 and 14 shared 94 labeled drug–dose conditions. Ninety-three had at least one line with 200 cells in both wells. The 5 µM L-thyroxine well on plate 14 had 430 cells in all, with a maximum of 37 in any one line. Exact matching by name, concentration and line gave 4,142 eligible pairs across 47 lines. Under the primary rule all 38 cross-plate disagreements involved COLO 205. In the 20 repeatedly discordant pairs described above, at least one Y-score crossing lay only 0.00017 log-score units above the threshold on plate 6. Plates 6 and 14 shared well positions. These pairs measure cross-plate reproducibility within the atlas. [Tahoe plate design audit](https://arxiv.org/abs/2611.01604)

## Experimental precedents

Combined decitabine and phenylbutyrate treatment induced XIST expression approximately sixfold in Caco-2 cells. A demethylation experiment in mouse–human somatic hybrid cells also induced XIST from an active human X. These studies show why stable binary calls within one screen cannot establish invariant transcription, chromosome dosage or inactivation. [Caco-2 treatment study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10199547/) · [Demethylation study](https://pubmed.ncbi.nlm.nih.gov/9611238/)

[Independent recount](independent_check.py) · [Numerical audit](independent_check.json) · [Control and bootstrap verification](independent_check_controls.py)
