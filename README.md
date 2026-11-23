# Sex-chromosome expression and annotated donor sex in drug-perturbation cell lines

Mateenah Jahan and Marouen Ben Guebila · Fatima Fellowship, Project 59

<a id="manuscript-and-figures"></a>
<a id="current-status"></a>

**[Manuscript (PDF)](manuscript/manuscript.pdf) · [Supplement (PDF)](manuscript/supplement.pdf) · [Figures and source data](manuscript/README.md) · [Reproduce the study](REPRODUCING.md)**

Can donor sex tell us which sex-chromosome genes a cultured cancer cell line still expresses? And would measuring those genes help explain the response to a drug? The paper tests these questions in two parts.

| | Data and comparison | Finding |
|---|---|---|
| **A. Measure the markers** | We measured Y-linked transcripts and XIST in Tahoe-100M vehicle cells, then compared each line with its donor annotation and independent RNA measurements. | In 27 of 49 annotated lines, the expressed markers differed from the combination expected from donor sex. Independent bulk RNA agreed with the Tahoe classification in 47 of 49 lines. [Expression calls](analysis/data/tahoe-emeraldbay-expression-states.tsv) · [Analysis](analysis/part_a/) |
| **B. Test drug response** | Within cell lineages, we tested whether Y or XIST expression predicted the transcriptional response to 114 compounds across 41 LINCS2020 lines. | Neither marker resolved a panel-wide difference (Y *p* = 0.7231, XIST *p* = 0.1127). The intervals allowed substantial effects, and XIST also showed an association in the vehicle control (*p* = 0.0323). [Drug and control results](verification/source_replay/part_b/) · [Analysis](analysis/lincs/) |

[![Study design, from measuring expression in vehicle cells to testing drug response within cell lineages.](manuscript/figures/fig1_design.png)](manuscript/figures/fig1_design.pdf)

<a id="how-the-analysis-developed"></a>
<a id="what-the-exploratory-screen-found"></a>
<a id="the-eight-fellowship-phases"></a>

## How the fellowship reached this question

The first experiment compared vorinostat profiles. Phases 2–6 expanded the screen to 50 compounds and changed how tissue enrichment was tested. In Phase 7, 17 female-donor lines and 18 male-donor lines contributed to the own-tissue category, yet 2,426 of its 2,656 signatures came from the female-donor lines. Repeated measurements of a few lines could dominate the comparison. The manuscript therefore returns to the cell line as the unit of evidence and measures the markers directly.

The [fellowship phase guide](fellowship/README.md) follows each experiment from question to result, including which calls changed on replay. Phase 8 examined adverse-event reporting for 18 drugs in FAERS. That analysis has a separate outcome and denominator.

<a id="marker-persistence-in-drug-treated-cultures"></a>

## Follow-up in treated Tahoe-100M cultures

In the separate marker-persistence study, the vehicle classification persisted in 98.71% of eligible treated line–wells when each cell line carried equal weight. COLO 205 accounted for 773 of 792 disagreements across 58,259 line–wells. Independent RNA and DNA assays showed stronger signal on the measured long arm of chromosome Y than on the short arm. [Follow-up paper and supplement](manuscript/followups/marker_stability/marker_stability.pdf) · [Code, data and authorship](analysis/marker_stability/) · [COLO 205 assays](analysis/marker_stability/independent_colo205/)

<a id="data-and-setup"></a>
<a id="reproducing-the-analyses"></a>

## Run the analyses

From a fresh clone, `python3 reproduce.py verify` checks the reported results using the files in this repository. For a source-data run, use `python3 reproduce.py all --work ../p59-run`. The [reproduction guide](REPRODUCING.md) lists installations, data sources, storage, individual stages and the limits of the historical reconstruction. [Code licence](LICENSE) · [Data licences and attribution](analysis/data/README.md)
