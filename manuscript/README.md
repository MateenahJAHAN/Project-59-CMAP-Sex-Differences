# Paper and figures

**Sex-chromosome expression and annotated donor sex in drug-perturbation cell lines**

Mateenah Jahan and Marouen Ben Guebila

[Paper](manuscript.pdf) · [Supplement](supplement.pdf) · [Manuscript source](manuscript.tex) · [Supplement source](supplement.tex)

Supplementary Table S1 contains 1,694 line–source records for 1,123 distinct cell lines. The [public atlas subset](../analysis/data/) contains 102 rows. The Part A command in the [reproduction guide](../REPRODUCING.md) writes the complete table to `results/supplement/supplementary_table_S1.tsv` under the work directory, using the original providers’ data.

| Figure | What it shows | Numerical source |
|---|---|---|
| [1. Study design](figures/fig1_design.pdf) | Expression-state measurement and the within-lineage drug comparison | [Schematic source](plotting/fig1_design.py) |
| [2. Donor annotation and expression](figures/fig2_label_vs_pattern.pdf) | Marker expression by cell line and donor label | [Tahoe tables](../verification/source_replay/part_a/tables/) |
| [3. Agreement across datasets](figures/fig3_agreement.pdf) | Independent RNA measurements and donor labels for the same lines | [Per-line comparisons](../verification/source_replay/figures/rebuilt_plot_data/fig3_agreement_plotted.tsv) |
| [4. Drug response](figures/fig4_partB.pdf) | Within-lineage permutations and pooled response-size intervals | [30,000 permutation draws](../verification/source_replay/part_b/figure4_nulls.tsv.gz) · [Intervals](../verification/source_replay/part_b/figure4_intervals.tsv) |
| [S1. Background RNA](figures/figS1_ambient.pdf) | Plate background and marker detection in individual cells | [Tahoe summary tables](../verification/source_replay/figures/source_data/) |
| [S2. A549](figures/figS2_a549.pdf) | Depth dependence and agreement of marker-negative fractions across cultures | [Depth curve](../verification/source_replay/part_a/tables/partA_local_depth_curve.tsv) · [Paired estimates](../verification/source_replay/part_a/tables/A3_pairs.tsv) |
| [S3. Sensitivity to planted effects](figures/figS3_power.pdf) | Detection rates and intervals under within-lineage relabelling | [Assignment records and recount](../verification/source_replay/power/) |

From the repository root, run `python3 reproduce.py paper --work ../p59-run`. This draws all seven figures and typesets both documents under `../p59-run/paper/`. It requires uv and Tectonic. Tectonic downloads its TeX bundle on first use. The release used Tectonic 0.15.0 and Arial for figures. Helvetica and DejaVu Sans are fallbacks on systems without Arial. To preserve the original fonts on such a system, run `python3 manuscript/build.py --use-released-figures --out ../p59-typeset` with Tectonic on your PATH.

With PyMuPDF installed, `python3 manuscript/compare_build.py --build ../p59-run/paper-from-run` compares the text, link targets and rendered pages of both PDFs and all seven figures against the release. The report identifies each changed page and records PDF byte equality separately. All 18 manuscript pages, 32 supplement pages and 81 link annotations matched at 144 dpi; all seven figure PDFs also matched byte-for-byte. [Full comparison](../verification/source_replay/manuscript_review/portable_build.json)

The [fellowship phases](../fellowship/README.md) document the analyses that preceded this paper. The [treated-culture follow-up](followups/marker_stability/marker_stability.pdf) has its own study design and authorship.
