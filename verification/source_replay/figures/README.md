# Manuscript figures rebuilt from the study data

The [source comparison](source_figure_comparison.json) checks every plotted field and compares each rebuilt PDF with the released file, extracted text and pixels at 144 dpi. The separate [baseline check](baseline_comparison.tsv) rendered all seven figures from archived plotting inputs.

| Figure | Evidence from the rerun |
| --- | --- |
| 1. Study design | The schematic has six labelled steps. The rebuilt PDF matches the release. |
| 2. Donor labels and expression | All 124 plotted rows and the PDF match. The rerun references its own source tables and its newly generated marker list. |
| 3. Agreement across sources | All 232 plotted rows reproduced, including every line and every contingency-table count. The six sci-Plex line–experiment rows are listed in a different order. The figure PDF is byte-identical. |
| 4. Drug response | The rebuilt figure, all 274 plotted rows, 30,000 permutation draws, and all 12 estimates and intervals matched. [Null distributions](../part_b/figure4_nulls_check.json) · [Estimates](../part_b/figure4_intervals.tsv) |
| S1. Background RNA | All 180 plotted rows and the PDF reproduced from the plate and cell counts. |
| S2. A549 | All 36 plotted rows and the PDF reproduced from the depth curve and cross-culture comparisons. |
| S3. Sensitivity to planted effects | The 1,000-assignment replay reproduced all 66 plotted coordinates and interval endpoints. The bytes of the figure PDF matched. Its revised caption distinguishes the fraction of assignments with any detection from the mean fraction of planted compounds detected. [Assignment results](../power/) · [Revised caption](../../../manuscript/supplement.pdf) |

The [archived](archived_plot_data/) and [rebuilt plotting tables](rebuilt_plot_data/) preserve the labels and numbers. Source files are identified by filename. To check newly generated figures with Python and PyMuPDF, pass the directory containing the released PDFs, the archived plotting tables, and the directory with your rebuilt PDFs and plotting tables to the [comparison script](verify_source_figures.py).

The [plotting programs](../../../manuscript/plotting/) build the figures from released numerical inputs. Figure S1 uses two [plate and cell-line summaries](source_data/) calculated by [the Tahoe background-RNA analysis](../../../analysis/tahoe/05_ambient_foreign.py) from the CC0 Tahoe-100M release. The [source-run build](../../../manuscript/from_run.py) reads the corresponding tables from a completed external working directory.
