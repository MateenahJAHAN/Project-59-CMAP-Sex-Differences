# Manuscript figures rebuilt from the study data

The [source comparison](source_figure_comparison.json) checks every plotted field and compares each rebuilt PDF with the released file, extracted text and pixels at 144 dpi. The separate [baseline check](baseline_comparison.tsv) rendered all seven figures from archived plotting inputs.

| Figure | Evidence from the rerun |
| --- | --- |
| 1. Study design | The schematic has six labelled steps. The rebuilt PDF matches the release. |
| 2. Donor labels and expression | All 124 plotted rows and the PDF match. The rerun references its own source tables and its newly generated marker list. |
| 3. Agreement across sources | All 232 plotted rows reproduced, including every line and every contingency-table count. The six sci-Plex line–experiment rows are listed in a different order. The figure PDF is byte-identical. |
| S1. Background RNA | All 180 plotted rows and the PDF reproduced from the plate and cell counts. |
| S2. A549 | All 36 plotted rows and the PDF reproduced from the depth curve and cross-culture comparisons. |

The [archived](archived_plot_data/) and [rebuilt plotting tables](rebuilt_plot_data/) preserve the labels and numbers. Source files are identified by filename. To check newly generated figures with Python and PyMuPDF, pass the directory containing the released PDFs, the archived plotting tables, and the directory with your rebuilt PDFs and plotting tables to the [comparison script](verify_source_figures.py).
