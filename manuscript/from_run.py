"""Collect the regenerated measurements, model fits and simulation results for typesetting."""
import argparse
from pathlib import Path
import shutil
import sys

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "plotting"))
from agreement import from_sources
sys.path.insert(0, str(HERE.parent))
from analysis.tahoe.common import Y8, YP, YQ


def copy(src, dst):
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)


def write(frame, out):
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, sep="\t", index=False)


def prepare(work):
    base = work / "results"
    tahoe, bulk, external = base / "tahoe/latest", base / "C", base / "F1"
    out = work / "plot-inputs"
    a = out / "part_a/tables"
    for suffix in ("states", "cutoffs", "label_audit_summary", "depth_curve"):
        file = "partA_local_" + suffix + ".tsv"
        copy(tahoe / file, a / file)
    for suffix in ("ambient", "foreign_cells"):
        file = "partA_local_" + suffix + ".tsv"
        copy(tahoe / file, out / "figures/source_data" / file)
    copy(base / "depth_matching/A3_pairs.tsv", a / "A3_pairs.tsv")
    copy(base / "depth_matching/A3_summary.json", out / "part_a/receipts/A3_summary.json")
    measured = pd.read_csv(tahoe / "partA_local_gene_level_scores.tsv", sep="\t", nrows=0)
    if not set(Y8).issubset(measured.columns):
        raise ValueError("The measured gene table is missing a Y-score marker")
    markers = [dict(element="gene_list", name=name, text=", ".join(genes))
               for name, genes in (("y", Y8), ("yp", YP), ("yq", YQ))]
    plotted = out / "figures/rebuilt_plot_data"
    write(pd.DataFrame(markers), plotted / "fig2_label_vs_pattern_plotted.tsv")
    str_dir = base / "str"
    counts = str_dir / "figure3_agreement_by_source.tsv"
    comparisons = from_sources({"tahoe": tahoe / "partA_local_states.tsv", "tahoe_cut": tahoe / "partA_local_cutoffs.tsv",
                                "c_labels": bulk / "C_cellline_sex_labels.tsv", "c_cut": bulk / "C_expression_cutoffs.tsv",
                                "emerald": external / "F1_emerald_states.tsv", "sciplex": external / "F1_sciplex_states.tsv",
                                "f1_cut": external / "F1_cutoffs.tsv", "str_by_source": counts,
                                "str": str_dir / "figure3_supp_table_amelogenin_by_line.tsv"})
    copy(counts, a / counts.name)
    write(comparisons["lines"].assign(element="line"), plotted / "fig3_agreement_plotted.tsv")

    b = out / "part_b"
    for file in ("figure4_intervals.tsv", "figure4_nulls.tsv.gz"):
        copy(base / "figure4" / file, b / file)
    copy(base / "figure4/plot_inputs/part_b/fig4_partB_plotted.tsv", b / "fig4_partB_plotted.tsv")
    copy(work / "power/results/figS3_power_values.tsv", out / "power/full_run/figS3_power_values.tsv")
    print(f"Prepared figure inputs from {work}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    work = args.work.resolve()
    if work.is_relative_to(HERE.parent):
        parser.error("Use the external analysis work directory")
    prepare(work)
