"""Render the published plots and typeset the main paper and supplementary methods."""
import argparse
from importlib import import_module
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
FIGURES = ("fig1_design", "fig2_label_vs_pattern", "fig3_agreement", "fig4_partB",
           "figS1_ambient", "figS2_a549", "figS3_power")


def draw(destination, inputs):
    sys.path.insert(0, str(HERE / "plotting"))
    import p59style

    for name in FIGURES:
        module = import_module(name)
        if name == "fig1_design":
            figure, _ = module.build()
            module.save(figure, destination / name)
        else:
            figure, _ = module.build(inputs)
            p59style.save(figure, str(destination / name))
        print(f"Rendered {name}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="Directory for generated PDFs and plots")
    parser.add_argument("--inputs", type=Path, default=HERE.parent / "verification/source_replay",
                        help="Root of the published numerical tables, or an equivalent reconstructed directory")
    choice = parser.add_mutually_exclusive_group()
    choice.add_argument("--figures-only", action="store_true", help="Render seven figures and skip typesetting")
    choice.add_argument("--use-released-figures", action="store_true", help="Typeset the papers with the supplied figure PDFs")
    args = parser.parse_args()
    target = args.out.resolve()
    if target == HERE or HERE in target.parents:
        parser.error("Choose a build directory outside the manuscript source.")
    tectonic = shutil.which("tectonic")
    if not args.figures_only and not tectonic:
        parser.error("Typesetting requires the Tectonic executable on PATH. See manuscript/README.md.")
    figures = target / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    if args.use_released_figures:
        for name in FIGURES:
            shutil.copy2(HERE / "figures" / f"{name}.pdf", figures)
    else:
        draw(figures, args.inputs.resolve())
    if not args.figures_only:
        for name in ("manuscript", "supplement"):
            source = target / f"{name}.tex"
            shutil.copy2(HERE / source.name, source)
            subprocess.run([tectonic, "--keep-logs", source.name], cwd=target, check=True)
            print(target / f"{name}.pdf", flush=True)


if __name__ == "__main__":
    main()
