#!/bin/sh
# Legacy entry point for an already prepared Tahoe DMSO count extraction and A/C inputs.
set -eu
: "${P59_DATA:?Set P59_DATA to the run's data directory}"
: "${P59_ANALYSIS:?Set P59_ANALYSIS to the directory containing inputs/}"
code=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
export P59_PARTA_OUT=${P59_PARTA_OUT:-"$P59_DATA/tahoe/partA_local"}
py=${P59_PYTHON:-"$P59_DATA/venv/bin/python"}
if [ ! -x "$py" ]; then py=python3; fi
mkdir -p "$P59_PARTA_OUT/logs" "$P59_PARTA_OUT/figures"
sh "$code/00_downloads.sh"
for stage in 01_gene_chrom test_calling_synthetic 02_call_states 03_sensitivity 04_agreement_audit \
             05_ambient_foreign test_mosaic_synthetic 06_mosaic 07_identity 08_mosaic_rule 09_figures 10_gene_level
do
    printf 'Running %s\n' "$stage"
    "$py" -I "$code/$stage.py"
done
