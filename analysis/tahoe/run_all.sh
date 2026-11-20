#!/bin/zsh
# Part A local caller (state-caller-local): reproduce every output from the DMSO extraction.
# Usage: zsh run_all.sh            (writes to $P59_DATA/tahoe/partA_local)
#        P59_PARTA_OUT=/some/new/folder zsh run_all.sh   (fresh folder; 02 refuses to overwrite a state table)
set -euo pipefail
CODE=${0:A:h}
PY=($P59_DATA/venv/bin/python -I)
export P59_PARTA_OUT=${P59_PARTA_OUT:-$P59_DATA/tahoe/partA_local}
mkdir -p $P59_PARTA_OUT/logs $P59_PARTA_OUT/figures
step() { local t0=$(date +%s); echo "== $1"; shift; "$@"; echo "   $(( $(date +%s) - t0 )) s" | tee -a $P59_PARTA_OUT/logs/run_all_runtime.txt; }
D=$P59_PARTA_OUT/downloads
if [[ -f $D/pie/descriptions.json && -f $D/cellosaurus/CVCL_C466.txt && -f $D/ensembl/Homo_sapiens.GRCh38.110.gtf.gz ]]; then
  echo "== downloads present; not fetched again"
else
  step 00_downloads zsh $CODE/00_downloads.sh
fi
step 01_gene_chrom          $PY $CODE/01_gene_chrom.py
step test_calling_synthetic $PY $CODE/test_calling_synthetic.py
step 02_call_states         $PY $CODE/02_call_states.py
step 03_sensitivity         $PY $CODE/03_sensitivity.py
step 04_agreement_audit     $PY $CODE/04_agreement_audit.py
step 05_ambient_foreign     $PY $CODE/05_ambient_foreign.py
step test_mosaic_synthetic  $PY $CODE/test_mosaic_synthetic.py
step 06_mosaic              $PY $CODE/06_mosaic.py
step 07_identity            $PY $CODE/07_identity.py
step 08_mosaic_rule         $PY $CODE/08_mosaic_rule.py
step 09_figures             $PY $CODE/09_figures.py
step 10_gene_level         $PY $CODE/10_gene_level.py
