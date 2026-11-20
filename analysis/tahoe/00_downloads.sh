#!/bin/zsh
# Part A local caller: external downloads (PIE context file, Cellosaurus entry for hTERT-HPNE,
# Ensembl gene annotation for chromosome of each gene). Each goes into its own directory.
# Writes a provenance TSV with url, local path, bytes, sha256 and UTC access time.
set -euo pipefail
OUT=${P59_PARTA_OUT:-$P59_DATA/tahoe/partA_local}
DL=$OUT/downloads
PROV=$OUT/provenance_downloads.tsv
mkdir -p $DL/pie $DL/cellosaurus $DL/ensembl
if [[ ! -f $PROV ]]; then
  printf "item\turl\tlocal_path\tbytes\tsha256\taccessed_utc\tnotes\n" > $PROV
fi
fetch() {
  local item=$1 url=$2 dest=$3 notes=$4
  local t=$(date -u +%Y-%m-%dT%H:%M:%SZ)
  curl -sSfL --retry 3 --retry-delay 5 -m 600 -o "$dest" "$url"
  local b=$(stat -f %z "$dest")
  local h=$(shasum -a 256 "$dest" | cut -d' ' -f1)
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$item" "$url" "$dest" "$b" "$h" "$t" "$notes" >> $PROV
  echo "$item $b $h"
}
fetch pie_descriptions "https://huggingface.co/datasets/arcinstitute/PIE_sources/resolve/main/context_text/descriptions.json" \
  $DL/pie/descriptions.json "rules expect 62,725 B, sha256 prefix 2b8875977f559faf, 88 contexts"
fetch cellosaurus_CVCL_C466 "https://api.cellosaurus.org/cell-line/CVCL_C466?format=txt" \
  $DL/cellosaurus/CVCL_C466.txt "Cellosaurus entry of hTERT-HPNE (no Task C row); SX line used as the annotation"
fetch ensembl_gtf_110 "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz" \
  $DL/ensembl/Homo_sapiens.GRCh38.110.gtf.gz "chromosome of each gene for test (f); not listed in the rules (deviation)"
fetch ensembl_gtf_110_checksums "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/CHECKSUMS" \
  $DL/ensembl/CHECKSUMS "Ensembl BSD sum checksums for the GTF"
