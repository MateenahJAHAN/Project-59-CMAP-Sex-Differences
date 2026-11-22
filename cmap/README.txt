CMAP DATA FOLDER
================

  This folder holds the CMAP (Connectivity Map / LINCS 2020) metadata files
used in the sex-differences drug-response pipeline.

    SMALL FILES (OK to store on GitHub):
    geneinfo_beta.txt        (~1.1 MB)  Entrez <-> Ensembl <-> Symbol map
      cellinfo_beta.txt        (~38 KB)   Cell line donor sex + lineage
      compoundinfo_beta.txt    (optional) Drug mechanism reference

LARGE FILES - DO NOT PUSH TO GITHUB (over 100 MB limit):
    siginfo_beta.txt                          (~444 MB) experiment metadata
  level5_beta_trt_cp_n720216x12328.gctx     (~33 GB)  drug z-score matrix
  -> keep these in Colab / Google Drive only, not in this repo.

SOURCE: CLUE / LINCS 2020 - https://clue.io
  s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/

  
This directory contains three archived metadata files.
  cellinfo_beta (2).txt        Cell-line donor-sex annotation and lineage
  geneinfo_beta (1).txt        CMap feature IDs, Ensembl IDs and symbols
  compoundinfo_beta (1).txt    Compound identifiers and annotations

The original Colab code expects these as cellinfo_beta.txt, geneinfo_beta.txt
and compoundinfo_beta.txt. Use copies with those names in a local data directory.

The full signature metadata and level-5 expression matrix are external inputs.
https://clue.io/data/CMap2020#LINCS2020
https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/siginfo_beta.txt
https://s3.amazonaws.com/macchiato.clue.io/builds/LINCS2020/level5/level5_beta_trt_cp_n720216x12328.gctx
