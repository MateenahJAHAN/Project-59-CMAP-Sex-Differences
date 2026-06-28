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

  
