CMap LINCS2020 inputs

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
