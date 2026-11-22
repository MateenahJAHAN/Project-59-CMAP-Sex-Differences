GTEx DATA FOLDER
================

This folder holds the GTEx (Genotype-Tissue Expression) v8 sex-biased gene data.

MAIN FILE:
  signif.sbgenes.txt   <- MAIN file (sex-biased genes, all 44 tissues)

OTHER FILES (extracted from GTEx_Analysis_v8_sbgenes.tar.gz):
  sbgenes.txt
  Lung.sbgenes.txt
  Breast_Mammary_Tissue.sbgenes.txt
  (other per-tissue .sbgenes.txt files)
  GTEx_Analysis_v8_sbgenes.tar.gz   (original archive)

SOURCE: GTEx Portal - https://gtexportal.org/home/datasets

NOTE: Add the actual data files into this folder.

Download GTEx_Analysis_v8_sbgenes.tar.gz and extract signif.sbgenes.txt into
your analysis data directory. The GTEx files are external inputs.

https://gtexportal.org/home/downloads/adult-gtex
https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/sb-eqtl/GTEx_Analysis_v8_sbgenes.tar.gz

The Phase 4 and Phase 6 split libraries use positive effsize for genes with
higher expression in females and negative effsize for males.
https://doi.org/10.1126/science.aba3066
