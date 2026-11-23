GTEx v8 tissue gene sets

The fellowship screens use the GTEx v8 sex-biased expression archive. The
fellowship runner downloads and checks the archive before extracting the
44 tissue gene sets. In the split libraries, positive effsize denotes higher
expression in females and negative effsize denotes higher expression in males.

From the repository root:
python3 reproduce.py fellowship --work ../p59-run

Source archive
https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/sb-eqtl/GTEx_Analysis_v8_sbgenes.tar.gz

Oliva et al., Science (2020)
https://doi.org/10.1126/science.aba3066
