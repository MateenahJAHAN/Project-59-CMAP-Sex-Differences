Broad LINCS2020 metadata used in the fellowship screens

cellinfo_beta (2).txt       Cell-line identifiers, donor annotation and lineage
geneinfo_beta (1).txt       CMap feature IDs, Ensembl IDs and gene symbols
compoundinfo_beta (1).txt   Compound identifiers and annotations

These filenames match the archived notebook files. From the repository root,
python3 reproduce.py fellowship --work ../p59-run downloads the pinned inputs
under their canonical names and reads the selected rows of the level-5 matrix.
See fellowship/README.md and verification/source_replay/phases2_5/source_manifest.json.

Broad LINCS2020 release
https://clue.io/data/CMap2020#LINCS2020
