# Replaying the 50-compound CMap–GTEx screens

Phase 2: 50/50 · Phase 3: 32/50 · Phase 4: 37/50 · Phase 5: 50/50. [Completion status](handoff_status.json) lists which full-cohort comparisons have finished. Each compressed CSV in [outputs](outputs/) contains the complete set of results for one compound. A partial collection can be resumed with the commands below.

Phase 2 reproduced all 50 published compound labels and counts. AS-605240 liver was significant with adjustment across the 20 tissues with overlap (q = 0.02471) and fell above threshold after adjustment over all 44 eligible tissues (q = 0.05437). Zero of 50 compounds had a hit under the latter rule. [All tissue tests](comparison/phase2_all_tissue_tests.csv)

Phase 5 reproduced all 39,996 compound–lineage–tissue categories, including the 620 tested comparisons. All 620 NES and FDR values matched using GSEApy 1.3.1. Version 1.2.1 changed 14 FDR values by at most 0.001 and left the calls intact. The alternate within-lineage panel in cell 60 reproduced the 1 A, 4 B and 45 C calls for all 50 archived names. Its selection tie at rank 50 admits either ZM-336372 or erastin. [Versions](comparison/phase5_software_versions.csv) · [Workbook checks](comparison/auxiliary_exports.json) · [Alternative panel](phase5_secondary/comparison.json) · [Phase 5 figure](figures/phase5_testable_tissues.png)

Fostamatinib's Phase 3 home-tissue result was sensitive to the ordering of exactly tied ranks. Applying the original notebook functions and cmapPy reader yielded FDR 0.0410 for the blood/lymphoid lineage against EBV-transformed lymphocytes. Reordering only exact ties gave 0.0336, 0.0757 and 0.0745. The two stored individual scores for this compound apply to different tissues. [All home-set scores](phase3_literal/home_tissue_scores.csv)

## Reproduce

Use Python 3.12 and the [pinned packages](requirements.txt). Set P59_PUBLIC to a checkout at the commit in [the source receipt](source_manifest.json). Set P59_META to the source directory with geneinfo_beta.txt, siginfo_beta.txt and cellinfo_beta.txt. Set P59_GTEX to the directory with the original GTEx signif.sbgenes.txt. Set P59_MATRIX to a directory with sig50.npy (signatures by genes), sig50.sigs.txt and row_ids.txt. The staged matrix used for these runs matched all 450,255,544 float32 values from the [provider GCTX](input_verification/comparison.json) and the [GTEx checksum](input_verification/gtex_comparison.json) matched the source archive.

```sh
export PYTHONHASHSEED=0 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 RAYON_NUM_THREADS=1
python replay.py --prepare
python replay.py --phases 2 3 4 5 --workers 2
python compare.py
python figures.py
```

The GSEA permutation seed was 42 throughout. Phase 3 used 1,000 permutations, Phase 4 used 200 and Phase 5 used 1,000. The input Level 5 signatures are the study's starting data. The assay bead files are upstream of these analyses.
