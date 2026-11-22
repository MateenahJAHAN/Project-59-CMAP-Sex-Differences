# Replaying the 50-compound CMap–GTEx screens

All four phases are complete. The categories reproduced for every compound in Phase 2 and for every tissue comparison in Phase 5. Two compounds changed Phase 3 category. Phase 4 retained most archived hits, with additional gains and losses at the enrichment threshold. [Full results and corrected figures](RESULTS.md) · [Per-compound outputs](outputs/)

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

For the changing calls, run `python literal_phase3.py`, `python literal_phase3.py --drug DL-PDMP` and `python literal_phase4.py`. `python verify_outputs.py` audits all 200 per-compound files and the three figure tables.
