# Replaying the 50-compound CMap–GTEx screens

All four phases are complete. The categories reproduced for every compound in Phase 2 and for every tissue comparison in Phase 5. Two compounds changed Phase 3 category. Phase 4 retained most archived hits, with additional gains and losses at the enrichment threshold. [Full results and corrected figures](RESULTS.md) · [Per-compound outputs](outputs/)

## Reproduce

The [fellowship runner](../../../fellowship/README.md#run-the-fellowship-analyses) retrieves and verifies the source data and writes into a separate results directory. Run these commands from the repository root with Python 3.12.

```sh
sh fellowship/install.sh "$PWD/../p59-env"
../p59-env/bin/python fellowship/run.py --data "$PWD/../p59-inputs" \
  --out "$PWD/../p59-results" --phases 2 3 4 5 --workers 2
```

Open `../p59-results/comparison/` for the new analyses and `../p59-results/figures/` for the regenerated plots. The prepared LINCS matrix matched all 450,255,544 float32 values in the [provider GCTX](input_verification/comparison.json), and the [GTEx checksum](input_verification/gtex_comparison.json) matched the source archive.

The GSEA permutation seed was 42 throughout. Phase 3 used 1,000 permutations, Phase 4 used 200 and Phase 5 used 1,000. The input Level 5 signatures are the study's starting data. The assay bead files are upstream of these analyses.

The [reported comparisons](RESULTS.md) include sensitivity analyses for the two changing Phase 3 calls and the Phase 4 NVP-BEZ235 example. `verify_outputs.py` audits the archived 200 per-compound files and the three figure tables.
