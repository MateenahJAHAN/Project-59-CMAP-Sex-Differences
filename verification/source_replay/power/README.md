# Sensitivity to planted effects

In the 41-line LINCS panel, increasing every response vector in the labelled group by 20% led the panel test to reject in 59 of 1,000 random assignments for Y and 28 for XIST. The 95% interval for the pooled, vehicle-referenced response excluded zero in 301 and 396 assignments, respectively. These rates describe sensitivity to that particular simulated change. [Figure S3](../../../manuscript/figures/figS3_power.pdf) · [Plotted values](recount/figS3_power_values.tsv)

The fresh simulation used 2,000 within-lineage permutations per test and reproduced all 26,000 scenario records, 120,000 planted-compound records and 5,750,000 interval entries. [Assignment-level comparison](checks/checkpoint_comparison.json)

Across the 4,000 simulated cases with one planted compound, 34 planted compounds reached BH q < 0.1. Four passed despite two to four permutation exceedances because the adjustment used all tested compounds. Detection rates at 10,000 permutations were not measured. [Four cases and matched records](checks/single_plant_comparison.json) · [Supplement, Note S17](../../../manuscript/supplement.pdf)

The source rebuild reproduced all 9,288,066 measured expression values at each scale, the 1,000 within-lineage assignments, and the 2,000 marker-specific model designs. [Input comparison](../power_input_check/)

[Independent recount of the archived assignments](recount_comparison.json) · [48 implementation checks](checks/harness_checks.tsv) · [Restart check](checks/checkpoint_restore.json)

To repeat the simulation, set `P59_G_REPO` to the task-G analysis checkout and `P59_G_RUN` to the production input directory containing `inputs/data/G` and `capture`. The scripts use the production model, design and original random seeds. From this directory run the following in the study Python environment.

```sh
export PYTHONHASHSEED=0
mkdir -p raw results logs
python -I scripts/check_harness.py
python -I scripts/power_checkpoint.py --B 2000 --R 1000 --workers 6 --tag power
python -I scripts/summarize_power.py --tag power
```
