# Sensitivity to planted effects

The checked-in recount uses the archived 26,000 scenario rows, 120,000 planted-compound records and saved interval estimates. It reproduces the source values of Figure S3. The independent simulation restarted from assignment zero after the initial process stopped. The restarted driver uses the original analysis and saves the full records for each completed assignment. A separate process restored the first saved assignment and reproduced all 26 scenario records, all 120 planted-compound records and all 5,750 entries of its interval array. [Restart verification](checks/checkpoint_restore.json)

[Recount comparisons](recount_comparison.json) · [Figure S3 values](recount/figS3_power_values.tsv) · [48 harness checks](checks/harness_checks.tsv)

To repeat the simulation, set `P59_G_REPO` to the task-G analysis checkout and `P59_G_RUN` to the production input directory containing `inputs/data/G` and `capture`. The scripts use the production model, design and original random seeds. From this directory run the following in the study Python environment.

```sh
mkdir -p raw results logs
python -I scripts/check_harness.py
python -I scripts/power_checkpoint.py --B 2000 --R 1000 --workers 6 --tag power
python -I scripts/summarize_power.py --tag power
```
