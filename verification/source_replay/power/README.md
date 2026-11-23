# Sensitivity to planted effects

In the 41-line LINCS panel, increasing every response vector in the labelled group by 20% led the panel test to reject in 59 of 1,000 random assignments for Y and 28 for XIST. The 95% interval for the pooled, vehicle-referenced response excluded zero in 301 and 396 assignments, respectively. These rates describe sensitivity to that particular simulated change. [Figure S3](../../../manuscript/figures/figS3_power.pdf) · [Plotted values](full_run/figS3_power_values.tsv)

The fresh simulation used 2,000 within-lineage permutations per test and reproduced all 26,000 scenario records, 120,000 planted-compound records and 5,750,000 interval entries. [Assignment-level comparison](checks/checkpoint_comparison.json)

Across the 4,000 simulated cases with one planted compound, 34 planted compounds reached BH q < 0.1. Four passed despite two to four permutation exceedances because the adjustment used all tested compounds. Detection rates at 10,000 permutations were not measured. [Four cases and matched records](checks/single_plant_comparison.json) · [Supplement, Note S17](../../../manuscript/supplement.pdf)

The separate calibration check reran 200 of the same assignments under two alternative designs for each marker. All 800 results and six summary rows matched. [Calibration results](full_run/size_diagnostic.tsv) · [Full comparison with the archive](full_run/comparison.json)

### Recalculate the figure

Download the [assignment records and verifier](replay_data/) together. With Python 3.11 or later, run the following from this directory. It recalculates all 66 plotted points and their interval endpoints, then checks all 4,000 single-compound cases.

```sh
python3 -I -S -B replay_data/replay_s3.py
```

[Power estimates, confidence-interval coverage and calibration](full_run/) · [Full-run verification code](release_tools/)

### Repeat the simulation

The source rebuild reproduced all 9,288,066 measured expression values at each scale, the 1,000 within-lineage assignments, and the 2,000 marker-specific model designs. [Input comparison](../power_input_check/)

[Independent recount of the archived assignments](recount_comparison.json) · [48 implementation checks](checks/harness_checks.tsv) · [Restart check](checks/checkpoint_restore.json)

From the repository root, complete the LINCS model and then run the sensitivity analysis on the same reconstructed tensors and labels.

```sh
python3 reproduce.py part-b --work ../p59-run
python3 reproduce.py power --work ../p59-run --workers 6
```

The simulation writes atomic checkpoints for the 1,000 assignments in the work directory. Repeating the same command on the same inputs resumes completed assignments. After the simulation, the runner compares the new power, coverage, calibration and Figure S3 tables with the release. The separate leave-one-line-out analysis uses the three flags reconstructed by the `transfer` stage. When run with Part B alone, it uses the identical [released sensitivity list](sensitivity_lines.tsv) and records that source.

If the assignment files are already complete and summarization was interrupted, run `verification/source_replay/power/run.py --work ../p59-run --summarize-saved` with the Part B interpreter. It reuses the scenario, planted-compound and interval files under `power/raw/`, reruns calibration diagnostics and compares all nine numerical output tables.
