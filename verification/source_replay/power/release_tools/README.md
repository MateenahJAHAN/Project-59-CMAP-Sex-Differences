# Verify the full power simulation

To recalculate Figure S3 from the public assignment records, use the [offline replay](../replay_data/). This directory contains the checks used to produce those records from the original simulation.

The full check needs four inputs.

- `study` is the fresh sensitivity run with `raw/` and `results/`.
- `original` is the archived sensitivity run with the same layout.
- `calibration` contains the registered model’s `G_calibration.tsv`.
- `checkpoints` is the fresh simulation directory containing `run.json` and the 1,000 individual assignment files.

With Python 3.11 or later, NumPy, pandas and SciPy, run these commands from this directory after setting the four paths. The first command checks the saved assignments against the assembled output and compares all scientific fields with the archive. The second independently calculates the figure endpoints, intervals, coverage and calibration. The last two build the small public dataset and repeat its checks using the Python standard library.

```sh
python3 -B package.py --run "$study" --reference "$original" --calibration "$calibration" --checkpoints "$checkpoints" --out full_run
python3 -B audit.py --run "$study" --calibration "$calibration" --packet full_run
python3 -B build_compact.py --run "$study" --final full_run --checkpoints "$checkpoints" --out replay_data
python3 -I -S -B replay_data/replay_s3.py
```

[Released aggregate results](../full_run/) · [Numerical checks](../full_run/comparison.json)
