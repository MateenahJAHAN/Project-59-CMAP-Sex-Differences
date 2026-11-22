# Data behind the planted-effect simulation

The frozen simulation and the independent Part B rebuild supplied the same 9,288,066 measured expression values at each of two scales. The primary ASG panel covered 85 cell lines, 115 drugs and 978 genes, with measurements for 9,497 of the 9,775 possible line–drug pairs. All assays were at 24 hours. The primary dose was the provider's nominal 1.11 µM label.

| Input | What reproduced |
|---|---|
| Level 3 expression | Subtracting the same plate's mean vehicle expression reproduced all 33,564 treated-well vectors across both doses. All 154 plates had at least two controls, with 2,033 vehicle wells in total. The 17,197 primary-dose wells reproduced every line–drug mean and count. |
| Level 5 expression | Of 18,892 consensus signatures in the two-dose panel, 9,724 carried the primary-dose label. Averaging those signatures reproduced the Level 5 means and all counts. After matching gene order and the model's float32 storage precision, every model input value agreed. |
| Cohort and covariates | Applying the frozen eligibility rules to the rebuilt labels reproduced 41 lines across 11 lineages. Ten expressed only the Y marker, 13 only XIST, and 18 neither. None expressed both. The observed Y, XIST and TP53 design columns, masks and per-lineage group counts matched. |
| Random assignments | Both inputs generated the same 1,000 within-lineage assignments, the same 2,000 Y and XIST designs and the same permutation keys. The simulation uses the first 2,000 of 10,000 precomputed keys. Unplanted response sizes, template norms, and the 29 strongest and 29 weakest compound assignments agreed at both scales. |

Two numerical differences were traceable to storage. Before float32 rounding, 111,372 Level 5 means differed from the stored model matrix by at most 4.77 × 10⁻⁷. Reading the covariate file at round-trip precision exposed seven ploidy differences, each at most 4.45 × 10⁻¹⁶. The simulation adjusts for TP53 and lineage. Its design matrices agreed exactly. The frozen loader also reads a Level 5 vehicle vector for each cell line. Those vectors matched, although the planted-effect analysis uses the treatment signatures directly at Level 5.

## Reproduce the comparison

Use the source-rebuilt Part B workspace, the corresponding checkpoint manifest and frozen task G inputs. The run root contains the archived `power/scripts` and rebuilt `part_b/run/{D,E,D_labels}` directories. The task G checkout supplies the registered model and configuration. Python 3.11 or later, NumPy, pandas, SciPy and the frozen engine's dependencies are required. For an exact replay, match the library versions recorded in the receipt.

```shell
PYTHONHASHSEED=0 python verify.py \
  --run-root "$REPRODUCTION_ROOT" \
  --engine "$TASK_G_CHECKOUT" \
  --frozen-g "$FROZEN_G_INPUTS" \
  --checkpoint-manifest "$CHECKPOINT_MANIFEST" \
  --output receipt.json
```

The [aggregate receipt](receipt.json) records file digests, dimensions, counts and maximum numerical differences. It establishes the inputs to Figure S3. The power estimates themselves require the separate planted-effect run.
