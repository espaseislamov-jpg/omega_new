# Neural integration-quality judge

The quality judge is a small multitask neural network. It does not change peak
boundaries and does not predict the manual omega value. It estimates only:

- expected absolute error of the current integrator;
- probability that the absolute error is greater than 0.3;
- probability that the absolute error is greater than 0.5.

## Leakage protection

All input features are available from the production integrator before a manual
answer is known. Sample IDs, manual values, errors, and batch identifiers are
excluded from the neural-network input. Validation is grouped by complete batch
date, so neighbouring chromatograms never appear on both sides of a train/test
split.

`14072026` is hard-coded as a sealed date in both export and training scripts.
The exporter rejects an explicit request for that date, and the trainer stops if
the date is found in an input dataset. Its manual workbook must not be committed.

## Cloud training

The `Train neural quality judge` GitHub Actions workflow performs the expensive
work on a GitHub-hosted runner:

1. Recalculate every labeled chromatogram sequentially with the current engine.
2. Export numeric peak geometry and integration-decision features.
3. Hold out each error-bearing batch in turn, together with several normal
   batches. This prevents an aggregate score from merely identifying one bad
   date instead of learning transferable integration geometry.
4. Train the final model on all non-sealed labeled batches.
5. Upload the dataset, out-of-fold predictions, validation report, and NumPy
   model artifact.

The workflow runs automatically when relevant training code or regression data
is pushed to the `work` branch. It can also be started manually from the Actions
tab. The produced artifact is not automatically accepted into the application:
the grouped validation report must be reviewed first.

## Lightweight local smoke test

Do not run a complete export on a laboratory computer. To validate the wiring on
one sample only:

```powershell
python scripts/export_quality_dataset.py --dates 06072026 --limit 1 --output artifacts/smoke.csv
```

The training script limits PyTorch to two CPU threads even when launched locally.
