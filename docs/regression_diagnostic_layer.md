# Regression diagnostic layer

The diagnostic layer is designed to explain **why** a calculated omega value misses
the manual reference, not to change the calculation itself.

## Generated outputs

A normal regression run now writes these additional artifacts next to the workbook:

- `omega_regression_sample_diagnostics.csv` — one row per evaluated sample with a
  diagnostic bucket and reason tags.
- `omega_regression_peak_diagnostics.csv` — one row per sample/target with RT,
  integration width, left/right width, asymmetry, area, and matching status.
- `omega_regression_issue_summary.csv` — grouped error statistics by diagnostic
  bucket.
- The Excel workbook also contains `Sample_diagnostics`, `Peak_diagnostics`, and
  `Issue_summary` sheets.
- The Markdown report includes a diagnostic issue summary and the highest-error
  diagnostic samples.

## Clinical thresholds

The report keeps two practical bands visible:

- `±0.3` — desirable inter-operator band.
- `±0.5` — maximum clinical tolerance band requested for the current workflow.

Samples outside `±0.5` are classified into directional buckets such as
`over_c22_cluster`, `under_c22_cluster`, or `under_c20_cluster` so the next tuning
step can target a specific mechanism rather than the whole algorithm.

## Reason tags

Reason tags are intentionally simple and inspectable. Examples:

- `c22_complex_boundaries` — C22 targets used tail/base expansion or fitting logic.
- `c20_complex_boundaries` — C20/EPA-neighbor region used local/base-expanded/fitted
  integration boundaries.
- `c18_complex_boundaries` — denominator-side C18 targets used valley/base-expanded
  boundaries.
- `high_dpa_to_c22_4_ratio` — DPA is large relative to C22:4 in the same C22 cluster.
- `c22_debit_applied` / `c22_credit_applied` — omega metric applied a C22 correction.
- `baseline_fallback` — the batch needed ASLS/arPLS fallback instead of the default
  Chebyshev baseline.
- `large_rt_error`, `wide_peak_window`, `asymmetric_peak_window` — target geometry
  is suspicious and should be reviewed on the chromatogram plot.

## How to use it

1. Open `omega_regression_issue_summary.csv` to see which mechanisms dominate the
   errors above `±0.5`.
2. Open `omega_regression_sample_diagnostics.csv` and filter by `diagnostic_bucket`.
3. For a specific sample, inspect matching rows in `omega_regression_peak_diagnostics.csv`.
4. If needed, rerun with `--debug-dir regression_debug --debug-threshold 0.5` to
   generate per-sample plots and matched-target CSVs for every miss above `±0.5`.

The current generated diagnostics show that the remaining misses above `±0.5` are
not random: they concentrate in C22 cluster buckets first, then C20 cluster cases.
