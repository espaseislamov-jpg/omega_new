# omega_new

Omega is a chromatogram-processing project for calculating Omega metrics from instrument CSV exports and manually reviewed reference workbooks.

## Conda environment

Use the Conda environment before running the GUI or regression scripts. The environment uses the conda-forge scientific stack; ChromatoPy is installed separately with `--no-deps` so it does not replace the conda-managed numerical packages.

```bash
bash scripts/setup_conda_env.sh
```

Manual equivalent:

```bash
conda env update --name omega --file environment.yml --prune
conda run -n omega python -m pip install --no-deps -r requirements-chromatopy.txt
conda run -n omega python scripts/verify_environment.py
```

## Run the GUI

```bash
conda run -n omega python New_idea.py
```


## Extract uploaded archive data

The uploaded `Desktop.part1.rar` / `Desktop.part2.rar` files are RAR5 split archives. Use `unar`; `unrar-free` may extract `02072026.CSV` only partially.

```bash
apt-get update && apt-get install -y unar
bash scripts/extract_desktop_archives.sh
```

After extraction, the regression harness will find `extracted_desktop/test_bigbatch_020726.xlsx`, `extracted_desktop/02072026.CSV`, `extracted_desktop/test_bigbatch_03072026.xlsx`, and `extracted_desktop/03072026.CSV` recursively.

## Run regression

Current production-like engine:

```bash
conda run -n omega python omega_regression.py --data-dir . --out regression_outputs/omega_regression_current.xlsx
```

Research-only oracle selector over available processing variants. This uses the manual reference to pick the closest variant and should only be used to estimate improvement potential, not for patient reporting:

```bash
conda run -n omega python omega_regression.py --data-dir . --selector-mode oracle --out regression_outputs/omega_regression_oracle.xlsx
```

Optional debug export for outliers. This writes CSV tables plus `plot.png` for each sample above the threshold so the baseline, integration boundaries, and matched targets can be reviewed visually:

```bash
conda run -n omega python omega_regression.py --data-dir . --debug-dir regression_debug --debug-threshold 0.5 --out regression_outputs/omega_regression_current.xlsx
```

## Regression outputs

The regression harness now writes:

- `omega_regression_current.xlsx` with Results, Summary, Outliers, Input_audit, Errors, and Variants sheets.
- `omega_regression_summary.csv` with aggregate metrics.
- `omega_regression_input_audit.csv` with reference/instrument matching counts.
- `omega_regression_errors.csv` with skipped or failed samples.
- `omega_regression_variants.csv` when candidate variants are evaluated.
- `omega_regression_report.md` for human review.

Large local archives, extracted data folders, and debug folders are intentionally ignored by Git.
