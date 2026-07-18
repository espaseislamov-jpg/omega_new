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



## Repository layout

- `omega_core/` — active modular processing engine.
- `New_idea.py` — current GUI entry point.
- `data/regression/` — committed historical regression CSV/XLSX pairs.
- `data/archives/` — committed split archive parts used to recover the newest large batches.
- `regression_outputs/` — local generated regression reports (ignored by Git).
- `docs/` — investigation notes, diagnostics, and integration-plan documents.

## Extract uploaded archive data

The uploaded `data/archives/Desktop.part1.rar` / `data/archives/Desktop.part2.rar` files are RAR5 split archives. Use `unar`; `unrar-free` may extract `02072026.CSV` only partially.

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

## Build omega_v2 executable

`omega_v2.py` is the bootstrap launcher for the GUI. When it is started from source on a machine with Python 3.10+, it checks the required Python packages and installs missing dependencies from `requirements.txt` before launching `New_idea.py`. Experimental judge-training tools use the separate `requirements-training.txt` file and are not bundled into the desktop application.

Build the distributable app with PyInstaller:

```bash
python scripts/build_omega_v2.py
```

On Windows, run:

```bat
scripts\build_omega_v2_windows.bat
```

The expected Windows output is:

```text
dist\omega_v2\omega_v2.exe
```

Note: PyInstaller builds for the current operating system. A true Windows `.exe` must be built on Windows; Linux produces `dist/omega_v2/omega_v2`.

## Download omega_v2 setup from GitHub

For a non-programmer workflow, use the GitHub Actions artifact:

1. Open the repository on GitHub.
2. Go to **Actions**.
3. Open **Build omega_v2 Windows app**.
4. Press **Run workflow** on branch `work`.
5. Wait until the run becomes green.
6. Open the finished run and download artifact **omega_v2_windows_setup**.
7. Unzip it and run `omega_v2.6_setup.exe`.

The workflow also uploads **omega_v2_windows_portable** if you want the raw portable folder without installer.
