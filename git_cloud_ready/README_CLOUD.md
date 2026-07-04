# Omega: cloud development bundle

This directory contains the source and support files required to continue
Omega development in Git. Generated installers, virtual environments, caches,
plots, and spreadsheet reports are intentionally excluded.

## Environment

Use Python 3.13.

### Windows

```powershell
py -3.13 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install --no-deps -r requirements-chromatopy.txt
python New_idea.py
```

### Linux / cloud workspace

Tkinter may need to be installed through the operating-system package manager
before creating the environment.

```bash
python3.13 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python -m pip install --no-deps -r requirements-chromatopy.txt
python New_idea.py
```

The GUI requires a graphical desktop. Headless workspaces can use
`omega_core`, `omega_regression.py`, and the experimental command-line scripts
without opening the GUI.

## Main files

- `New_idea.py`: current Tkinter application.
- `omega_chromatopy_clean.py`: processing engine currently used by the GUI.
- `omega_core/`: modular alternative/legacy production pipeline.
- `Omega_cluster_engine.py`: experimental cluster-deconvolution engine.
- `omega_regression.py`: regression harness.
- `reference_targets_reverted_c22fixed.json`: fatty-acid target definitions.
- `omega_path_compat.py`: Windows Cyrillic-path compatibility bootstrap.
- `test_data/`: one compact diagnostic snapshot for smoke testing.

## ChromatoPy compatibility

ChromatoPy 2.0.0 pins `hdbscan==0.8.40`. That release has no ready-made
Python 3.13 wheel on Windows and requires Microsoft C++ Build Tools.
The tested workaround is to install the compatible dependency set first and
then install ChromatoPy with `--no-deps`, as shown above.

## Regression data

The full regression corpus is not present in the original project directory.
`omega_regression.py` currently defaults to:

```text
C:\Users\marat\Desktop\CSV_Omega
```

Pass the real corpus explicitly when it becomes available:

```powershell
python omega_regression.py --data-dir PATH_TO_CORPUS
```

