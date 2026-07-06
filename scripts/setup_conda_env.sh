#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${OMEGA_CONDA_ENV:-omega}"
CONDA_BIN="${CONDA_BIN:-conda}"

if ! command -v "$CONDA_BIN" >/dev/null 2>&1; then
  echo "Conda executable not found: $CONDA_BIN" >&2
  echo "Install Miniconda/Mambaforge or set CONDA_BIN=/path/to/conda." >&2
  exit 1
fi

"$CONDA_BIN" env update --name "$ENV_NAME" --file environment.yml --prune
"$CONDA_BIN" run -n "$ENV_NAME" python -m pip install --no-deps -r requirements-chromatopy.txt
"$CONDA_BIN" run -n "$ENV_NAME" python scripts/verify_environment.py
