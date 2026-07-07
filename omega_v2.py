from __future__ import annotations

import importlib.util
import os
import runpy
import subprocess
import sys
from pathlib import Path

APP_NAME = "omega_v2"
MIN_PYTHON = (3, 10)
REQUIREMENT_FILES = ("requirements.txt", "requirements-chromatopy.txt")
REQUIRED_MODULES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "matplotlib": "matplotlib",
    "openpyxl": "openpyxl",
    "pybaselines": "pybaselines",
    "lmfit": "lmfit",
    "sklearn": "scikit-learn",
    "hdbscan": "hdbscan",
    "chromatopy": "chromatopy",
}
OPTIONAL_MODULES = {
    "pyopenms": "pyopenms",
}


def _app_dir() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parent


def _resource_path(relative: str) -> Path:
    base = Path(getattr(sys, "_MEIPASS", _app_dir()))
    return base / relative


def _missing_modules() -> list[str]:
    missing: list[str] = []
    for module_name, package_name in {**REQUIRED_MODULES, **OPTIONAL_MODULES}.items():
        if importlib.util.find_spec(module_name) is None:
            missing.append(package_name)
    return missing


def _run_pip(args: list[str]) -> None:
    command = [sys.executable, "-m", "pip", *args]
    print("[omega_v2]", " ".join(command), flush=True)
    subprocess.check_call(command)


def _install_requirements_if_needed() -> None:
    if getattr(sys, "frozen", False):
        return
    missing = _missing_modules()
    if not missing:
        return

    print("[omega_v2] Missing Python packages:", ", ".join(missing), flush=True)
    _run_pip(["install", "--upgrade", "pip"])
    for requirements_name in REQUIREMENT_FILES:
        requirements_path = _resource_path(requirements_name)
        if requirements_path.exists():
            _run_pip(["install", "-r", str(requirements_path)])

    still_missing = _missing_modules()
    if still_missing:
        raise RuntimeError(
            "Could not install required packages automatically: "
            + ", ".join(still_missing)
        )


def main() -> int:
    if sys.version_info < MIN_PYTHON:
        raise RuntimeError(
            f"{APP_NAME} requires Python {MIN_PYTHON[0]}.{MIN_PYTHON[1]}+; "
            f"current Python is {sys.version.split()[0]}"
        )

    os.environ.setdefault("OMEGA_APP_NAME", APP_NAME)
    _install_requirements_if_needed()
    runpy.run_module("New_idea", run_name="__main__")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
