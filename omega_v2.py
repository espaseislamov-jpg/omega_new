from __future__ import annotations

import importlib.util
import importlib.metadata
import os
import runpy
import subprocess
import sys
import logging
import json
from logging.handlers import RotatingFileHandler
from pathlib import Path

from omega_path_compat import configure_windows_path_compat
from omega_version import APP_NAME

MIN_PYTHON = (3, 14)
REQUIREMENT_FILES = ("requirements.txt",)
REQUIRED_MODULES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "matplotlib": "matplotlib",
    "openpyxl": "openpyxl",
    "pybaselines": "pybaselines",
    "lmfit": "lmfit",
    "pyopenms": "pyopenms",
}
OPTIONAL_MODULES = {
    "chromatopy": "chromatopy",
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
    path = _resource_path("requirements.txt")
    pins = dict(line.strip().split("==", 1) for line in path.read_text(encoding="utf-8").splitlines()
                if "==" in line and not line.lstrip().startswith("#")) if path.exists() else {}
    for module_name, package_name in REQUIRED_MODULES.items():
        if importlib.util.find_spec(module_name) is None:
            missing.append(package_name)
        elif package_name in pins:
            try:
                if importlib.metadata.version(package_name) != pins[package_name]:
                    missing.append(package_name)
            except importlib.metadata.PackageNotFoundError:
                missing.append(package_name)
    return missing


def _missing_optional_modules() -> list[str]:
    missing: list[str] = []
    for module_name, package_name in OPTIONAL_MODULES.items():
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

    optional_missing = _missing_optional_modules()
    if optional_missing:
        print(
            "[omega_v2] Optional Python packages are not installed: "
            + ", ".join(optional_missing),
            flush=True,
        )


def main() -> int:
    log_dir = Path(os.environ.get("LOCALAPPDATA", str(Path.home()))) / "Omega" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    handler = RotatingFileHandler(log_dir / "startup.log", maxBytes=2_000_000, backupCount=2, encoding="utf-8")
    logging.getLogger("omega").addHandler(handler)
    logging.getLogger("omega").setLevel(logging.INFO)
    logging.getLogger("omega").info("Starting Omega from %s", sys.executable)
    if sys.version_info < MIN_PYTHON:
        raise RuntimeError(
            f"{APP_NAME} requires Python {MIN_PYTHON[0]}.{MIN_PYTHON[1]}+; "
            f"current Python is {sys.version.split()[0]}"
    )

    os.environ.setdefault("OMEGA_APP_NAME", APP_NAME)
    configure_windows_path_compat()
    _install_requirements_if_needed()
    try:
        from omega_core import runtime
        runtime.verify_versions()
        (log_dir / "runtime.json").write_text(json.dumps(runtime.snapshot(), indent=2), encoding="utf-8")
        if "--verify-csv" in sys.argv:
            from omega_core.verification import main as verify
            return verify(sys.argv[1:])
        runpy.run_module("New_idea", run_name="__main__")
    except Exception:
        logging.getLogger("omega").exception("Application startup failed")
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
