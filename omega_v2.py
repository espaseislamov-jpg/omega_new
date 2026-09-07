from __future__ import annotations

import importlib.util
import os
import runpy
import subprocess
import sys
import traceback
from pathlib import Path

from omega_path_compat import configure_windows_path_compat

APP_NAME = "omega_v2"
MIN_PYTHON = (3, 10)
REQUIREMENT_FILES = ("requirements.txt",)
REQUIRED_MODULES = {
    "numpy": "numpy",
    "pandas": "pandas",
    "scipy": "scipy",
    "matplotlib": "matplotlib",
    "openpyxl": "openpyxl",
    "pybaselines": "pybaselines",
    "lmfit": "lmfit",
}
OPTIONAL_MODULES = {
    "pyopenms": "pyopenms",
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
    for module_name, package_name in REQUIRED_MODULES.items():
        if importlib.util.find_spec(module_name) is None:
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
    if sys.version_info < MIN_PYTHON:
        raise RuntimeError(
            f"{APP_NAME} requires Python {MIN_PYTHON[0]}.{MIN_PYTHON[1]}+; "
            f"current Python is {sys.version.split()[0]}"
    )

    os.environ.setdefault("OMEGA_APP_NAME", APP_NAME)
    configure_windows_path_compat()
    _install_requirements_if_needed()
    if "--smoke-test" in sys.argv:
        import tkinter as tk
        from New_idea import ChromatogramApp, process_chromatogram_batch
        import omega_core
        root = tk.Tk()
        try:
            app = ChromatogramApp(root)
            root.update()
            if "--smoke-csv" in sys.argv:
                csv_path = Path(sys.argv[sys.argv.index("--smoke-csv") + 1])
                batches = omega_core.load_batches(csv_path)
                if not batches:
                    raise RuntimeError("Smoke test CSV contains no batches")
                result = process_chromatogram_batch(batches[0]["dataframe"], app.reference_targets)
                if result["matched_targets_df"].empty:
                    raise RuntimeError("Smoke test produced no target rows")
            print("OMEGA_SMOKE_OK", flush=True)
        finally:
            root.destroy()
    else:
        runpy.run_module("New_idea", run_name="__main__")
    return 0


if __name__ == "__main__":
    # A windowed executable has no console. Keep startup failures observable.
    log_root = Path(os.environ.get("LOCALAPPDATA", str(_app_dir()))) / "Omega" / "logs"
    log_root.mkdir(parents=True, exist_ok=True)
    log_path = log_root / "startup.log"
    with log_path.open("w", encoding="utf-8", buffering=1) as log:
        old_stdout, old_stderr = sys.stdout, sys.stderr
        if getattr(sys, "frozen", False):
            sys.stdout = sys.stderr = log
        try:
            print("Starting Omega", sys.executable, flush=True)
            exit_code = main()
        except Exception:
            traceback.print_exc(file=log)
            log.flush()
            if "--smoke-test" not in sys.argv:
                try:
                    import ctypes
                    ctypes.windll.user32.MessageBoxW(None, "Не удалось запустить Omega. Подробности:\n" + str(log_path), "Ошибка запуска Omega", 16)
                except Exception:
                    pass
            exit_code = 1
        finally:
            sys.stdout, sys.stderr = old_stdout, old_stderr
    raise SystemExit(exit_code)
