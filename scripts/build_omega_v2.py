from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
SPEC_PATH = PROJECT_DIR / "omega_v2.spec"
REQUIREMENT_FILES = [PROJECT_DIR / "requirements.txt", PROJECT_DIR / "requirements-chromatopy.txt"]


def run(command: list[str]) -> None:
    print("[build_omega_v2]", " ".join(command), flush=True)
    subprocess.check_call(command, cwd=PROJECT_DIR)


def install_runtime_dependencies() -> None:
    run([sys.executable, "-m", "pip", "install", "--upgrade", "pip"])
    for requirements_path in REQUIREMENT_FILES:
        if requirements_path.exists():
            run([sys.executable, "-m", "pip", "install", "-r", str(requirements_path)])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build the omega_v2 PyInstaller app.")
    parser.add_argument("--skip-deps", action="store_true", help="Do not install runtime dependencies before building.")
    args = parser.parse_args(argv)

    if not args.skip_deps:
        install_runtime_dependencies()
    if importlib.util.find_spec("PyInstaller") is None:
        run([sys.executable, "-m", "pip", "install", "pyinstaller>=6,<7"])
    run([sys.executable, "-m", "PyInstaller", "--clean", "--noconfirm", str(SPEC_PATH)])
    print("\nBuild finished.")
    print("Windows: dist/omega_v2/omega_v2.exe")
    print("Linux/macOS: dist/omega_v2/omega_v2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
