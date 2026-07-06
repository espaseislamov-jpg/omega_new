from __future__ import annotations

import importlib
import sys

REQUIRED_IMPORTS = [
    "chromatopy",
    "hdbscan",
    "lmfit",
    "matplotlib",
    "numpy",
    "openpyxl",
    "pandas",
    "pybaselines",
    "pyopenms",
    "scipy",
    "sklearn",
]


def main() -> int:
    missing: list[str] = []
    for module_name in REQUIRED_IMPORTS:
        try:
            importlib.import_module(module_name)
        except Exception as exc:  # pragma: no cover - diagnostic script
            missing.append(f"{module_name}: {type(exc).__name__}: {exc}")
    if missing:
        print("Omega environment check failed:", file=sys.stderr)
        for item in missing:
            print(f"- {item}", file=sys.stderr)
        return 1
    print("Omega environment check passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
