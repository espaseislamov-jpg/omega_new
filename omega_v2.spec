# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_submodules, copy_metadata, collect_data_files
from pathlib import Path
import hashlib
import json
import runpy

release = runpy.run_path('omega_version.py')

digest = hashlib.sha256()
for path in sorted(Path("omega_core").glob("*.py")):
    digest.update(path.name.encode())
    digest.update(path.read_text(encoding="utf-8").encode())
manifest = Path("build/runtime_manifest.json")
manifest.parent.mkdir(exist_ok=True)
manifest.write_text(json.dumps({"engine_sha256": digest.hexdigest()}), encoding="utf-8")
runtime_data = [(str(manifest), "."), ("requirements.txt", ".")]
for package in ("numpy", "pandas", "scipy", "pybaselines", "lmfit", "pyopenms"):
    runtime_data += copy_metadata(package)
runtime_data += collect_data_files("pyopenms")

hiddenimports = [
    "New_idea",
    "omega_path_compat",
    "matplotlib.backends.backend_tkagg",
    "scipy.optimize",
    "scipy.signal",
    "scipy.stats",
    "scipy.integrate",
    "openpyxl",
    "pybaselines",
    "lmfit",
    "pyopenms",
]
hiddenimports += collect_submodules("omega_core")

excludes = [
    "IPython",
    "jedi",
    "pytest",
    "tkinter.test",
    "numpy.tests",
    "pandas.tests",
    "scipy.tests",
    "matplotlib.tests",
    "PyQt5",
    "PyQt6",
    "PySide2",
    "PySide6",
]

block_cipher = None

a = Analysis(
    ["omega_v2.py"],
    pathex=[],
    binaries=[],
    datas=[
        ("reference_targets_reverted_c22fixed.json", "."),
        ("chebyshev_coefficients.csv", "."),
        ("omega_default_profiles.json", "."),
    ] + runtime_data,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=excludes,
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="Omega",
    version="installer/version_info.txt",
    contents_directory=release['RUNTIME_DIR'],
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="Omega_V3.0",
)
