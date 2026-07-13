# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_submodules

hiddenimports = [
    "New_idea",
    "omega_chromatopy_clean",
    "omega_path_compat",
    "matplotlib.backends.backend_tkagg",
    "scipy.optimize",
    "scipy.signal",
    "scipy.stats",
    "scipy.integrate",
    "openpyxl",
    "pybaselines",
    "lmfit",
    "sklearn",
    "hdbscan",
    "chromatopy",
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
    "sklearn.tests",
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
        ("requirements.txt", "."),
        ("requirements-chromatopy.txt", "."),
    ],
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
pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="omega_v2",
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
    a.zipfiles,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name="omega_v2",
)
