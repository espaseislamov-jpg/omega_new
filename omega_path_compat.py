from __future__ import annotations

import os
import sys
from pathlib import Path


def _contains_non_ascii(value: str) -> bool:
    try:
        value.encode("ascii")
    except UnicodeEncodeError:
        return True
    return False


def get_windows_short_path(path: str | os.PathLike[str]) -> Path:
    """Return an existing path through its ASCII-safe Windows 8.3 alias."""
    original = Path(path).expanduser()
    if os.name != "nt" or not original.exists():
        return original

    try:
        import ctypes

        kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
        get_short_path = kernel32.GetShortPathNameW
        get_short_path.argtypes = [
            ctypes.c_wchar_p,
            ctypes.c_wchar_p,
            ctypes.c_uint32,
        ]
        get_short_path.restype = ctypes.c_uint32

        required = get_short_path(str(original), None, 0)
        if required == 0:
            return original
        buffer = ctypes.create_unicode_buffer(required)
        if get_short_path(str(original), buffer, required) == 0:
            return original
        return Path(buffer.value)
    except (AttributeError, ImportError, OSError, ValueError):
        return original


def native_safe_path(path: str | os.PathLike[str]) -> Path:
    """Use a short alias only when a native Windows library may dislike Unicode."""
    original = Path(path)
    if os.name == "nt" and _contains_non_ascii(str(original)):
        return get_windows_short_path(original)
    return original


def configure_windows_path_compat() -> Path | None:
    """Route caches and temporary files through an ASCII-only Windows path."""
    if os.name != "nt":
        return None

    temp_value = os.environ.get("TEMP") or os.environ.get("TMP")
    if temp_value:
        safe_temp = get_windows_short_path(temp_value)
        if safe_temp.exists():
            os.environ["TEMP"] = str(safe_temp)
            os.environ["TMP"] = str(safe_temp)

    local_value = os.environ.get("LOCALAPPDATA") or temp_value
    if not local_value:
        return None

    safe_local = get_windows_short_path(local_value)
    cache_root = safe_local / "Omega" / "cache"
    try:
        cache_root.mkdir(parents=True, exist_ok=True)
    except OSError:
        if not temp_value:
            return None
        cache_root = get_windows_short_path(temp_value) / "Omega" / "cache"
        cache_root.mkdir(parents=True, exist_ok=True)

    cache_paths = {
        "MPLCONFIGDIR": cache_root / "matplotlib",
        "JOBLIB_TEMP_FOLDER": cache_root / "joblib",
    }
    for variable, directory in cache_paths.items():
        directory.mkdir(parents=True, exist_ok=True)
        os.environ.setdefault(variable, str(directory))

    app_dir = (
        Path(sys.executable).resolve().parent
        if getattr(sys, "frozen", False)
        else Path(__file__).resolve().parent
    )
    os.environ.setdefault("OMEGA_NATIVE_APP_DIR", str(native_safe_path(app_dir)))
    return cache_root
