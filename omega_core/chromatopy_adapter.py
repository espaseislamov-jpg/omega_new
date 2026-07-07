from __future__ import annotations

import os
import importlib.metadata
import importlib.util
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from .signal import _get_x_column_name, _robust_sigma


ENABLE_CHROMATOPY_INTEGRATION = os.environ.get("OMEGA_USE_CHROMATOPY_INTEGRATION", "0").strip() == "1"
CHROMATOPY_FIT_MODE = os.environ.get("OMEGA_CHROMATOPY_FIT_MODE", "single").strip().lower()
CHROMATOPY_GAUS_ITERATIONS = int(os.environ.get("OMEGA_CHROMATOPY_GAUS_ITERATIONS", "800"))
CHROMATOPY_PK_SNS = float(os.environ.get("OMEGA_CHROMATOPY_PK_SNS", "0.001"))
CHROMATOPY_SMOOTHING_WINDOW = int(os.environ.get("OMEGA_CHROMATOPY_SMOOTHING_WINDOW", "7"))
CHROMATOPY_SMOOTHING_POLYORDER = int(os.environ.get("OMEGA_CHROMATOPY_SMOOTHING_POLYORDER", "3"))
CHROMATOPY_MATCH_TOLERANCE = float(os.environ.get("OMEGA_CHROMATOPY_MATCH_TOLERANCE", "0.055"))
CHROMATOPY_MAX_WIDTH = float(os.environ.get("OMEGA_CHROMATOPY_MAX_WIDTH", "0.180"))
CHROMATOPY_MIN_AREA_RATIO = float(os.environ.get("OMEGA_CHROMATOPY_MIN_AREA_RATIO", "0.45"))
CHROMATOPY_MAX_AREA_RATIO = float(os.environ.get("OMEGA_CHROMATOPY_MAX_AREA_RATIO", "2.50"))
CHROMATOPY_TARGET_CODES_TEXT = os.environ.get("OMEGA_CHROMATOPY_TARGET_CODES", "C16:1N7,C18:3N6,C20:5,C22:6")
CHROMATOPY_TARGET_CODES = {
    item.strip().upper()
    for item in CHROMATOPY_TARGET_CODES_TEXT.split(",")
    if item.strip()
}

_CHROMATOPY_FUNCTIONS = None


def _chromatopy_fid_module_path() -> Path:
    try:
        distribution = importlib.metadata.distribution("chromatopy")
        module_path = Path(distribution.locate_file("chromatopy/FID/FID_Integration_functions.py"))
        if module_path.exists():
            return module_path
    except importlib.metadata.PackageNotFoundError:
        pass

    package_spec = importlib.util.find_spec("chromatopy")
    if package_spec is not None and package_spec.submodule_search_locations:
        for package_dir in package_spec.submodule_search_locations:
            module_path = Path(package_dir) / "FID" / "FID_Integration_functions.py"
            if module_path.exists():
                return module_path

    raise ImportError(
        "Cannot locate ChromatoPy FID_Integration_functions.py. "
        "Install chromatopy or rebuild omega_v2 with chromatopy bundled."
    )


def _load_chromatopy_functions():
    global _CHROMATOPY_FUNCTIONS
    if _CHROMATOPY_FUNCTIONS is not None:
        return _CHROMATOPY_FUNCTIONS

    try:
        from chromatopy.FID.FID_Integration_functions import fit_gaussians, smoother
    except Exception:
        module_path = _chromatopy_fid_module_path()
        spec = importlib.util.spec_from_file_location("_omega_chromatopy_fid_functions", module_path)
        if spec is None or spec.loader is None:
            raise ImportError(f"Cannot load ChromatoPy FID functions from {module_path}")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        fit_gaussians = module.fit_gaussians
        smoother = module.smoother

    _CHROMATOPY_FUNCTIONS = (fit_gaussians, smoother)
    return _CHROMATOPY_FUNCTIONS


def _build_chromatopy_peak_index(x: np.ndarray, y: np.ndarray, smoother_func) -> tuple[np.ndarray, np.ndarray]:
    smoothing = [CHROMATOPY_SMOOTHING_WINDOW, CHROMATOPY_SMOOTHING_POLYORDER]
    y_smooth = np.asarray(smoother_func(y, smoothing[0], smoothing[1]), dtype=float)
    y_smooth = np.clip(y_smooth, 0.0, None)
    noise = max(_robust_sigma(y), 1.0)
    prominence_floor = max(noise * 1.8, float(np.quantile(y_smooth, 0.75)) * 0.025, 5.0)
    height_floor = max(noise * 1.2, float(np.quantile(y_smooth, 0.60)))
    dx = float(np.median(np.diff(x))) if len(x) > 1 else 0.001
    min_distance = max(1, int(round(0.012 / max(dx, 1e-9))))
    peak_indices, _ = find_peaks(
        y_smooth,
        height=height_floor,
        prominence=prominence_floor,
        distance=min_distance,
    )
    return peak_indices.astype(int), y_smooth


def _nearest_peak_index(x: np.ndarray, peak_indices: np.ndarray, target_rt: float) -> int | None:
    if peak_indices.size == 0 or not np.isfinite(target_rt):
        return None
    distances = np.abs(x[peak_indices] - float(target_rt))
    best_pos = int(np.argmin(distances))
    if float(distances[best_pos]) > CHROMATOPY_MATCH_TOLERANCE:
        return None
    return int(peak_indices[best_pos])


def _choose_neighbor_peaks(x: np.ndarray, peak_indices: np.ndarray, peak_idx: int, max_neighbors: int = 3) -> np.ndarray:
    if peak_indices.size == 0:
        return np.array([peak_idx], dtype=int)
    ordered = sorted(
        {int(idx) for idx in peak_indices.tolist() + [int(peak_idx)]},
        key=lambda idx: abs(float(x[idx]) - float(x[peak_idx])),
    )
    return np.array(sorted(ordered[:max_neighbors]), dtype=int)


def _fit_one_peak(
    x: np.ndarray,
    y_smooth: np.ndarray,
    peak_indices: np.ndarray,
    peak_idx: int,
) -> dict | None:
    fit_gaussians, _ = _load_chromatopy_functions()
    x_series = pd.Series(x)
    y_series = pd.Series(y_smooth)
    neighbors = _choose_neighbor_peaks(x, peak_indices, peak_idx, max_neighbors=3)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit_x, fit_y, _area_smooth, area_ensemble, model = fit_gaussians(
                x_series,
                y_series,
                int(peak_idx),
                neighbors,
                [CHROMATOPY_SMOOTHING_WINDOW, CHROMATOPY_SMOOTHING_POLYORDER],
                CHROMATOPY_PK_SNS,
                gi=CHROMATOPY_GAUS_ITERATIONS,
                mode=CHROMATOPY_FIT_MODE,
            )
    except Exception:
        return None

    fit_x = np.asarray(fit_x, dtype=float)
    fit_y = np.asarray(fit_y, dtype=float)
    area_ensemble = np.asarray(area_ensemble, dtype=float)
    area_ensemble = area_ensemble[np.isfinite(area_ensemble)]
    if fit_x.size < 2 or fit_y.size < 2 or area_ensemble.size == 0:
        return None

    peak_pos = int(np.nanargmax(fit_y))
    area = float(np.nanmedian(area_ensemble))
    if not np.isfinite(area) or area <= 0:
        return None
    return {
        "found_rt": float(fit_x[peak_pos]),
        "area": area,
        "integration_start_x": float(np.nanmin(fit_x)),
        "integration_end_x": float(np.nanmax(fit_x)),
        "model": str(model.get("name", CHROMATOPY_FIT_MODE)) if isinstance(model, dict) else CHROMATOPY_FIT_MODE,
    }


def apply_chromatopy_target_integration(
    processed: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets.copy()
    if not ENABLE_CHROMATOPY_INTEGRATION or processed is None or processed.empty or out is None or out.empty:
        return out

    try:
        _, smoother_func = _load_chromatopy_functions()
    except Exception:
        return out

    x_col = _get_x_column_name(processed)
    x = processed[x_col].to_numpy(dtype=float)
    y = np.clip(processed["y_corrected"].to_numpy(dtype=float), 0.0, None)
    if len(x) < 16 or not np.any(y > 0):
        return out

    peak_indices, y_smooth = _build_chromatopy_peak_index(x, y, smoother_func)
    if peak_indices.size == 0:
        return out

    for row_idx, row in out.iterrows():
        code = str(row.get("code", "")).strip().upper()
        if "ALL" not in CHROMATOPY_TARGET_CODES and code not in CHROMATOPY_TARGET_CODES:
            continue

        target_rt = pd.to_numeric(pd.Series([row.get("found_rt")]), errors="coerce").iloc[0]
        current_area = pd.to_numeric(pd.Series([row.get("area")]), errors="coerce").iloc[0]
        if not np.isfinite(target_rt) or not np.isfinite(current_area) or current_area <= 0:
            continue

        peak_idx = _nearest_peak_index(x, peak_indices, float(target_rt))
        if peak_idx is None:
            continue

        fit = _fit_one_peak(x, y_smooth, peak_indices, peak_idx)
        if fit is None:
            continue

        width = float(fit["integration_end_x"] - fit["integration_start_x"])
        area_ratio = float(fit["area"] / current_area)
        if width <= 0 or width > CHROMATOPY_MAX_WIDTH:
            continue
        if not (CHROMATOPY_MIN_AREA_RATIO <= area_ratio <= CHROMATOPY_MAX_AREA_RATIO):
            continue

        out.at[row_idx, "found_rt"] = float(fit["found_rt"])
        out.at[row_idx, "area"] = float(fit["area"])
        out.at[row_idx, "integration_start_x"] = float(fit["integration_start_x"])
        out.at[row_idx, "integration_end_x"] = float(fit["integration_end_x"])
        out.at[row_idx, "status"] = f"{row.get('status', '')}_chromatopy_{fit['model']}".strip("_")

    total_area = float(pd.to_numeric(out["area"], errors="coerce").fillna(0.0).sum())
    if total_area > 0:
        out["percent_area"] = 100.0 * pd.to_numeric(out["area"], errors="coerce") / total_area
    return out
