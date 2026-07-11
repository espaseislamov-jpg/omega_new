from __future__ import annotations

import argparse
import csv
import importlib.metadata
import importlib.util
import json
import math
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from omega_path_compat import configure_windows_path_compat

configure_windows_path_compat()

import numpy as np
import pandas as pd
from scipy.integrate import simpson
from scipy.signal import find_peaks, savgol_filter
from scipy.stats import median_abs_deviation


PROJECT_DIR = Path(__file__).resolve().parent
DEFAULT_REFERENCE_PATH = PROJECT_DIR / "reference_targets_reverted_c22fixed.json"
DEFAULT_CUTOFF_MINUTES = 4.0

OMEGA_CODES = {"C20:5", "C22:5", "C22:6"}
NEGATIVE_VALLEY_INTEGRATION_CODES = {
    "C18:3N6",
    "C18:2N6C",
    "C18:1N9C",
    "C18:3N3",
    "C18:0",
    "C20:4N6",
    "C20:5",
    "C20:3N8",
    "C22:6",
    "C22:5",
    "C22:4",
}
OMEGA_STRICT_WINDOWS = {
    "C20:5": (8.398, 8.438),
    "C22:6": (9.225, 9.268),
    "C22:5": (9.268, 9.300),
}
ENABLE_C22_DPA_SPLIT_RULE = True
C22_DPA_HIGH_RATIO_TRIGGER = 1.15
ENABLE_C20_SHOULDER_VALLEY_MODE = True
C20_SHOULDER_EPA_TO_C20_3_RANGE = (0.28, 0.85)
C20_SHOULDER_EPA_TO_C20_4_RANGE = (0.045, 0.13)
SAMPLE_NAME_RE = re.compile(r"\bO\d+_[A-Za-z0-9._-]+\b")

DEFAULT_TARGET_RTS = {
    "C16:1N7": 6.600,
    "C16:0": 6.700,
    "C18:3N6": 7.500,
    "C18:2N6C": 7.590,
    "C18:1N9C": 7.625,
    "C18:3N3": 7.653,
    "C18:0": 7.755,
    "C20:4N6": 8.384,
    "C20:5": 8.414,
    "C20:3N8": 8.469,
    "C22:6": 9.250,
    "C22:5": 9.283,
    "C22:4": 9.313,
    "C24:1N9": 10.302,
    "C24:0": 10.395,
}

ANCHOR_CODES = ("C16:1N7", "C16:0", "C24:1N9", "C24:0")
TARGET_WINDOWS = {
    "C16:1N7": (6.54, 6.66),
    "C16:0": (6.66, 6.75),
    "C18:3N6": (7.45, 7.55),
    "C18:2N6C": (7.55, 7.61),
    "C18:1N9C": (7.61, 7.645),
    "C18:3N3": (7.64, 7.70),
    "C18:0": (7.70, 7.82),
    "C20:4N6": (8.34, 8.402),
    "C20:5": (8.398, 8.440),
    "C20:3N8": (8.435, 8.515),
    "C22:6": (9.22, 9.268),
    "C22:5": (9.265, 9.298),
    "C22:4": (9.295, 9.34),
    "C24:1N9": (10.24, 10.35),
    "C24:0": (10.35, 10.43),
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


def load_chromatopy_fid_functions():
    global _CHROMATOPY_FUNCTIONS
    if _CHROMATOPY_FUNCTIONS is not None:
        return _CHROMATOPY_FUNCTIONS

    module_path = _chromatopy_fid_module_path()
    spec = importlib.util.spec_from_file_location("_omega_clean_chromatopy_fid_functions", module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load ChromatoPy FID functions from {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    _CHROMATOPY_FUNCTIONS = {
        "baseline": module.baseline,
        "calculate_boundaries": module.calculate_boundaries,
        "fit_gaussians": module.fit_gaussians,
        "smoother": module.smoother,
    }
    return _CHROMATOPY_FUNCTIONS


@dataclass(frozen=True)
class IntegrationConfig:
    cutoff_minutes: float = DEFAULT_CUTOFF_MINUTES
    smoothing_window: int = 7
    smoothing_polyorder: int = 3
    fit_mode: str = "single"
    gaussian_iterations: int = 500
    derivative_sensitivity: float = 0.001
    match_tolerance: float = 0.060
    peak_prominence_sigma: float = 1.8
    peak_height_sigma: float = 1.0
    min_peak_distance_minutes: float = 0.010
    max_fit_width_minutes: float = 0.180
    use_chromatopy_fit: bool = False


def robust_sigma(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return 0.0
    sigma = float(median_abs_deviation(arr, scale="normal", nan_policy="omit"))
    if sigma <= 0:
        sigma = float(np.std(arr))
    return sigma


def omega_sample_sort_key(name: str) -> tuple[int, int | str, str]:
    text = name or ""
    match = re.search(r"\bO(\d+)\b", text) or re.search(r"\bO(\d+)", text)
    if match:
        return (0, int(match.group(1)), text)
    return (1, text, text)


def load_reference_targets(reference_path: Path = DEFAULT_REFERENCE_PATH) -> pd.DataFrame:
    with Path(reference_path).open("r", encoding="utf-8") as handle:
        raw_targets = json.load(handle)

    rows = []
    for idx, item in enumerate(raw_targets, start=1):
        code = str(item.get("code", item.get("component", f"target_{idx}"))).strip()
        expected_rt = DEFAULT_TARGET_RTS.get(code, item.get("expected_rt", item.get("target_rt")))
        rows.append({
            "component": item.get("component", code),
            "code": code,
            "display_name": item.get("display_name", code),
            "order_index": int(item.get("order_index", idx)),
            "expected_rt": expected_rt,
            "rt_reliable": bool(item.get("rt_reliable", False)) or code in ANCHOR_CODES,
        })

    out = pd.DataFrame(rows)
    out["expected_rt"] = pd.to_numeric(out["expected_rt"], errors="coerce")
    return out.sort_values("order_index").reset_index(drop=True)


def finalize_chromatogram_dataframe(df: pd.DataFrame, cutoff_minutes: float) -> pd.DataFrame:
    out = df.copy()
    out["x"] = pd.to_numeric(out["x"], errors="coerce")
    out["y"] = pd.to_numeric(out["y"], errors="coerce")
    out = out.dropna(subset=["x", "y"]).reset_index(drop=True)
    if out.empty:
        raise ValueError("Chromatogram contains no numeric points.")

    unique_x = np.sort(out["x"].unique())
    if len(unique_x) < 2:
        raise ValueError("Not enough unique x values.")
    step = float(np.median(np.diff(unique_x)))

    corrected_x: list[float] = []
    x_values = out["x"].to_numpy(dtype=float)
    index = 0
    while index < len(out):
        current_x = x_values[index]
        end = index + 1
        while end < len(out) and x_values[end] == current_x:
            end += 1
        group_size = end - index
        offsets = np.linspace(0.0, step * (group_size - 1) / max(group_size, 1), group_size)
        corrected_x.extend((current_x + offsets).tolist())
        index = end

    out["x"] = corrected_x
    return out[out["x"] >= cutoff_minutes].reset_index(drop=True)


def extract_sample_name(text: str, fallback: str) -> str:
    match = SAMPLE_NAME_RE.search(text or "")
    return match.group(0) if match else fallback


def is_chromtab_file(file_path: Path) -> bool:
    try:
        with Path(file_path).open("r", encoding="utf-8", errors="ignore") as handle:
            first = next(handle, "").strip()
            next(handle, "")
            third = next(handle, "").strip()
    except OSError:
        return False
    return first.startswith('"Path","File","Date Acquired"') and third.startswith('"Signal: ')


def load_chromtab_batches(file_path: Path, cutoff_minutes: float) -> list[dict[str, Any]]:
    batches: list[dict[str, Any]] = []
    current_meta: dict[str, str] | None = None
    current_rows: list[tuple[float, float]] = []

    def flush_current() -> None:
        nonlocal current_meta, current_rows
        if current_meta is None or not current_rows:
            current_rows = []
            return
        raw_df = pd.DataFrame(current_rows, columns=["x", "y"])
        chrom_df = finalize_chromatogram_dataframe(raw_df, cutoff_minutes=cutoff_minutes)
        file_name = current_meta.get("file_name", "")
        sample_name = extract_sample_name(
            " ".join(filter(None, [current_meta.get("signal_name", ""), file_name])),
            fallback=Path(file_name).stem if file_name else f"batch_{len(batches) + 1}",
        )
        batches.append({
            "sample_name": sample_name,
            "file_name": file_name,
            "dataframe": chrom_df,
        })
        current_meta = None
        current_rows = []

    with Path(file_path).open("r", encoding="utf-8", errors="ignore") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('"Path","File","Date Acquired"'):
                flush_current()
                current_meta = {"file_name": "", "signal_name": ""}
                continue
            if current_meta is not None and not current_meta["file_name"] and line.startswith('"'):
                parsed = next(csv.reader([line]))
                current_meta["file_name"] = parsed[1] if len(parsed) > 1 else ""
                continue
            if current_meta is not None and line.startswith('"Signal: '):
                current_meta["signal_name"] = line.strip('"')
                continue
            if current_meta is None:
                continue
            parsed = next(csv.reader([line]), [])
            if len(parsed) < 2:
                continue
            try:
                current_rows.append((float(parsed[0]), float(parsed[1])))
            except ValueError:
                continue

    flush_current()
    if not batches:
        raise ValueError(f"No batches parsed from {file_path}")
    return sorted(batches, key=lambda item: omega_sample_sort_key(item["sample_name"]))


def load_batches(file_path: Path, cutoff_minutes: float = DEFAULT_CUTOFF_MINUTES) -> list[dict[str, Any]]:
    file_path = Path(file_path)
    if is_chromtab_file(file_path):
        return load_chromtab_batches(file_path, cutoff_minutes=cutoff_minutes)

    try:
        with file_path.open("r", encoding="utf-8", errors="ignore") as handle:
            header = " ".join(next(handle, "") for _ in range(3))
    except OSError:
        header = ""
    raw_df = pd.read_csv(file_path, skiprows=3, header=None, names=["x", "y"])
    return [{
        "sample_name": extract_sample_name(header, fallback=file_path.stem),
        "file_name": file_path.name,
        "dataframe": finalize_chromatogram_dataframe(raw_df, cutoff_minutes=cutoff_minutes),
    }]


def preprocess_signal(df: pd.DataFrame, config: IntegrationConfig) -> pd.DataFrame:
    chromatopy_functions = load_chromatopy_fid_functions()
    chromatopy_smoother = chromatopy_functions["smoother"]
    chromatopy_baseline = chromatopy_functions["baseline"]
    x = df["x"].to_numpy(dtype=float)
    y = df["y"].to_numpy(dtype=float)
    y_smooth_for_base = chromatopy_smoother(y, config.smoothing_window, config.smoothing_polyorder)
    base, min_peak_amp = chromatopy_baseline(x, np.asarray(y_smooth_for_base, dtype=float), deg=5, max_it=1000, tol=1e-4)
    y_corrected = np.clip(y - np.asarray(base, dtype=float), 0.0, None)
    y_smooth = savgol_filter(
        y_corrected,
        window_length=min(config.smoothing_window if config.smoothing_window % 2 else config.smoothing_window + 1, len(y_corrected) - (1 - len(y_corrected) % 2)),
        polyorder=min(config.smoothing_polyorder, max(1, len(y_corrected) - 2)),
        mode="interp",
    ) if len(y_corrected) > 8 else y_corrected
    out = df.copy()
    out["baseline"] = base
    out["y_corrected"] = y_corrected
    out["y_smooth"] = np.clip(y_smooth, 0.0, None)
    out.attrs["chromatopy_min_peak_amp"] = float(min_peak_amp)
    return out


def detect_peaks(processed: pd.DataFrame, config: IntegrationConfig) -> tuple[np.ndarray, dict]:
    x = processed["x"].to_numpy(dtype=float)
    y = processed["y_smooth"].to_numpy(dtype=float)
    dx = float(np.median(np.diff(x))) if len(x) > 1 else 0.001
    noise = max(robust_sigma(processed["y_corrected"].to_numpy(dtype=float)), 1.0)
    min_distance = max(1, int(round(config.min_peak_distance_minutes / max(dx, 1e-9))))
    height = max(noise * config.peak_height_sigma, float(np.quantile(y, 0.60)))
    prominence = max(noise * config.peak_prominence_sigma, float(np.quantile(y, 0.75)) * 0.03, 5.0)
    peaks, props = find_peaks(y, height=height, prominence=prominence, distance=min_distance)
    return peaks.astype(int), props


def detect_negative_valleys(processed: pd.DataFrame, peaks: np.ndarray, config: IntegrationConfig) -> np.ndarray:
    """Mark valleys as negative peaks so boundaries are explicit objects."""
    x = processed["x"].to_numpy(dtype=float)
    y = processed["y_smooth"].to_numpy(dtype=float)
    if len(x) < 5:
        return np.array([], dtype=int)

    dx = float(np.median(np.diff(x))) if len(x) > 1 else 0.001
    min_distance = max(1, int(round(config.min_peak_distance_minutes / max(dx, 1e-9))))
    noise = max(robust_sigma(processed["y_corrected"].to_numpy(dtype=float)), 1.0)
    prominence = max(noise * 0.35, float(np.nanmax(y)) * 0.003, 1.0)

    valley_set: set[int] = set()
    negative_peaks, _ = find_peaks(-y, prominence=prominence, distance=min_distance)
    valley_set.update(int(idx) for idx in negative_peaks)

    ordered_peaks = np.asarray(sorted(set(int(idx) for idx in peaks)), dtype=int)
    for left_peak, right_peak in zip(ordered_peaks[:-1], ordered_peaks[1:]):
        if right_peak - left_peak < 3:
            continue
        valley_idx = int(np.argmin(y[left_peak:right_peak + 1]) + left_peak)
        valley_set.add(valley_idx)

    valleys = np.asarray(sorted(valley_set), dtype=int)
    valleys = valleys[(valleys > 0) & (valleys < len(x) - 1)]
    return valleys.astype(int)


def negative_valley_area(
    x: np.ndarray,
    y: np.ndarray,
    peak_idx: int,
    valley_indices: np.ndarray,
    left_limit_x: float,
    right_limit_x: float,
) -> dict[str, float] | None:
    if valley_indices is None or len(valley_indices) == 0:
        return None
    peak_idx = int(peak_idx)
    left_limit_idx = int(np.searchsorted(x, left_limit_x, side="left"))
    right_limit_idx = int(np.searchsorted(x, right_limit_x, side="right") - 1)
    left_limit_idx = max(0, min(left_limit_idx, peak_idx))
    right_limit_idx = min(len(x) - 1, max(right_limit_idx, peak_idx))

    valleys = np.asarray(valley_indices, dtype=int)
    left_candidates = valleys[(valleys < peak_idx) & (valleys >= left_limit_idx)]
    right_candidates = valleys[(valleys > peak_idx) & (valleys <= right_limit_idx)]
    if len(left_candidates) == 0 or len(right_candidates) == 0:
        return None

    left_idx = int(left_candidates[-1])
    right_idx = int(right_candidates[0])
    if right_idx <= left_idx or left_idx >= peak_idx or right_idx <= peak_idx:
        return None
    if float(x[right_idx] - x[left_idx]) < 0.006:
        return None

    y_seg = np.clip(y[left_idx:right_idx + 1], 0.0, None)
    area = float(np.trapezoid(y_seg, x[left_idx:right_idx + 1]))
    if not np.isfinite(area) or area <= 0:
        return None
    return {
        "start_x": float(x[left_idx]),
        "end_x": float(x[right_idx]),
        "area": area,
    }


def estimate_rt_shift(reference: pd.DataFrame, x: np.ndarray, peaks: np.ndarray) -> float:
    shifts = []
    if peaks.size == 0:
        return 0.0
    for code in ANCHOR_CODES:
        row = reference[reference["code"] == code]
        if row.empty:
            continue
        expected = float(row.iloc[0]["expected_rt"])
        distances = np.abs(x[peaks] - expected)
        best_pos = int(np.argmin(distances))
        if float(distances[best_pos]) <= 0.060:
            shifts.append(float(x[peaks[best_pos]] - expected))
    return float(np.median(shifts)) if shifts else 0.0


def target_peak_index(
    x: np.ndarray,
    y: np.ndarray,
    peaks: np.ndarray,
    target_rt: float,
    tolerance: float,
    used_peaks: set[int],
    left_limit_x: float,
    right_limit_x: float,
) -> int | None:
    if not np.isfinite(target_rt):
        return None

    available_peaks = np.array([
        int(peak)
        for peak in peaks
        if int(peak) not in used_peaks and left_limit_x <= float(x[int(peak)]) <= right_limit_x
    ], dtype=int)
    if available_peaks.size > 0:
        distances = np.abs(x[available_peaks] - target_rt)
        best_pos = int(np.argmin(distances))
        if float(distances[best_pos]) <= tolerance:
            return int(available_peaks[best_pos])

    window_mask = (x >= left_limit_x) & (x <= right_limit_x)
    candidate_indices = np.flatnonzero(window_mask)
    candidate_indices = np.array([int(idx) for idx in candidate_indices if int(idx) not in used_peaks], dtype=int)
    if candidate_indices.size == 0:
        return None
    return int(candidate_indices[int(np.argmax(y[candidate_indices]))])


def strict_omega_peak_index(
    x: np.ndarray,
    y: np.ndarray,
    peaks: np.ndarray,
    code: str,
    target_rt: float,
    rt_shift: float,
    used_peaks: set[int],
) -> tuple[int | None, str]:
    if code not in OMEGA_STRICT_WINDOWS:
        return None, "not_omega_strict"

    left, right = OMEGA_STRICT_WINDOWS[code]
    left += rt_shift
    right += rt_shift
    candidates = [
        int(peak)
        for peak in peaks
        if int(peak) not in used_peaks and left <= float(x[int(peak)]) <= right
    ]

    if candidates:
        best = int(max(candidates, key=lambda idx: (float(y[idx]), -abs(float(x[idx]) - target_rt))))
        return best, "strict_peak"

    window_mask = (x >= left) & (x <= right)
    indices = np.flatnonzero(window_mask)
    indices = np.array([int(idx) for idx in indices if int(idx) not in used_peaks], dtype=int)
    if indices.size == 0:
        return None, "strict_empty_window"
    best = int(indices[int(np.argmax(y[indices]))])
    if y[best] <= max(robust_sigma(y[indices]) * 0.5, 1.0):
        return None, "strict_no_signal"
    return best, "strict_local_max"


def bounded_valley_area(
    x: np.ndarray,
    y: np.ndarray,
    peak_idx: int,
    left_limit_x: float,
    right_limit_x: float,
) -> dict[str, float]:
    left_limit_idx = int(np.searchsorted(x, left_limit_x, side="left"))
    right_limit_idx = int(np.searchsorted(x, right_limit_x, side="right") - 1)
    left_limit_idx = max(0, min(left_limit_idx, int(peak_idx)))
    right_limit_idx = min(len(x) - 1, max(right_limit_idx, int(peak_idx)))

    left_idx = int(np.argmin(y[left_limit_idx:peak_idx + 1]) + left_limit_idx)
    right_idx = int(np.argmin(y[peak_idx:right_limit_idx + 1]) + peak_idx)
    if right_idx <= left_idx:
        left_idx = max(0, peak_idx - 4)
        right_idx = min(len(x) - 1, peak_idx + 4)
    area = float(np.trapezoid(y[left_idx:right_idx + 1], x[left_idx:right_idx + 1]))
    return {
        "start_x": float(x[left_idx]),
        "end_x": float(x[right_idx]),
        "area": max(area, 0.0),
    }


def chromatopy_numeric_area(
    x: np.ndarray,
    y: np.ndarray,
    peak_idx: int,
    left_limit_x: float,
    right_limit_x: float,
    config: IntegrationConfig,
) -> dict[str, float]:
    functions = load_chromatopy_fid_functions()
    calculate_boundaries_func = functions["calculate_boundaries"]

    left_limit_idx = int(np.searchsorted(x, left_limit_x, side="left"))
    right_limit_idx = int(np.searchsorted(x, right_limit_x, side="right") - 1)
    left_limit_idx = max(0, min(left_limit_idx, int(peak_idx)))
    right_limit_idx = min(len(x) - 1, max(right_limit_idx, int(peak_idx)))
    if right_limit_idx <= left_limit_idx:
        return bounded_valley_area(x, y, peak_idx, left_limit_x, right_limit_x)

    x_local = x[left_limit_idx:right_limit_idx + 1]
    y_local = y[left_limit_idx:right_limit_idx + 1]
    local_peak_idx = int(peak_idx) - left_limit_idx
    try:
        local_left, local_right = calculate_boundaries_func(
            pd.Series(x_local),
            pd.Series(y_local),
            local_peak_idx,
            [config.smoothing_window, config.smoothing_polyorder],
            config.derivative_sensitivity,
        )
        start_idx = left_limit_idx + int(local_left)
        end_idx = left_limit_idx + int(local_right)
    except Exception:
        return bounded_valley_area(x, y, peak_idx, left_limit_x, right_limit_x)

    start_idx = max(left_limit_idx, min(start_idx, int(peak_idx)))
    end_idx = min(right_limit_idx, max(end_idx, int(peak_idx)))
    if end_idx <= start_idx or start_idx >= int(peak_idx) or end_idx <= int(peak_idx):
        return bounded_valley_area(x, y, peak_idx, left_limit_x, right_limit_x)
    if x[end_idx] - x[start_idx] < 0.010:
        return bounded_valley_area(x, y, peak_idx, left_limit_x, right_limit_x)

    # ChromatoPy's derivative crossing can occur on a shoulder while the peak is
    # still descending.  Treat it as an initial boundary and complete both sides
    # to the local signal minima inside the target RT corridor.  This preserves
    # separation from neighbouring targets while preventing visibly truncated
    # integrations.
    valley = bounded_valley_area(x, y, peak_idx, left_limit_x, right_limit_x)
    valley_start_idx = int(np.searchsorted(x, float(valley["start_x"]), side="left"))
    valley_end_idx = int(np.searchsorted(x, float(valley["end_x"]), side="right") - 1)
    start_idx = max(left_limit_idx, min(start_idx, valley_start_idx))
    end_idx = min(right_limit_idx, max(end_idx, valley_end_idx))

    x_seg = x[start_idx:end_idx + 1]
    y_seg = np.clip(y[start_idx:end_idx + 1], 0.0, None)
    if len(x_seg) >= 3:
        area = float(simpson(y_seg, x=x_seg))
    else:
        area = float(np.trapezoid(y_seg, x_seg))
    if not np.isfinite(area) or area <= 0:
        return bounded_valley_area(x, y, peak_idx, left_limit_x, right_limit_x)
    return {
        "start_x": float(x[start_idx]),
        "end_x": float(x[end_idx]),
        "area": area,
    }


def chromatopy_fit_area(
    x: np.ndarray,
    y: np.ndarray,
    peaks: np.ndarray,
    peak_idx: int,
    config: IntegrationConfig,
) -> dict[str, float] | None:
    fit_gaussians = load_chromatopy_fid_functions()["fit_gaussians"]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fit_x, fit_y, _area_smooth, area_ensemble, model = fit_gaussians(
                pd.Series(x),
                pd.Series(y),
                int(peak_idx),
                np.array([int(peak_idx)], dtype=int),
                [config.smoothing_window, config.smoothing_polyorder],
                config.derivative_sensitivity,
                gi=config.gaussian_iterations,
                mode=config.fit_mode,
            )
    except Exception:
        return None

    fit_x = np.asarray(fit_x, dtype=float)
    fit_y = np.asarray(fit_y, dtype=float)
    area_ensemble = np.asarray(area_ensemble, dtype=float)
    area_ensemble = area_ensemble[np.isfinite(area_ensemble)]
    if fit_x.size < 2 or fit_y.size < 2 or area_ensemble.size == 0:
        return None
    width = float(np.nanmax(fit_x) - np.nanmin(fit_x))
    if width <= 0 or width > config.max_fit_width_minutes:
        return None
    apex_pos = int(np.nanargmax(fit_y))
    return {
        "found_rt": float(fit_x[apex_pos]),
        "start_x": float(np.nanmin(fit_x)),
        "end_x": float(np.nanmax(fit_x)),
        "area": float(np.nanmedian(area_ensemble)),
        "model": str(model.get("name", config.fit_mode)) if isinstance(model, dict) else config.fit_mode,
    }


def build_target_limits(reference: pd.DataFrame, rt_shift: float) -> dict[str, tuple[float, float]]:
    targets = reference[["code", "expected_rt"]].dropna().copy()
    targets["expected_rt"] = pd.to_numeric(targets["expected_rt"], errors="coerce")
    targets = targets.dropna(subset=["expected_rt"]).sort_values("expected_rt").reset_index(drop=True)
    limits: dict[str, tuple[float, float]] = {}
    for pos, row in targets.iterrows():
        code = str(row["code"])
        rt = float(row["expected_rt"])
        left = rt - 0.070
        right = rt + 0.070
        if pos > 0:
            left = 0.5 * (float(targets.iloc[pos - 1]["expected_rt"]) + rt)
        if pos < len(targets) - 1:
            right = 0.5 * (rt + float(targets.iloc[pos + 1]["expected_rt"]))
        if code in TARGET_WINDOWS:
            window_left, window_right = TARGET_WINDOWS[code]
            left = max(left, window_left)
            right = min(right, window_right)
        if right <= left:
            left, right = rt - 0.020, rt + 0.020
        limits[code] = (left + rt_shift, right + rt_shift)
    return limits


def add_omega_diagnostics(matched: pd.DataFrame) -> pd.DataFrame:
    out = matched.copy()
    out["omega_diagnostic"] = ""

    area_by_code = {
        str(row["code"]): float(row["area"])
        for _, row in out.dropna(subset=["area"]).iterrows()
        if np.isfinite(float(row["area"]))
    }
    epa = area_by_code.get("C20:5", np.nan)
    c20_3 = area_by_code.get("C20:3N8", np.nan)
    dha = area_by_code.get("C22:6", np.nan)
    dpa = area_by_code.get("C22:5", np.nan)
    c22_4 = area_by_code.get("C22:4", np.nan)

    diagnostics: dict[str, list[str]] = {code: [] for code in ["C20:5", "C22:6", "C22:5", "C22:4"]}
    if np.isfinite(epa) and np.isfinite(c20_3) and c20_3 > 0:
        ratio = epa / c20_3
        if ratio < 0.35:
            diagnostics["C20:5"].append(f"epa_low_vs_c20_3={ratio:.2f}")
        if ratio > 2.80:
            diagnostics["C20:5"].append(f"epa_high_vs_c20_3={ratio:.2f}")

    if np.isfinite(dpa) and np.isfinite(c22_4) and c22_4 > 0:
        ratio = dpa / c22_4
        if ratio > 1.15:
            diagnostics["C22:5"].append(f"dpa_high_vs_c22_4={ratio:.2f}")
        if ratio < 0.20:
            diagnostics["C22:5"].append(f"dpa_low_vs_c22_4={ratio:.2f}")

    if np.isfinite(dha) and np.isfinite(dpa) and dpa > 0:
        ratio = dha / dpa
        if ratio < 1.20:
            diagnostics["C22:6"].append(f"dha_low_vs_dpa={ratio:.2f}")
            diagnostics["C22:5"].append(f"dpa_high_vs_dha={1.0 / max(ratio, 1e-9):.2f}")

    for code, notes in diagnostics.items():
        if not notes:
            continue
        mask = out["code"] == code
        out.loc[mask, "omega_diagnostic"] = ";".join(notes)
        out.loc[mask, "status"] = out.loc[mask, "status"].astype(str) + "_omega_check"
    return out


def _trapezoid_cumulative(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    if len(x) <= 1:
        return np.zeros(len(x), dtype=float)
    increments = 0.5 * (y[:-1] + y[1:]) * np.diff(x)
    return np.concatenate([[0.0], np.cumsum(increments)])


def _area_between_x(processed: pd.DataFrame, start_x: float, end_x: float) -> float:
    if not (np.isfinite(start_x) and np.isfinite(end_x)) or end_x <= start_x:
        return np.nan
    x = processed["x"].to_numpy(dtype=float)
    y = np.clip(processed["y_corrected"].to_numpy(dtype=float), 0.0, None)
    start_idx = int(np.searchsorted(x, float(start_x), side="left"))
    end_idx = int(np.searchsorted(x, float(end_x), side="right") - 1)
    start_idx = max(0, min(start_idx, len(x) - 1))
    end_idx = max(0, min(end_idx, len(x) - 1))
    if end_idx <= start_idx:
        return np.nan
    return float(np.trapezoid(y[start_idx:end_idx + 1], x[start_idx:end_idx + 1]))


def _set_row_interval(out: pd.DataFrame, row_idx: int, processed: pd.DataFrame, start_x: float, end_x: float, suffix: str) -> None:
    area = _area_between_x(processed, start_x, end_x)
    if not np.isfinite(area) or area <= 0:
        return
    out.at[row_idx, "integration_start_x"] = float(start_x)
    out.at[row_idx, "integration_end_x"] = float(end_x)
    out.at[row_idx, "area"] = float(area)
    out.at[row_idx, "status"] = f"{out.at[row_idx, 'status']}_{suffix}"


def _local_valley_x(processed: pd.DataFrame, left_rt: float, right_rt: float) -> float:
    if not (np.isfinite(left_rt) and np.isfinite(right_rt)) or right_rt <= left_rt:
        return np.nan
    x = processed["x"].to_numpy(dtype=float)
    y = processed["y_smooth"].to_numpy(dtype=float)
    mask = (x >= float(left_rt)) & (x <= float(right_rt))
    indices = np.flatnonzero(mask)
    if indices.size < 3:
        return np.nan
    local_idx = int(indices[int(np.argmin(y[indices]))])
    return float(x[local_idx])


def _last_preceding_peak_x(
    processed: pd.DataFrame,
    found_rt: float,
    left_limit: float,
    min_prominence_fraction: float = 0.08,
    apex_exclusion_minutes: float = 0.012,
) -> float:
    if not (np.isfinite(found_rt) and np.isfinite(left_limit)) or found_rt <= left_limit:
        return np.nan
    x = processed["x"].to_numpy(dtype=float)
    y = processed["y_smooth"].to_numpy(dtype=float)
    mask = (x >= float(left_limit)) & (x < float(found_rt) - float(apex_exclusion_minutes))
    indices = np.flatnonzero(mask)
    if indices.size < 5:
        return np.nan
    y_local = y[indices]
    peak_height = float(y[int(np.argmin(np.abs(x - found_rt)))])
    prominence = max(peak_height * min_prominence_fraction, robust_sigma(y_local) * 0.45, 1.0)
    peaks, _ = find_peaks(y_local, prominence=prominence, distance=3)
    if len(peaks) == 0:
        max_pos = int(np.argmax(y_local))
        if float(y_local[max_pos]) < max(peak_height * min_prominence_fraction, 1.0):
            return np.nan
        return float(x[indices[max_pos]])
    return float(x[indices[int(peaks[-1])]])


def apply_interpeak_boundary_guards(processed: pd.DataFrame, matched: pd.DataFrame, rt_shift: float) -> pd.DataFrame:
    out = matched.copy()
    if processed is None or processed.empty or out.empty:
        return out

    # C22:6 must not absorb a separate pre-DHA peak; split at the valley after that pre-peak.
    dha_rows = out.index[out["code"] == "C22:6"].tolist()
    if dha_rows:
        dha_idx = dha_rows[0]
        dha_rt = pd.to_numeric(pd.Series([out.at[dha_idx, "found_rt"]]), errors="coerce").iloc[0]
        dha_start = pd.to_numeric(pd.Series([out.at[dha_idx, "integration_start_x"]]), errors="coerce").iloc[0]
        dha_end = pd.to_numeric(pd.Series([out.at[dha_idx, "integration_end_x"]]), errors="coerce").iloc[0]
        pre_peak_x = _last_preceding_peak_x(processed, float(dha_rt), 9.18 + float(rt_shift))
        if np.isfinite(pre_peak_x):
            valley_x = _local_valley_x(processed, pre_peak_x, float(dha_rt))
            x_values = processed["x"].to_numpy(dtype=float)
            y_values = processed["y_smooth"].to_numpy(dtype=float)
            pre_y = float(y_values[int(np.argmin(np.abs(x_values - pre_peak_x)))])
            valley_y = float(y_values[int(np.argmin(np.abs(x_values - valley_x)))]) if np.isfinite(valley_x) else np.nan
            dha_y = float(y_values[int(np.argmin(np.abs(x_values - float(dha_rt))))])
            separate_prepeak = (
                np.isfinite(valley_x)
                and np.isfinite(pre_y)
                and np.isfinite(valley_y)
                and np.isfinite(dha_y)
                and pre_y >= dha_y * 0.18
                and valley_y <= pre_y * 0.72
            )
            if separate_prepeak and np.isfinite(dha_start) and np.isfinite(dha_end) and dha_start < valley_x < dha_end:
                _set_row_interval(out, dha_idx, processed, valley_x, float(dha_end), "prepeak_guard")

    # Do not cap an EPA interval merely because its apex is near the edge of the
    # expected RT window.  RT windows identify peaks; they are not integration
    # boundaries.  The previous ``apex + 0.014`` cap cut asymmetric peaks through
    # a still-descending tail (notably the RT 8.3893 field case).  The numeric
    # ChromatoPy/valley boundary selected above must remain authoritative.

    return out


def _target_dpa_to_c22_4_ratio(dpa_to_c22_4: float, dha_to_dpa: float) -> float:
    if dha_to_dpa < 2.20:
        return 0.55
    if dha_to_dpa < 3.00:
        return 0.75
    if dpa_to_c22_4 > 1.35:
        return 0.85
    return 1.00


def apply_c22_dpa_split_rule(processed: pd.DataFrame, matched: pd.DataFrame) -> pd.DataFrame:
    out = matched.copy()
    if not ENABLE_C22_DPA_SPLIT_RULE or processed is None or processed.empty or out.empty:
        return out

    rows = {}
    for code in ["C22:6", "C22:5", "C22:4"]:
        row = out[out["code"] == code]
        if row.empty:
            return out
        rows[code] = row.iloc[0]

    dpa_area = float(pd.to_numeric(pd.Series([rows["C22:5"].get("area")]), errors="coerce").iloc[0])
    c22_4_area = float(pd.to_numeric(pd.Series([rows["C22:4"].get("area")]), errors="coerce").iloc[0])
    dha_area = float(pd.to_numeric(pd.Series([rows["C22:6"].get("area")]), errors="coerce").iloc[0])
    if not (np.isfinite(dpa_area) and np.isfinite(c22_4_area) and np.isfinite(dha_area)):
        return out
    if dpa_area <= 0 or c22_4_area <= 0:
        return out

    dpa_to_c22_4 = dpa_area / c22_4_area
    dha_to_dpa = dha_area / dpa_area
    if dpa_to_c22_4 <= C22_DPA_HIGH_RATIO_TRIGGER:
        return out

    target_ratio = _target_dpa_to_c22_4_ratio(dpa_to_c22_4, dha_to_dpa)
    if target_ratio >= dpa_to_c22_4:
        return out

    dpa_idx = out.index[out["code"] == "C22:5"][0]
    c22_4_idx = out.index[out["code"] == "C22:4"][0]
    dpa_start = float(rows["C22:5"].get("integration_start_x"))
    c22_4_end = float(rows["C22:4"].get("integration_end_x"))
    dpa_rt = float(rows["C22:5"].get("found_rt"))
    c22_4_rt = float(rows["C22:4"].get("found_rt"))
    if not all(np.isfinite(value) for value in [dpa_start, c22_4_end, dpa_rt, c22_4_rt]):
        return out
    if not (dpa_start < dpa_rt < c22_4_rt < c22_4_end):
        return out

    x = processed["x"].to_numpy(dtype=float)
    y = np.clip(processed["y_corrected"].to_numpy(dtype=float), 0.0, None)
    mask = (x >= dpa_start) & (x <= c22_4_end)
    if int(np.count_nonzero(mask)) < 5:
        return out
    x_seg = x[mask]
    y_seg = y[mask]
    cumulative = _trapezoid_cumulative(x_seg, y_seg)
    combined_area = float(cumulative[-1])
    if not np.isfinite(combined_area) or combined_area <= 0:
        return out

    desired_dpa_area = combined_area * target_ratio / (1.0 + target_ratio)
    if desired_dpa_area >= dpa_area:
        return out

    split_pos = int(np.searchsorted(cumulative, desired_dpa_area, side="left"))
    split_pos = max(1, min(split_pos, len(x_seg) - 2))
    valley_x = _local_valley_x(processed, dpa_rt, c22_4_rt)
    if np.isfinite(valley_x):
        valley_pos = int(np.searchsorted(x_seg, valley_x, side="left"))
        valley_pos = max(1, min(valley_pos, len(x_seg) - 2))
        split_pos = max(split_pos, valley_pos)
    split_x = float(x_seg[split_pos])
    if not (dpa_rt < split_x < c22_4_rt):
        return out

    new_dpa_area = float(np.trapezoid(y_seg[:split_pos + 1], x_seg[:split_pos + 1]))
    new_c22_4_area = float(np.trapezoid(y_seg[split_pos:], x_seg[split_pos:]))
    if not (np.isfinite(new_dpa_area) and np.isfinite(new_c22_4_area)):
        return out
    if new_dpa_area <= 0 or new_c22_4_area <= 0:
        return out

    out.at[dpa_idx, "area"] = new_dpa_area
    out.at[c22_4_idx, "area"] = new_c22_4_area
    out.at[dpa_idx, "integration_end_x"] = split_x
    out.at[c22_4_idx, "integration_start_x"] = split_x
    suffix = f"_dpa_split_cap_{target_ratio:.2f}"
    if np.isfinite(valley_x) and abs(split_x - float(valley_x)) <= 0.0015:
        suffix += "_valley_limited"
    out.at[dpa_idx, "status"] = str(out.at[dpa_idx, "status"]) + suffix
    out.at[c22_4_idx, "status"] = str(out.at[c22_4_idx, "status"]) + suffix
    return out


def should_use_c20_shoulder_valley_mode(matched: pd.DataFrame) -> bool:
    if not ENABLE_C20_SHOULDER_VALLEY_MODE or matched is None or matched.empty:
        return False

    area_by_code = {
        str(row["code"]): float(row["area"])
        for _, row in matched.dropna(subset=["area"]).iterrows()
        if np.isfinite(float(row["area"]))
    }
    epa = area_by_code.get("C20:5", np.nan)
    c20_3 = area_by_code.get("C20:3N8", np.nan)
    c20_4 = area_by_code.get("C20:4N6", np.nan)
    if not all(np.isfinite(value) and value > 0 for value in [epa, c20_3, c20_4]):
        return False

    epa_to_c20_3 = epa / c20_3
    epa_to_c20_4 = epa / c20_4
    min_203, max_203 = C20_SHOULDER_EPA_TO_C20_3_RANGE
    min_204, max_204 = C20_SHOULDER_EPA_TO_C20_4_RANGE
    return (
        min_203 <= epa_to_c20_3 <= max_203
        and min_204 <= epa_to_c20_4 <= max_204
    )


def _mark_c20_shoulder_valley_mode(matched: pd.DataFrame) -> pd.DataFrame:
    out = matched.copy()
    mask = out["code"].isin(["C20:4N6", "C20:5", "C20:3N8"])
    out.loc[mask, "status"] = out.loc[mask, "status"].astype(str) + "_c20_shoulder_valley_mode"
    return out


def _integrate_batch_core(
    dataframe: pd.DataFrame,
    reference: pd.DataFrame,
    config: IntegrationConfig,
    boundary_mode: str = "chromatopy",
) -> dict[str, Any]:
    processed = preprocess_signal(dataframe, config)
    x = processed["x"].to_numpy(dtype=float)
    y = processed["y_corrected"].to_numpy(dtype=float)
    y_smooth = processed["y_smooth"].to_numpy(dtype=float)
    peaks, props = detect_peaks(processed, config)
    rt_shift = estimate_rt_shift(reference, x, peaks)
    target_limits = build_target_limits(reference, rt_shift)

    rows = []
    used_peaks: set[int] = set()
    for target in reference.itertuples(index=False):
        code = str(target.code)
        expected_rt = float(target.expected_rt) if np.isfinite(target.expected_rt) else DEFAULT_TARGET_RTS.get(code, np.nan)
        target_rt = expected_rt + rt_shift if np.isfinite(expected_rt) else np.nan
        left_limit, right_limit = target_limits.get(code, (target_rt - 0.050, target_rt + 0.050))
        omega_match_status = ""
        if code in OMEGA_STRICT_WINDOWS:
            peak_idx, omega_match_status = strict_omega_peak_index(
                x=x,
                y=y_smooth,
                peaks=peaks,
                code=code,
                target_rt=target_rt,
                rt_shift=rt_shift,
                used_peaks=used_peaks,
            )
        else:
            peak_idx = target_peak_index(
                x,
                y_smooth,
                peaks,
                target_rt,
                config.match_tolerance,
                used_peaks,
                left_limit,
                right_limit,
            )
        if peak_idx is None or peak_idx in used_peaks:
            rows.append({
                "code": code,
                "display_name": target.display_name,
                "expected_rt": expected_rt,
                "target_rt": target_rt,
                "found_rt": np.nan,
                "area": np.nan,
                "integration_start_x": np.nan,
                "integration_end_x": np.nan,
                "status": omega_match_status or "not_found",
            })
            continue
        used_peaks.add(peak_idx)

        fit = (
            chromatopy_fit_area(x, y_smooth, peaks, peak_idx, config)
            if config.use_chromatopy_fit and boundary_mode == "chromatopy"
            else None
        )
        if fit is None:
            if boundary_mode == "valley_all":
                fallback = bounded_valley_area(x, np.clip(y, 0.0, None), peak_idx, left_limit, right_limit)
                status = "matched_valley_numeric"
            else:
                fallback = chromatopy_numeric_area(x, y, peak_idx, left_limit, right_limit, config)
                status = "matched_chromatopy_numeric"
            found_rt = float(x[peak_idx])
            area = fallback["area"]
            start_x = fallback["start_x"]
            end_x = fallback["end_x"]
            if omega_match_status:
                status = f"{status}_{omega_match_status}"
        else:
            found_rt = fit["found_rt"]
            area = fit["area"]
            start_x = fit["start_x"]
            end_x = fit["end_x"]
            status = f"matched_chromatopy_{fit['model']}"

        rows.append({
            "code": code,
            "display_name": target.display_name,
            "expected_rt": expected_rt,
            "target_rt": target_rt,
            "found_rt": found_rt,
            "area": area,
            "integration_start_x": start_x,
            "integration_end_x": end_x,
            "status": status,
        })

    matched = pd.DataFrame(rows)
    matched = apply_interpeak_boundary_guards(processed, matched, rt_shift)
    matched = apply_c22_dpa_split_rule(processed, matched)
    total_area = float(pd.to_numeric(matched["area"], errors="coerce").fillna(0.0).sum())
    matched["percent_area"] = (
        100.0 * pd.to_numeric(matched["area"], errors="coerce") / total_area
        if total_area > 0 else np.nan
    )
    matched = add_omega_diagnostics(matched)
    omega_area = float(matched[matched["code"].isin(OMEGA_CODES)]["area"].fillna(0.0).sum())
    omega_value = 100.0 * omega_area / total_area if total_area > 0 else np.nan
    return {
        "processed_df": processed,
        "matched_targets_df": matched,
        "peaks": peaks,
        "peak_properties": props,
        "rt_shift": rt_shift,
        "omega3_trio": omega_value,
        "total_area": total_area,
        "boundary_mode": boundary_mode,
    }


def integrate_batch(dataframe: pd.DataFrame, reference: pd.DataFrame, config: IntegrationConfig) -> dict[str, Any]:
    result = _integrate_batch_core(dataframe, reference, config, boundary_mode="chromatopy")
    if not should_use_c20_shoulder_valley_mode(result["matched_targets_df"]):
        return result

    valley_result = _integrate_batch_core(dataframe, reference, config, boundary_mode="valley_all")
    valley_result["matched_targets_df"] = _mark_c20_shoulder_valley_mode(valley_result["matched_targets_df"])
    valley_result["boundary_mode"] = "c20_shoulder_valley"
    return valley_result


def process_file(
    file_path: Path,
    reference_path: Path = DEFAULT_REFERENCE_PATH,
    config: IntegrationConfig | None = None,
    sample_filter: str | None = None,
) -> list[dict[str, Any]]:
    config = config or IntegrationConfig()
    reference = load_reference_targets(reference_path)
    batches = load_batches(file_path, cutoff_minutes=config.cutoff_minutes)
    if sample_filter:
        needle = sample_filter.strip().lower()
        exact_match = re.fullmatch(r"o\d+", needle)
        if exact_match:
            pattern = re.compile(rf"\b{re.escape(needle)}(?:\b|_)", re.IGNORECASE)
            batches = [batch for batch in batches if pattern.search(batch["sample_name"])]
        else:
            batches = [batch for batch in batches if needle in batch["sample_name"].lower()]
    results = []
    for batch in batches:
        result = integrate_batch(batch["dataframe"], reference, config)
        results.append({**batch, **result})
    return results


def main() -> int:
    parser = argparse.ArgumentParser(description="Clean ChromatoPy-based Omega integration experiment.")
    parser.add_argument("csv_path", type=Path)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE_PATH)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--sample", default=None, help="Substring filter for sample name, e.g. O11.")
    parser.add_argument("--fit", action="store_true", help="Enable ChromatoPy Gaussian fit. Slow; use only for targeted experiments.")
    parser.add_argument("--no-fit", action="store_true", help="Compatibility flag; bounded local-valley mode is already the default.")
    parser.add_argument("--fit-mode", choices=["single", "multi", "both"], default="single")
    parser.add_argument("--gaus-iterations", type=int, default=500)
    args = parser.parse_args()

    config = IntegrationConfig(
        fit_mode=args.fit_mode,
        gaussian_iterations=args.gaus_iterations,
        use_chromatopy_fit=bool(args.fit and not args.no_fit),
    )
    results = process_file(args.csv_path, args.reference, config, sample_filter=args.sample)
    if not results:
        raise SystemExit(f"No samples matched filter: {args.sample}")
    summary_rows = []
    for result in results:
        summary_rows.append({
            "sample_name": result["sample_name"],
            "omega3_trio": result["omega3_trio"],
            "total_area": result["total_area"],
            "rt_shift": result["rt_shift"],
            "n_matched": int(result["matched_targets_df"]["area"].notna().sum()),
        })
        print(
            f"{result['sample_name']}: omega={result['omega3_trio']:.4f} "
            f"matched={summary_rows[-1]['n_matched']} rt_shift={result['rt_shift']:+.4f}"
        )

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(args.out) as writer:
            pd.DataFrame(summary_rows).to_excel(writer, sheet_name="Summary", index=False)
            for idx, result in enumerate(results, start=1):
                sheet_name = re.sub(r"[^A-Za-z0-9_]", "_", result["sample_name"])[:25] or f"sample_{idx}"
                result["matched_targets_df"].to_excel(writer, sheet_name=sheet_name, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
