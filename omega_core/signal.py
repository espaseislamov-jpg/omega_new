from __future__ import annotations

import math
import os
import threading
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import find_peaks, peak_widths, savgol_filter
from scipy.stats import median_abs_deviation

from . import rt_profile
from .io import get_runtime_app_dir

try:
    from pybaselines import Baseline
except Exception:
    Baseline = None

PEAK_RECORD_COLUMNS = [
    "peak_id",
    "start_idx",
    "apex_idx",
    "end_idx",
    "start_x",
    "apex_x",
    "end_x",
    "height",
    "prominence",
    "width_points",
    "area",
    "percent_area",
]

PEAK_INTEGRATION_REL_HEIGHT = 0.71
ARPLS_BASELINE_LAM = 1e8
WRITE_CHEBYSHEV_COEFFICIENTS = False
ENABLE_ARPLS_BASELINE_FALLBACK = True
ENABLE_ASLS_SHAPE_FALLBACK = True
ASLS_SHAPE_FALLBACK_LAM = 1e8
ASLS_SHAPE_FALLBACK_P = 0.01
CLUSTER_QUALITY_COMPLETE_SCORE = 50.0

BASELINE_KWARGS = {
    "degree": None,
    "n_bins": 300,
    "lower_quantile": 0.08,
    "n_iter": 10,
    "sigma_threshold": 0.7,
}
SAVGOL_POLYORDER = 3
SAVGOL_CANDIDATE_WINDOWS = [11, 15, 21, 31, 41, 51, 61, 81, 101, 151]
SAVGOL_MAX_SELECTED_WINDOW = 101
PEAK_BOUNDARY_SIGNAL_SMOOTH_WEIGHT = 0.70
PEAK_DETECTION_HEIGHT_SIGMA = 1.5
PEAK_DETECTION_PROMINENCE_SIGMA = 2.0
PEAK_BOUNDARY_MODE = os.environ.get("OMEGA_PEAK_BOUNDARY_MODE", "rel_height").strip().lower()
PEAK_PROMINENCE_BASE_MAX_WIDTH = 0.180
ENABLE_PYOPENMS_PEAK_ASSIST = True
PYOPENMS_GAUSS_WIDTH_SECONDS = 7.2
PYOPENMS_SIGNAL_TO_NOISE = 0.2
PYOPENMS_SN_WIN_LEN = 50.0
PYOPENMS_MIN_PROMINENCE_SIGMA = 0.75
PYOPENMS_MIN_PROMINENCE_FLOOR = 20.0

# ChemStation integrator defaults observed in the field setup.  Keep these as
# lower bounds for Omega's detector so small shoulders below the operator's
# threshold are not promoted into standalone C20/C22 target peaks.
CHEMSTATION_INITIAL_AREA_REJECT = 1.0
CHEMSTATION_INITIAL_PEAK_WIDTH = 0.016
CHEMSTATION_SHOULDER_DETECTION = False
CHEMSTATION_INITIAL_THRESHOLD = 12.9

_PYOPENMS_PEAK_PICKER = None
_PYOPENMS_IMPORT_ATTEMPTED = False
oms = None
_PYOPENMS_LOCK = threading.Lock()
_PYOPENMS_THREAD = threading.local()


def _robust_sigma(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return 0.0
    sigma = float(median_abs_deviation(arr, scale="normal", nan_policy="omit"))
    if sigma <= 0:
        sigma = float(np.std(arr))
    return sigma


def _get_x_column_name(df: pd.DataFrame) -> str:
    if "x_corrected" in df.columns:
        return "x_corrected"
    if "x" in df.columns:
        return "x"
    raise KeyError("DataFrame must contain 'x_corrected' or 'x'.")


def _extract_peak_geometry(df: pd.DataFrame, apex_idx: int, max_half_window_points: int = 160):
    x_col = _get_x_column_name(df)
    x = df[x_col].to_numpy(dtype=float)
    y_smooth = df["y_smooth"].to_numpy(dtype=float)
    y_corrected = df["y_corrected"].to_numpy(dtype=float)
    dy = df["dy"].to_numpy(dtype=float)
    boundary_signal = (
        PEAK_BOUNDARY_SIGNAL_SMOOTH_WEIGHT * y_smooth
        + (1.0 - PEAK_BOUNDARY_SIGNAL_SMOOTH_WEIGHT) * y_corrected
    )

    left_idx = int(apex_idx)
    steps = 0
    while left_idx > 1 and steps < max_half_window_points:
        if dy[left_idx - 1] <= 0 < dy[left_idx]:
            break
        left_idx -= 1
        steps += 1
    left_slice = slice(max(0, left_idx - 2), apex_idx + 1)
    left_local = np.argmin(boundary_signal[left_slice]) + left_slice.start
    left_idx = int(left_local)

    right_idx = int(apex_idx)
    steps = 0
    while right_idx < len(x) - 2 and steps < max_half_window_points:
        if dy[right_idx] < 0 <= dy[right_idx + 1]:
            break
        right_idx += 1
        steps += 1
    right_slice = slice(apex_idx, min(len(x), right_idx + 3))
    right_local = np.argmin(boundary_signal[right_slice]) + right_slice.start
    right_idx = int(right_local)

    if right_idx <= left_idx:
        left_idx = max(0, int(apex_idx) - 3)
        right_idx = min(len(x) - 1, int(apex_idx) + 3)
        if right_idx <= left_idx:
            return None

    local_floor = max(float(boundary_signal[left_idx]), float(boundary_signal[right_idx]))
    prominence = float(y_smooth[apex_idx] - local_floor)
    area = float(np.trapezoid(np.clip(y_corrected[left_idx:right_idx + 1], 0.0, None), x[left_idx:right_idx + 1]))
    return {
        "start_idx": left_idx,
        "apex_idx": int(apex_idx),
        "end_idx": right_idx,
        "start_x": float(x[left_idx]),
        "apex_x": float(x[apex_idx]),
        "end_x": float(x[right_idx]),
        "height": float(y_smooth[apex_idx]),
        "prominence": prominence,
        "width_points": float(right_idx - left_idx),
        "area": area,
    }


def _merge_peak_records(peaks_df: pd.DataFrame, extra_records) -> pd.DataFrame:
    if peaks_df is None or peaks_df.empty:
        base_records = []
    else:
        base_records = peaks_df.reindex(columns=PEAK_RECORD_COLUMNS).to_dict("records")
    if not extra_records:
        if not base_records:
            return pd.DataFrame(columns=PEAK_RECORD_COLUMNS)
        out = pd.DataFrame(base_records, columns=PEAK_RECORD_COLUMNS)
        return out.sort_values("apex_x").reset_index(drop=True)

    # Restore the pre-regression selection of complete detector records.
    # Keeping the primary unconditionally discarded the valley-based geometry.
    # A duplicate without valid geometry must never displace a usable peak.
    merged_records = list(base_records)
    for record in extra_records:
        item = {column: record.get(column, np.nan) for column in PEAK_RECORD_COLUMNS}
        apex_x = float(item["apex_x"])
        if not np.isfinite(apex_x):
            continue
        duplicates = [base for base in base_records
                      if abs(apex_x - float(base["apex_x"])) <= 0.006]
        if duplicates:
            start, end, area = (float(item[key]) for key in ("start_x", "end_x", "area"))
            if not (np.isfinite([start, end, area]).all() and start < apex_x < end and area > 0):
                continue
            if not any(start < float(base["apex_x"]) < end for base in duplicates):
                continue
            if any(start <= float(base["apex_x"]) <= end
                   for base in base_records
                   if abs(apex_x - float(base["apex_x"])) > 0.006):
                continue
        merged_records.append(item)

    merged_records.sort(key=lambda row: (float(row["apex_x"]), -float(row["area"])))
    deduped = []
    last_apex = None
    for row in merged_records:
        apex_x = float(row["apex_x"])
        if last_apex is not None and abs(apex_x - last_apex) <= 0.006:
            continue
        deduped.append(row)
        last_apex = apex_x

    if not deduped:
        return pd.DataFrame(columns=PEAK_RECORD_COLUMNS)

    out = pd.DataFrame(deduped, columns=PEAK_RECORD_COLUMNS)
    out["peak_id"] = np.arange(1, len(out) + 1)
    total_area = float(pd.to_numeric(out["area"], errors="coerce").fillna(0.0).sum())
    out["percent_area"] = 100.0 * pd.to_numeric(out["area"], errors="coerce") / total_area if total_area > 0 else np.nan
    return out


def _find_targeted_peak_candidate(
    df: pd.DataFrame,
    target_x: float,
    search_radius: float,
    min_prominence: float,
    min_area: float,
):
    x_col = _get_x_column_name(df)
    x = df[x_col].to_numpy(dtype=float)
    dy = df["dy"].to_numpy(dtype=float)
    min_prominence = max(float(min_prominence), CHEMSTATION_INITIAL_THRESHOLD)
    min_area = max(float(min_area), CHEMSTATION_INITIAL_AREA_REJECT)

    if len(x) < 3:
        return None

    zero_crossings = []
    for i in range(1, len(x)):
        if not (target_x - search_radius <= x[i] <= target_x + search_radius):
            continue
        if dy[i - 1] > 0 >= dy[i]:
            apex_idx = i - 1 if x[i - 1] >= target_x - search_radius else i
            zero_crossings.append(int(apex_idx))

    best = None
    for apex_idx in zero_crossings:
        geom = _extract_peak_geometry(df, apex_idx)
        if geom is None:
            continue
        if not CHEMSTATION_SHOULDER_DETECTION and geom["prominence"] < CHEMSTATION_INITIAL_THRESHOLD:
            continue
        if geom["prominence"] < min_prominence or geom["area"] < min_area:
            continue
        distance = abs(geom["apex_x"] - target_x)
        score = geom["prominence"] + 0.25 * geom["area"] - 2500.0 * distance
        if best is None or score > best["score"]:
            geom["score"] = float(score)
            best = geom
    return best


def augment_targeted_cluster_peaks(
    df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    reference_targets: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if df is None or df.empty:
        return peaks_df

    noise = max(_robust_sigma(df["y_corrected"].to_numpy(dtype=float)), 1.0)
    codes = [
        "C18:1N9C", "C18:3N3", "C18:0",
        "C20:4N6", "C20:5", "C20:3N8",
        "C22:6", "C22:5", "C22:4",
    ]
    fallback_centers = [7.623, 7.650, 7.750, 8.381, 8.410, 8.467, 9.252, 9.285, 9.316]
    centers = rt_profile.target_centers(
        reference_targets, codes, fallback_centers, corrected=False
    ) if reference_targets is not None else fallback_centers
    settings = [
        (0.018, max(8.0 * noise, 1200.0), 120.0),
        (0.020, max(6.0 * noise, 900.0), 120.0),
        (0.020, max(8.0 * noise, 1200.0), 120.0),
        (0.018, max(8.0 * noise, 1200.0), 120.0),
        (0.018, max(5.0 * noise, 700.0), 80.0),
        (0.020, max(5.0 * noise, 700.0), 80.0),
        (0.018, max(2.0 * noise, 250.0), 20.0),
        (0.018, max(2.0 * noise, 220.0), 20.0),
        (0.018, max(2.0 * noise, 220.0), 20.0),
    ]
    specs = [
        {"target_x": center, "search_radius": radius, "min_prominence": prominence, "min_area": area}
        for center, (radius, prominence, area) in zip(centers, settings)
    ]

    extra_records = []
    existing = peaks_df.copy() if peaks_df is not None else pd.DataFrame()
    for spec in specs:
        if not existing.empty and (existing["apex_x"] - spec["target_x"]).abs().min() <= 0.008:
            continue
        candidate = _find_targeted_peak_candidate(df, **spec)
        if candidate is not None:
            extra_records.append(candidate)
    return _merge_peak_records(peaks_df, extra_records)


def _get_pyopenms_peak_picker():
    cached = getattr(_PYOPENMS_THREAD, "picker", None)
    if cached is not None:
        return cached
    pyopenms = _load_pyopenms()
    if pyopenms is None:
        raise RuntimeError("PyOpenMS не загружен.")
    picker = pyopenms.PeakPickerChromatogram()
    params = picker.getParameters()
    params.setValue(b"gauss_width", float(PYOPENMS_GAUSS_WIDTH_SECONDS))
    params.setValue(b"signal_to_noise", float(PYOPENMS_SIGNAL_TO_NOISE))
    params.setValue(b"sn_win_len", float(PYOPENMS_SN_WIN_LEN))
    params.setValue(b"use_gauss", b"true")
    params.setValue(b"remove_overlapping_peaks", b"false")
    picker.setParameters(params)
    _PYOPENMS_THREAD.picker = picker
    return picker


def _load_pyopenms():
    global _PYOPENMS_IMPORT_ATTEMPTED, oms
    with _PYOPENMS_LOCK:
        if oms is not None:
            return oms
        try:
            import pyopenms as loaded_oms
        except Exception as exc:
            raise RuntimeError("PyOpenMS не загрузился. Расчёт остановлен: отключать детектор молча нельзя.") from exc
        oms = loaded_oms
        _PYOPENMS_IMPORT_ATTEMPTED = True
        return oms


def detect_peaks_with_pyopenms(df: pd.DataFrame) -> pd.DataFrame:
    if not ENABLE_PYOPENMS_PEAK_ASSIST or df is None or df.empty:
        return pd.DataFrame(columns=PEAK_RECORD_COLUMNS)
    pyopenms = _load_pyopenms()
    if pyopenms is None:
        raise RuntimeError("PyOpenMS не загружен.")

    try:
        x_col = _get_x_column_name(df)
        x = df[x_col].to_numpy(dtype=float)
        y_corrected = np.clip(df["y_corrected"].to_numpy(dtype=float), 0.0, None)
        if x.size < 20 or not np.any(y_corrected > 0):
            return pd.DataFrame(columns=PEAK_RECORD_COLUMNS)

        chromatogram = pyopenms.MSChromatogram()
        chromatogram.set_peaks((x * 60.0, y_corrected))
        picked = pyopenms.MSChromatogram()
        picker = _get_pyopenms_peak_picker()
        if picker is None:
            raise RuntimeError("Детектор PyOpenMS не создан.")
        picker.pickChromatogram(chromatogram, picked)
        peak_rts, _ = picked.get_peaks()
    except Exception as exc:
        raise RuntimeError("Ошибка детектора PyOpenMS. Пустой список пиков не является результатом расчёта.") from exc

    prominence_floor = max(
        PYOPENMS_MIN_PROMINENCE_FLOOR,
        _robust_sigma(y_corrected) * PYOPENMS_MIN_PROMINENCE_SIGMA,
    )
    extra_records = []
    for peak_rt_seconds in peak_rts:
        peak_rt = float(peak_rt_seconds) / 60.0
        apex_idx = int(np.argmin(np.abs(x - peak_rt)))
        geom = _extract_peak_geometry(df, apex_idx)
        if geom is None:
            continue
        if float(geom["prominence"]) < prominence_floor or float(geom["area"]) <= 0.0:
            continue
        extra_records.append({
            "start_idx": int(geom["start_idx"]),
            "apex_idx": int(geom["apex_idx"]),
            "end_idx": int(geom["end_idx"]),
            "start_x": float(geom["start_x"]),
            "apex_x": float(geom["apex_x"]),
            "end_x": float(geom["end_x"]),
            "height": float(geom["height"]),
            "prominence": float(geom["prominence"]),
            "width_points": float(geom["width_points"]),
            "area": float(geom["area"]),
        })
    return _merge_peak_records(pd.DataFrame(columns=PEAK_RECORD_COLUMNS), extra_records)


def _fit_chebyshev_baseline(
    x: np.ndarray,
    y: np.ndarray,
    degree: int,
    n_bins: int,
    lower_quantile: float,
    n_iter: int,
    sigma_threshold: float,
) -> tuple[np.ndarray, np.ndarray]:
    if x.size < 8:
        baseline = np.full_like(y, np.quantile(y, lower_quantile))
        return baseline, np.array([float(np.median(y))], dtype=float)

    x_scaled = np.interp(x, (x.min(), x.max()), (-1.0, 1.0))
    bin_edges = np.linspace(x.min(), x.max(), num=max(16, min(n_bins, x.size // 4)) + 1)
    anchor_x = []
    anchor_y = []
    for left, right in zip(bin_edges[:-1], bin_edges[1:]):
        mask = (x >= left) & (x < right if right < bin_edges[-1] else x <= right)
        if not mask.any():
            continue
        x_bin = x[mask]
        y_bin = y[mask]
        anchor_x.append(float(np.median(x_bin)))
        anchor_y.append(float(np.quantile(y_bin, lower_quantile)))

    if len(anchor_x) < degree + 1:
        step = max(1, x.size // max(degree + 2, 16))
        anchor_x = x[::step].tolist()
        anchor_y = np.quantile(y.reshape(-1, 1), lower_quantile, axis=1)[::step].tolist()

    anchor_x = np.asarray(anchor_x, dtype=float)
    anchor_y = np.asarray(anchor_y, dtype=float)
    fit_mask = np.ones_like(y, dtype=bool)
    coeffs = np.zeros(degree + 1, dtype=float)
    baseline = np.full_like(y, np.median(anchor_y) if anchor_y.size else np.median(y))

    for _ in range(max(1, int(n_iter))):
        fit_x = np.concatenate([anchor_x, x[fit_mask]])
        fit_y = np.concatenate([anchor_y, y[fit_mask]])
        if fit_x.size < degree + 1:
            break
        fit_x_scaled = np.interp(fit_x, (x.min(), x.max()), (-1.0, 1.0))
        coeffs = np.polynomial.chebyshev.chebfit(fit_x_scaled, fit_y, deg=degree)
        baseline = np.polynomial.chebyshev.chebval(x_scaled, coeffs)
        residual = y - baseline
        sigma = _robust_sigma(residual)
        if sigma <= 0:
            break
        fit_mask = residual <= sigma_threshold * sigma

    return baseline, coeffs


def add_baseline(
    df: pd.DataFrame,
    degree=None,
    n_bins: int = 300,
    lower_quantile: float = 0.08,
    n_iter: int = 10,
    sigma_threshold: float = 0.7,
) -> pd.DataFrame:
    out = df.copy()
    x = out["x_corrected"].to_numpy(dtype=float)
    y = out["y"].to_numpy(dtype=float)
    if x.size < 8:
        raise ValueError("Not enough points for baseline correction.")

    resolved_degree = 6 if degree is None else int(degree)
    baseline, coeffs = _fit_chebyshev_baseline(x, y, resolved_degree, n_bins, lower_quantile, n_iter, sigma_threshold)
    out["baseline"] = baseline
    out["y_corrected"] = y - baseline

    if WRITE_CHEBYSHEV_COEFFICIENTS:
        coeff_path = Path(get_runtime_app_dir()) / "chebyshev_coefficients.csv"
        pd.DataFrame({
            "coefficient_index": np.arange(len(coeffs), dtype=int),
            "coefficient": coeffs,
        }).to_csv(coeff_path, index=False)
    return out


def add_arpls_baseline(df: pd.DataFrame, lam: float = ARPLS_BASELINE_LAM) -> pd.DataFrame:
    if Baseline is None:
        return add_baseline(df, **BASELINE_KWARGS)

    out = df.copy()
    y = out["y"].to_numpy(dtype=float)
    if y.size < 8:
        raise ValueError("Not enough points for baseline correction.")

    baseline, _ = Baseline().arpls(y, lam=float(lam))
    baseline = np.asarray(baseline, dtype=float)
    out["baseline"] = baseline
    out["y_corrected"] = y - baseline
    return out


def add_asls_baseline(
    df: pd.DataFrame,
    lam: float = ASLS_SHAPE_FALLBACK_LAM,
    p: float = ASLS_SHAPE_FALLBACK_P,
) -> pd.DataFrame:
    if Baseline is None:
        return add_baseline(df, **BASELINE_KWARGS)

    out = df.copy()
    y = out["y"].to_numpy(dtype=float)
    if y.size < 8:
        raise ValueError("Not enough points for baseline correction.")

    baseline, _ = Baseline().asls(y, lam=float(lam), p=float(p))
    baseline = np.asarray(baseline, dtype=float)
    out["baseline"] = baseline
    out["y_corrected"] = y - baseline
    return out


def add_smoothing_and_derivatives(
    df: pd.DataFrame,
    polyorder: int = SAVGOL_POLYORDER,
    candidate_windows=None,
) -> tuple[pd.DataFrame, int]:
    """Derivatives in physical time; smoothing is internal, never the display signal."""
    from .peak_geometry import derivative_signal
    out = df.copy()
    x = out[_get_x_column_name(out)].to_numpy(dtype=float)
    y = out['y_corrected'].to_numpy(dtype=float)
    gx, _, smooth, dy, d2y, window = derivative_signal(x, y)
    for name, values in (('y_smooth', smooth), ('dy', dy), ('d2y', d2y)):
        out[name] = np.interp(x, gx, values)
    return out, window


def detect_peak_candidates(
    df: pd.DataFrame,
    best_window=None,
    reference_targets: pd.DataFrame | None = None,
    height_sigma: float = PEAK_DETECTION_HEIGHT_SIGMA,
    prominence_sigma: float = PEAK_DETECTION_PROMINENCE_SIGMA,
    rel_height: float = PEAK_INTEGRATION_REL_HEIGHT,
) -> pd.DataFrame:
    """Detect geometry independently of target names and legacy integration widths."""
    from .peak_geometry import detect_geometry
    if df is None or df.empty:
        return pd.DataFrame(columns=PEAK_RECORD_COLUMNS)
    return detect_geometry(df[_get_x_column_name(df)].to_numpy(dtype=float),
                           df['y_corrected'].to_numpy(dtype=float))


def _augment_detected_peaks(df, peaks_df, reference_targets):
    """Recovery must also run when the primary detector rejects every peak."""
    peaks_df = augment_targeted_cluster_peaks(df, peaks_df, reference_targets=reference_targets)
    pyopenms_peaks_df = detect_peaks_with_pyopenms(df)
    if not pyopenms_peaks_df.empty:
        peaks_df = _merge_peak_records(peaks_df, pyopenms_peaks_df.to_dict("records"))
    return peaks_df
