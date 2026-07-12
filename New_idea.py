import csv
import itertools
import json
import math
import re
import shutil
import sys
import warnings
from pathlib import Path

from omega_path_compat import configure_windows_path_compat

configure_windows_path_compat()

import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog, ttk, messagebox

import matplotlib
matplotlib.use("TkAgg")

from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from scipy.optimize import least_squares
from scipy.signal import find_peaks, peak_widths, savgol_filter
from scipy.stats import median_abs_deviation

import omega_core
import omega_chromatopy_clean
from omega_core import io as core_io
from omega_core import metrics as core_metrics
from omega_core import signal as core_signal

try:
    import pyopenms as oms
except Exception:
    oms = None

try:
    from pybaselines import Baseline
except Exception:
    Baseline = None

try:
    from lmfit.models import LinearModel, PseudoVoigtModel
except Exception:
    LinearModel = None
    PseudoVoigtModel = None


_PYOPENMS_PEAK_PICKER = None


DEFAULT_REFERENCE_TARGETS = [
    {"component": "C16:1N7", "code": "C16:1N7", "display_name": "Пальмитолениновая", "order_index": 1, "expected_rt": 6.6, "rt_reliable": True, "historical_area": 77164.071881, "historical_percent": 0.623009, "notes": ""},
    {"component": "C16:0", "code": "C16:0", "display_name": "Пальмитиновая", "order_index": 2, "expected_rt": 6.7, "rt_reliable": True, "historical_area": 2523674.606557, "historical_percent": 20.922966, "notes": ""},
    {"component": "C18:3N6", "code": "C18:3N6", "display_name": "γ-Линоленовая", "order_index": 3, "expected_rt": 7.5, "rt_reliable": True, "historical_area": 85833.037901, "historical_percent": 0.201383, "notes": ""},
    {"component": "C18:2N6C", "code": "C18:2N6C", "display_name": "Линолевая", "order_index": 4, "expected_rt": 7.6, "rt_reliable": True, "historical_area": 2625136.957945, "historical_percent": 19.869624, "notes": ""},
    {"component": "C18:1N9C", "code": "C18:1N9C", "display_name": "Олеиновая", "order_index": 5, "expected_rt": 7.7, "rt_reliable": False, "historical_area": 1881004.860367, "historical_percent": 14.606759, "notes": "RT использовать мягко; рядом C18:3n3."},
    {"component": "C18:3N3", "code": "C18:3N3", "display_name": "Линоленовая", "order_index": 6, "expected_rt": None, "rt_reliable": False, "historical_area": 0.0, "historical_percent": 1.379643, "notes": "RT отсутствует; идентификация по порядку."},
    {"component": "C18:0", "code": "C18:0", "display_name": "Стеариновая", "order_index": 7, "expected_rt": 7.8, "rt_reliable": False, "historical_area": 1688262.321032, "historical_percent": 13.215531, "notes": "RT дублируется; идентификация order-first."},
    {"component": "C20:4N6", "code": "C20:4N6", "display_name": "Арахидоновая", "order_index": 8, "expected_rt": 8.4, "rt_reliable": False, "historical_area": 2003832.525827, "historical_percent": 13.012492, "notes": "RT дублируется; идентификация order-first."},
    {"component": "C20:5", "code": "C20:5", "display_name": "Эйкозопентаеновая", "order_index": 9, "expected_rt": 8.4, "rt_reliable": False, "historical_area": 134069.083312, "historical_percent": 0.993616, "notes": "RT дублируется; идентификация order-first."},
    {"component": "C20:3N8", "code": "C20:3N8", "display_name": "цис-8,11,14-эйкозатриеновая 320", "order_index": 10, "expected_rt": 7.8, "rt_reliable": False, "historical_area": 184300.772567, "historical_percent": 1.016786, "notes": "Фактически идет после C20:5; RT ненадёжен."},
    {"component": "C22:6", "code": "C22:6", "display_name": "Докозагексаеновая", "order_index": 11, "expected_rt": 9.2, "rt_reliable": False, "historical_area": 609757.540544, "historical_percent": 4.262376, "notes": "Группа C22; брать 2-й пик в группе."},
    {"component": "C22:5", "code": "C22:5", "display_name": "Docosapentaenoic acid", "order_index": 12, "expected_rt": 9.2, "rt_reliable": False, "historical_area": 186274.058539, "historical_percent": 1.302286, "notes": "Группа C22; DPA идёт вторым после C22:6, брать 3-й пик в группе."},
    {"component": "C22:4", "code": "C22:4", "display_name": "Докозатетраеновая", "order_index": 13, "expected_rt": 9.2, "rt_reliable": False, "historical_area": 218113.146615, "historical_percent": 1.505486, "notes": "Группа C22; идёт после C22:5, брать 4-й пик в группе."},
    {"component": "C24:1N9", "code": "C24:1N9", "display_name": "Нервоновая", "order_index": 14, "expected_rt": 10.3, "rt_reliable": True, "historical_area": 482256.655642, "historical_percent": 3.400287, "notes": ""},
    {"component": "C24:0", "code": "C24:0", "display_name": "Лигноцериновая", "order_index": 15, "expected_rt": 10.4, "rt_reliable": True, "historical_area": 544993.843627, "historical_percent": 3.687757, "notes": ""},
]

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

PEAK_INTEGRATION_REL_HEIGHT = 0.68
PEAK_SUPPORT_THRESHOLD_SIGMA = 0.80
PEAK_SUPPORT_THRESHOLD_FRACTION = 0.020
PEAK_SUPPORT_CONSECUTIVE_POINTS = 4
RELIABLE_RT_WINDOW = 0.035
RELIABLE_RT_DOMINANT_DISTANCE_MAX = 0.025
RELIABLE_RT_DOMINANT_AREA_MULTIPLIER = 8.0
RELIABLE_RT_DOMINANT_AREA_MIN_DELTA = 400.0
C22_OVERLAP_TRIGGER_OMEGA_MIN = 4.0
C22_OVERLAP_RATIO_OFFSET = 1.25
C22_OVERLAP_RATIO_SLOPE = 1.0
C22_OVERLAP_FRACTION_CAP = 0.95
ENABLE_DATA_DRIVEN_C22_OVERLAP_MODEL = True
C22_OVERLAP_MODEL_BLEND = 0.60
C22_OVERLAP_MODEL_APPLY_FRACTION_MIN = 0.90
C22_OVERLAP_WIDE_CLUSTER_MEAN_WIDTH = 0.030
C22_OVERLAP_WIDE_CLUSTER_SCALE = 0.65
ENABLE_C22_TAIL_TIGHTENING = True
C22_TAIL_TIGHTENING_MEAN_WIDTH = 0.036
C22_TAIL_TIGHTENING_WIDTH_SCALE = 0.95
C22_TAIL_TIGHTENING_AREA_RATIO_MIN = 0.92
C22_TAIL_TIGHTENING_AREA_RATIO_MAX = 0.995
C22_TAIL_TIGHTENING_DPA_RATIO_TRIGGER = 0.90
C22_TAIL_TIGHTENING_DHA_WIDTH_TRIGGER = 0.040
ENABLE_OVERWIDE_C22_PVFIT_REFINEMENT = True
C22_PVFIT_OVERWIDE_MEAN_WIDTH_MIN = 0.032
C22_PVFIT_OVERWIDE_DHA_WIDTH_MIN = 0.039
C22_PVFIT_OVERWIDE_C22_4_WIDTH_MIN = 0.033
C22_PVFIT_AREA_RATIO_MIN = 0.92
C22_PVFIT_AREA_RATIO_MAX = 0.995
C18_DENOMINATOR_DOMINANCE_RATIO = 1.60
C18_DENOMINATOR_SMALL_N3_FRACTION = 0.08
C18_DENOMINATOR_AREA_SCALE = 0.90
C18_DENOMINATOR_EXTREME_RATIO = 2.20
C18_DENOMINATOR_EXTREME_SMALL_N3_FRACTION = 0.03
C18_DENOMINATOR_EXTREME_WIDTH_MIN = 0.043
C18_DENOMINATOR_EXTREME_STRICT_MAX = 6.0
C18_DENOMINATOR_EXTREME_AREA_SCALE = 0.70
C22_OVERLAP_MODEL_SCALES = np.asarray([
    1.9632271836116177,
    0.542462045184553,
    0.5456616773085775,
    0.007115503837232517,
    0.0038617345289501124,
    0.010517989312295157,
    0.36083424193869684,
    0.38394711984974766,
], dtype=float)
C22_OVERLAP_MODEL_PARAMS = np.asarray([
    5.1693037865961315,
    0.06019002814508915,
    -2.076049325959199,
    0.5332510263009543,
    -1.0580072082286593,
    -0.1927639442183294,
    -0.9528310466284664,
    -0.5952360000557718,
    0.27922158744448305,
], dtype=float)
ENABLE_DATA_DRIVEN_C20_EPA_MODEL = True
C20_EPA_MODEL_GATE_RATIO_MAX = 0.70
C20_EPA_MODEL_BLEND = 0.25
C20_EPA_OVERLAP_WIDE_NEIGHBOR_RATIO = 1.30
C20_EPA_OVERLAP_EXTRA_SCALE = 1.60
C20_EPA_UNDERFIT_RATIO_MAX = 0.13
C20_EPA_UNDERFIT_WIDTH_RATIO = 1.80
C20_EPA_UNDERFIT_STRICT_MAX = 4.20
C20_EPA_UNDERFIT_CREDIT_MIN = 20.0
C20_EPA_UNDERFIT_EXTRA_SCALE = 1.50
C20_EPA_MODEL_SCALES = np.asarray([
    0.9895635778699382,
    0.5615594221089807,
    0.3069731180293255,
    0.3211440874041738,
    0.003591532739879512,
    0.002105803315168974,
    0.005357045460285181,
], dtype=float)
C20_EPA_MODEL_PARAMS = np.asarray([
    -29.84201291974811,
    -0.5346083467980479,
    -2.321206829880846,
    1.5284999270125061,
    -0.8310768530097283,
    1.6002603944153417,
    0.6483089815285562,
    -0.49967132476694437,
], dtype=float)
SMALL_PEAK_SHARP_SEARCH_HALF_WINDOW = 0.070
SMALL_PEAK_SHARP_APEX_SEARCH_RADIUS = 0.012
SMALL_PEAK_SHARP_THRESHOLD_FRACTION = 0.08
SMALL_PEAK_SHARP_THRESHOLD_SIGMA = 1.15
SMALL_PEAK_SHARP_MAX_HALF_WIDTH = 0.024
SMALL_PEAK_SHARP_MAX_ASYMMETRY = 1.35
SMALL_PEAK_SHARP_MIN_AREA_RATIO = 0.55
SMALL_PEAK_SHARP_MAX_AREA_RATIO = 0.98
SMALL_PEAK_SHARP_MAX_PERCENT_AREA = 1.50
SMALL_PEAK_SHARP_SPECS = {
    "C18:3N6": {
        "mode": "isolated",
        "max_percent": SMALL_PEAK_SHARP_MAX_PERCENT_AREA,
        "min_area_ratio": SMALL_PEAK_SHARP_MIN_AREA_RATIO,
        "max_area_ratio": SMALL_PEAK_SHARP_MAX_AREA_RATIO,
        "max_width_ratio": 0.96,
        "threshold_fraction": SMALL_PEAK_SHARP_THRESHOLD_FRACTION,
        "threshold_sigma": SMALL_PEAK_SHARP_THRESHOLD_SIGMA,
    },
    "C20:5": {
        "mode": "bounded",
        "max_percent": 1.80,
        "min_area_ratio": 0.82,
        "max_area_ratio": 0.98,
        "max_width_ratio": 0.90,
        "threshold_fraction": 0.11,
        "threshold_sigma": 1.15,
        "min_asymmetry": 2.80,
    },
    "C22:5": {
        "mode": "bounded",
        "max_percent": 1.20,
        "min_area_ratio": 0.92,
        "max_area_ratio": 0.99,
        "max_width_ratio": 0.92,
        "threshold_fraction": 0.10,
        "threshold_sigma": 1.10,
        "min_asymmetry": 1.15,
    },
}
C18_OVERLAP_START_TOLERANCE = 0.001
C20_LOCAL_AREA_RATIO_TRIGGER = 1.05
C20_LOCAL_BOUNDARY_EXTENSION = 0.004
C20_FIT_EPA_AREA_MAX = 450.0
C20_FIT_EPA_PROMINENCE_MAX = 1500.0
ENABLE_PYOPENMS_PEAK_ASSIST = True
PYOPENMS_GAUSS_WIDTH_SECONDS = 7.2
PYOPENMS_SIGNAL_TO_NOISE = 0.2
PYOPENMS_SN_WIN_LEN = 50.0
PYOPENMS_MIN_PROMINENCE_SIGMA = 0.75
PYOPENMS_MIN_PROMINENCE_FLOOR = 20.0
ENABLE_LMFIT_LOCAL_PSEUDOVOIGT = True
LMFIT_LOCAL_PSEUDOVOIGT_MIN_R2 = 0.82
LMFIT_LOCAL_AREA_RATIO_MIN = 0.65
LMFIT_LOCAL_AREA_RATIO_MAX = 1.35
ENABLE_LMFIT_C18_RECOVERY = False
ENABLE_SPLIT_PSEUDOVOIGT_CLUSTER_FIT = False
SPLIT_PSEUDOVOIGT_MIN_R2 = 0.84
SPLIT_PSEUDOVOIGT_AREA_RATIO_MIN = 0.72
SPLIT_PSEUDOVOIGT_AREA_RATIO_MAX = 1.42
SPLIT_PSEUDOVOIGT_ASYMMETRY_MIN = 0.55
SPLIT_PSEUDOVOIGT_ASYMMETRY_MAX = 1.90
SPLIT_PSEUDOVOIGT_BOUNDARY_SNAP_WINDOW = 0.012
SPLIT_PSEUDOVOIGT_OUTER_SUPPORT_WIDTH_FACTOR = 2.65
SPLIT_PSEUDOVOIGT_VALLEY_WEIGHT_BOOST = 0.90
SPLIT_PSEUDOVOIGT_FOOT_WEIGHT_BOOST = 0.55
SPLIT_PSEUDOVOIGT_EDGE_WEIGHT_BOOST = 0.20
ENABLE_ARPLS_BASELINE_FALLBACK = True
ARPLS_BASELINE_LAM = 1e8
WRITE_CHEBYSHEV_COEFFICIENTS = False
CLUSTER_QUALITY_COMPLETE_SCORE = 50.0
PREVIEW_WINDOWS = [
    ("6.0-7.3", 6.0, 7.3),
    ("7.4-7.7", 7.4, 7.7),
    ("8.3-8.7", 8.3, 8.7),
    ("9.1-9.4", 9.1, 9.4),
]
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
CLUSTER_BOUNDARY_METRIC_SMOOTH_WEIGHT = 0.60
PEAK_DETECTION_HEIGHT_SIGMA = 1.5
PEAK_DETECTION_PROMINENCE_SIGMA = 2.0


def get_runtime_resource_dir() -> Path:
    if getattr(sys, "frozen", False):
        meipass = getattr(sys, "_MEIPASS", None)
        if meipass:
            return Path(meipass)
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parent


def get_runtime_app_dir() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return Path(__file__).resolve().parent


def ensure_runtime_file(file_name: str) -> Path:
    app_dir = get_runtime_app_dir()
    app_dir.mkdir(parents=True, exist_ok=True)
    target_path = app_dir / file_name
    if target_path.exists():
        return target_path

    source_path = get_runtime_resource_dir() / file_name
    if source_path.exists():
        try:
            shutil.copy2(source_path, target_path)
        except Exception:
            pass
    return target_path


def init_reference_json(reference_json_path: Path):
    if not reference_json_path.exists():
        with reference_json_path.open("w", encoding="utf-8") as f:
            json.dump(DEFAULT_REFERENCE_TARGETS, f, ensure_ascii=False, indent=2)


def load_reference_json(reference_json_path: Path) -> pd.DataFrame:
    init_reference_json(reference_json_path)
    with reference_json_path.open("r", encoding="utf-8") as f:
        raw_targets = json.load(f)

    if not isinstance(raw_targets, list):
        raise ValueError("Reference JSON must contain a list of targets.")

    records = []
    for idx, item in enumerate(raw_targets, start=1):
        if not isinstance(item, dict):
            continue
        record = {
            "component": item.get("component", item.get("code", f"target_{idx}")),
            "code": item.get("code", item.get("component", f"target_{idx}")),
            "display_name": item.get("display_name", item.get("component", item.get("code", f"Target {idx}"))),
            "order_index": item.get("order_index", idx),
            "expected_rt": item.get("expected_rt", item.get("target_rt")),
            "rt_reliable": bool(item.get("rt_reliable", False)),
            "historical_area": item.get("historical_area", 0.0),
            "historical_percent": item.get("historical_percent", 0.0),
            "notes": item.get("notes", ""),
        }
        records.append(record)

    df = pd.DataFrame(records)
    if df.empty:
        raise ValueError("Reference target list is empty.")

    order_index = pd.to_numeric(df["order_index"], errors="coerce")
    order_index = order_index.where(order_index.notna(), pd.Series(np.arange(1, len(df) + 1), index=df.index, dtype=float))
    df["order_index"] = order_index.astype(int)
    df["expected_rt"] = pd.to_numeric(df["expected_rt"], errors="coerce")
    df["historical_area"] = pd.to_numeric(df["historical_area"], errors="coerce").fillna(0.0)
    df["historical_percent"] = pd.to_numeric(df["historical_percent"], errors="coerce").fillna(0.0)
    df["rt_reliable"] = df["rt_reliable"].fillna(False).astype(bool)
    return df.sort_values("order_index").reset_index(drop=True)


def _robust_sigma(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return 0.0
    sigma = float(median_abs_deviation(arr, scale="normal", nan_policy="omit"))
    if sigma <= 0:
        sigma = float(np.std(arr))
    return float(sigma)


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


def _extend_boundary_to_support(
    signal: np.ndarray,
    start_idx: int,
    limit_idx: int,
    direction: int,
    threshold: float,
    consecutive_points: int,
) -> int:
    if direction < 0:
        start_idx = int(max(limit_idx, start_idx))
        for idx in range(start_idx, int(limit_idx) - 1, -1):
            seg_start = max(int(limit_idx), idx - int(consecutive_points) + 1)
            if np.all(signal[seg_start:idx + 1] <= threshold):
                return int(idx)
        return int(start_idx)

    start_idx = int(min(limit_idx, start_idx))
    for idx in range(start_idx, int(limit_idx) + 1):
        seg_end = min(int(limit_idx) + 1, idx + int(consecutive_points))
        if np.all(signal[idx:seg_end] <= threshold):
            return int(idx)
    return int(start_idx)


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
        if geom["prominence"] < min_prominence or geom["area"] < min_area:
            continue
        distance = abs(geom["apex_x"] - target_x)
        score = geom["prominence"] + 0.25 * geom["area"] - 2500.0 * distance
        if best is None or score > best["score"]:
            geom["score"] = float(score)
            best = geom
    return best


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

    merged_records = list(base_records)
    for record in extra_records:
        item = {column: record.get(column, np.nan) for column in PEAK_RECORD_COLUMNS}
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


def _recompute_matched_percent_area(matched_targets_df: pd.DataFrame) -> pd.DataFrame:
    out = matched_targets_df.copy()
    total_area = float(pd.to_numeric(out["area"], errors="coerce").fillna(0.0).sum())
    if total_area > 0:
        out["percent_area"] = 100.0 * pd.to_numeric(out["area"], errors="coerce") / total_area
    else:
        out["percent_area"] = np.nan
    return out


def _collect_local_cluster_peak_geometries(
    df: pd.DataFrame,
    window_left: float,
    window_right: float,
    min_prominence: float,
    min_area: float,
    dedupe_distance: float = 0.004,
) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame()

    x_col = _get_x_column_name(df)
    x = df[x_col].to_numpy(dtype=float)
    dy = df["dy"].to_numpy(dtype=float)
    records = []
    for i in range(1, len(x) - 1):
        if x[i] < window_left or x[i] > window_right:
            continue
        if dy[i - 1] > 0 >= dy[i]:
            geom = _extract_peak_geometry(df, i)
            if geom is None:
                continue
            if geom["prominence"] < min_prominence or geom["area"] < min_area:
                continue
            records.append(geom)

    if not records:
        return pd.DataFrame()

    ordered = sorted(records, key=lambda row: row["apex_x"])
    deduped = []
    for row in ordered:
        if deduped and abs(float(row["apex_x"]) - float(deduped[-1]["apex_x"])) <= dedupe_distance:
            if float(row["prominence"]) > float(deduped[-1]["prominence"]):
                deduped[-1] = row
            continue
        deduped.append(row)
    return pd.DataFrame(deduped)


def _select_ordered_cluster_peaks(
    candidates_df: pd.DataFrame,
    target_apexes,
    max_distances,
):
    if candidates_df is None or candidates_df.empty:
        return None

    candidates = candidates_df.sort_values("apex_x").reset_index(drop=True)
    target_apexes = [float(value) for value in target_apexes]
    max_distances = [float(value) for value in max_distances]
    if len(candidates) < len(target_apexes):
        return None

    best_choice = None
    for combo in itertools.combinations(range(len(candidates)), len(target_apexes)):
        chosen = candidates.iloc[list(combo)].copy().reset_index(drop=True)
        distances = [abs(float(chosen.iloc[i]["apex_x"]) - target_apexes[i]) for i in range(len(target_apexes))]
        if any(distance > max_distances[i] for i, distance in enumerate(distances)):
            continue
        score = (
            float(sum(distances))
            - 1e-6 * float(chosen["prominence"].sum())
            - 1e-7 * float(chosen["area"].sum())
        )
        if best_choice is None or score < best_choice[0]:
            best_choice = (score, chosen)

    return None if best_choice is None else best_choice[1]


def _cluster_has_integration_overlap(matched_targets_df: pd.DataFrame, cluster_codes) -> bool:
    cluster = matched_targets_df[matched_targets_df["code"].isin(cluster_codes)].copy()
    if cluster.empty:
        return False

    cluster["found_rt"] = pd.to_numeric(cluster["found_rt"], errors="coerce")
    cluster["integration_start_x"] = pd.to_numeric(cluster["integration_start_x"], errors="coerce")
    cluster["integration_end_x"] = pd.to_numeric(cluster["integration_end_x"], errors="coerce")
    cluster = cluster.dropna(subset=["found_rt"]).sort_values("found_rt")

    previous_end = None
    for _, row in cluster.iterrows():
        start_x = row.get("integration_start_x")
        end_x = row.get("integration_end_x")
        if previous_end is not None and np.isfinite(start_x) and start_x < previous_end - C18_OVERLAP_START_TOLERANCE:
            return True
        if np.isfinite(end_x):
            previous_end = float(end_x)
    return False


def _cluster_has_duplicate_peak_ids(matched_targets_df: pd.DataFrame, cluster_codes) -> bool:
    cluster = matched_targets_df[matched_targets_df["code"].isin(cluster_codes)].copy()
    if cluster.empty:
        return False
    peak_ids = pd.to_numeric(cluster["matched_peak_id"], errors="coerce").dropna().astype(int)
    return bool(peak_ids.duplicated().any())


def _attach_local_peak_records(peaks_df: pd.DataFrame, chosen_geometries: pd.DataFrame) -> pd.DataFrame:
    if chosen_geometries is None or chosen_geometries.empty:
        return peaks_df.copy()

    existing = peaks_df.copy()
    extra_records = []
    for _, geom in chosen_geometries.iterrows():
        if not existing.empty and (existing["apex_x"] - float(geom["apex_x"])).abs().min() <= 0.006:
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
    return _merge_peak_records(existing, extra_records)


def _assign_local_geometry_to_row(out: pd.DataFrame, row_idx: int, geom: pd.Series, status: str):
    out.at[row_idx, "found_rt"] = float(geom["apex_x"])
    out.at[row_idx, "area"] = float(geom["area"])
    out.at[row_idx, "integration_start_x"] = float(geom["start_x"])
    out.at[row_idx, "integration_end_x"] = float(geom["end_x"])
    out.at[row_idx, "status"] = status
    out.at[row_idx, "match_score"] = np.nan
    out.at[row_idx, "matched_peak_id"] = np.nan


def _assign_local_geometry_bounds_to_row(out: pd.DataFrame, row_idx: int, geom: pd.Series, status: str):
    out.at[row_idx, "found_rt"] = float(geom["apex_x"])
    out.at[row_idx, "integration_start_x"] = float(geom["start_x"])
    out.at[row_idx, "integration_end_x"] = float(geom["end_x"])
    out.at[row_idx, "status"] = status
    out.at[row_idx, "match_score"] = np.nan


def _append_status_suffix(status_value, suffix: str) -> str:
    status_text = str(status_value or "").strip()
    if not status_text:
        return suffix
    if status_text.endswith(f"_{suffix}") or status_text == suffix:
        return status_text
    return f"{status_text}_{suffix}"


def _estimate_local_linear_baseline(
    x: np.ndarray,
    y: np.ndarray,
    start_idx: int,
    end_idx: int,
    edge_fraction: float = 0.16,
):
    start_idx = int(max(0, start_idx))
    end_idx = int(min(len(x) - 1, end_idx))
    if end_idx <= start_idx:
        return np.zeros(0, dtype=float)

    x_seg = np.asarray(x[start_idx:end_idx + 1], dtype=float)
    y_seg = np.asarray(y[start_idx:end_idx + 1], dtype=float)
    if x_seg.size <= 4:
        edge_line = np.linspace(float(y_seg[0]), float(y_seg[-1]), len(y_seg))
        return np.asarray(edge_line, dtype=float)

    edge_count = max(3, min(len(x_seg) // 2, int(round(len(x_seg) * edge_fraction))))
    left_x = x_seg[:edge_count]
    right_x = x_seg[-edge_count:]
    left_y = y_seg[:edge_count]
    right_y = y_seg[-edge_count:]

    left_anchor_x = float(np.mean(left_x))
    right_anchor_x = float(np.mean(right_x))
    left_anchor_y = float(np.quantile(left_y, 0.20))
    right_anchor_y = float(np.quantile(right_y, 0.20))
    if right_anchor_x <= left_anchor_x + 1e-9:
        return np.full(len(x_seg), min(left_anchor_y, right_anchor_y), dtype=float)
    slope = (right_anchor_y - left_anchor_y) / (right_anchor_x - left_anchor_x)
    intercept = left_anchor_y - slope * left_anchor_x
    return intercept + slope * x_seg


def _find_preferred_minimum_index(
    metric: np.ndarray,
    start_idx: int,
    end_idx: int,
    target_idx=None,
):
    start_idx = int(max(0, start_idx))
    end_idx = int(min(len(metric) - 1, end_idx))
    if end_idx <= start_idx:
        return start_idx

    local_candidates = []
    for idx in range(start_idx + 1, end_idx):
        if metric[idx - 1] >= metric[idx] <= metric[idx + 1]:
            local_candidates.append(idx)
    if not local_candidates:
        local_candidates = list(range(start_idx, end_idx + 1))

    if target_idx is None:
        return int(min(local_candidates, key=lambda idx: float(metric[idx])))

    span = max(end_idx - start_idx, 1)
    scale = max(float(np.nanmax(metric[start_idx:end_idx + 1])), 1.0)
    target_idx = float(target_idx)

    def score(idx: int):
        value_score = float(metric[idx]) / scale
        distance_score = 0.12 * abs(float(idx) - target_idx) / span
        return value_score + distance_score

    return int(min(local_candidates, key=score))


def _build_cluster_local_metric(
    x: np.ndarray,
    y_corrected: np.ndarray,
    y_smooth: np.ndarray,
    start_idx: int,
    end_idx: int,
):
    baseline = _estimate_local_linear_baseline(x, y_corrected, start_idx, end_idx)
    if baseline.size == 0:
        return None

    baseline_full = np.zeros(len(x), dtype=float)
    baseline_full[start_idx:end_idx + 1] = baseline
    corrected_local = np.zeros(len(x), dtype=float)
    smooth_local = np.zeros(len(x), dtype=float)
    corrected_local[start_idx:end_idx + 1] = np.clip(
        y_corrected[start_idx:end_idx + 1] - baseline,
        0.0,
        None,
    )
    smooth_local[start_idx:end_idx + 1] = np.clip(
        y_smooth[start_idx:end_idx + 1] - baseline,
        0.0,
        None,
    )
    metric = (
        CLUSTER_BOUNDARY_METRIC_SMOOTH_WEIGHT * smooth_local
        + (1.0 - CLUSTER_BOUNDARY_METRIC_SMOOTH_WEIGHT) * corrected_local
    )
    return baseline_full, corrected_local, smooth_local, metric


def _extract_sharp_isolated_peak_geometry(
    df: pd.DataFrame,
    target_rt: float,
    search_half_window: float = SMALL_PEAK_SHARP_SEARCH_HALF_WINDOW,
    apex_search_radius: float = SMALL_PEAK_SHARP_APEX_SEARCH_RADIUS,
):
    if df is None or df.empty or not np.isfinite(target_rt):
        return None

    x_col = _get_x_column_name(df)
    x = df[x_col].to_numpy(dtype=float)
    y_corrected_raw = df["y_corrected"].to_numpy(dtype=float)
    y_smooth = df["y_smooth"].to_numpy(dtype=float)
    if len(x) < 5:
        return None

    left_x = float(target_rt - search_half_window)
    right_x = float(target_rt + search_half_window)
    left_idx = int(np.searchsorted(x, left_x, side="left"))
    right_idx = int(np.searchsorted(x, right_x, side="right") - 1)
    left_idx = max(0, min(left_idx, len(x) - 2))
    right_idx = max(left_idx + 1, min(right_idx, len(x) - 1))

    local_metric_pack = _build_cluster_local_metric(
        x=x,
        y_corrected=y_corrected_raw,
        y_smooth=y_smooth,
        start_idx=left_idx,
        end_idx=right_idx,
    )
    if local_metric_pack is None:
        return None
    _, _, _, boundary_metric = local_metric_pack
    positive_metric = np.clip(boundary_metric, 0.0, None)

    apex_left = int(np.searchsorted(x, float(target_rt - apex_search_radius), side="left"))
    apex_right = int(np.searchsorted(x, float(target_rt + apex_search_radius), side="right") - 1)
    apex_left = max(left_idx, min(apex_left, right_idx))
    apex_right = max(apex_left, min(apex_right, right_idx))
    candidate_apices = []
    for idx in range(max(apex_left + 1, 1), min(apex_right, len(x) - 2)):
        if positive_metric[idx - 1] <= positive_metric[idx] >= positive_metric[idx + 1]:
            candidate_apices.append(idx)
    if candidate_apices:
        apex_idx = int(min(
            candidate_apices,
            key=lambda idx: (abs(float(x[idx]) - float(target_rt)), -float(positive_metric[idx])),
        ))
    else:
        apex_idx = int(np.argmin(np.abs(x[apex_left:apex_right + 1] - float(target_rt)))) + apex_left
    apex_height = float(positive_metric[apex_idx])
    if apex_height <= 0:
        return None

    local_noise = max(_robust_sigma(y_corrected_raw[left_idx:right_idx + 1]), 1.0)
    threshold = max(apex_height * SMALL_PEAK_SHARP_THRESHOLD_FRACTION, local_noise * SMALL_PEAK_SHARP_THRESHOLD_SIGMA)

    start_idx = apex_idx
    while start_idx > left_idx and positive_metric[start_idx] > threshold:
        start_idx -= 1
    end_idx = apex_idx
    while end_idx < right_idx and positive_metric[end_idx] > threshold:
        end_idx += 1

    if start_idx > left_idx:
        refine_left = slice(max(left_idx, start_idx - 2), min(apex_idx + 1, start_idx + 3))
        start_idx = int(refine_left.start + np.argmin(positive_metric[refine_left]))
    if end_idx < right_idx:
        refine_right = slice(max(apex_idx, end_idx - 2), min(right_idx + 1, end_idx + 3))
        end_idx = int(refine_right.start + np.argmin(positive_metric[refine_right]))

    max_half_width_idx = max(2, int(round(SMALL_PEAK_SHARP_MAX_HALF_WIDTH / max(float(np.median(np.diff(x))), 1e-6))))
    if apex_idx - start_idx > max_half_width_idx:
        target_idx = max(left_idx, apex_idx - max_half_width_idx)
        start_idx = int(target_idx + np.argmin(positive_metric[target_idx:apex_idx + 1]))
    if end_idx - apex_idx > max_half_width_idx:
        target_idx = min(right_idx, apex_idx + max_half_width_idx)
        end_idx = int(apex_idx + np.argmin(positive_metric[apex_idx:target_idx + 1]))

    left_width = float(x[apex_idx] - x[start_idx])
    right_width = float(x[end_idx] - x[apex_idx])
    if left_width > 0 and right_width > 0:
        if left_width > right_width * SMALL_PEAK_SHARP_MAX_ASYMMETRY:
            target_x = float(x[apex_idx] - right_width * SMALL_PEAK_SHARP_MAX_ASYMMETRY)
            target_idx = int(np.argmin(np.abs(x[left_idx:apex_idx + 1] - target_x))) + left_idx
            start_idx = int(target_idx + np.argmin(positive_metric[target_idx:apex_idx + 1]))
        elif right_width > left_width * SMALL_PEAK_SHARP_MAX_ASYMMETRY:
            target_x = float(x[apex_idx] + left_width * SMALL_PEAK_SHARP_MAX_ASYMMETRY)
            target_idx = int(np.argmin(np.abs(x[apex_idx:right_idx + 1] - target_x))) + apex_idx
            end_idx = int(apex_idx + np.argmin(positive_metric[apex_idx:target_idx + 1]))

    if end_idx <= start_idx:
        return None

    area = float(np.trapezoid(np.clip(y_corrected_raw[start_idx:end_idx + 1], 0.0, None), x[start_idx:end_idx + 1]))
    return {
        "start_idx": int(start_idx),
        "apex_idx": int(apex_idx),
        "end_idx": int(end_idx),
        "start_x": float(x[start_idx]),
        "apex_x": float(x[apex_idx]),
        "end_x": float(x[end_idx]),
        "area": area,
    }


def _extract_sharp_peak_geometry_within_bounds(
    df: pd.DataFrame,
    target_rt: float,
    start_x: float,
    end_x: float,
    threshold_fraction: float,
    threshold_sigma: float,
):
    if (
        df is None or df.empty
        or not np.isfinite(target_rt)
        or not np.isfinite(start_x)
        or not np.isfinite(end_x)
        or end_x <= start_x
    ):
        return None

    x_col = _get_x_column_name(df)
    x = df[x_col].to_numpy(dtype=float)
    y_corrected_raw = df["y_corrected"].to_numpy(dtype=float)
    y_smooth = df["y_smooth"].to_numpy(dtype=float)
    if len(x) < 5:
        return None

    left_idx = int(np.searchsorted(x, float(start_x), side="left"))
    right_idx = int(np.searchsorted(x, float(end_x), side="right") - 1)
    left_idx = max(0, min(left_idx, len(x) - 2))
    right_idx = max(left_idx + 1, min(right_idx, len(x) - 1))

    local_metric_pack = _build_cluster_local_metric(
        x=x,
        y_corrected=y_corrected_raw,
        y_smooth=y_smooth,
        start_idx=left_idx,
        end_idx=right_idx,
    )
    if local_metric_pack is None:
        return None
    _, _, _, boundary_metric = local_metric_pack
    positive_metric = np.clip(boundary_metric, 0.0, None)

    candidate_apices = []
    for idx in range(max(left_idx + 1, 1), min(right_idx, len(x) - 2)):
        if positive_metric[idx - 1] <= positive_metric[idx] >= positive_metric[idx + 1]:
            candidate_apices.append(idx)
    if candidate_apices:
        apex_idx = int(min(
            candidate_apices,
            key=lambda idx: (abs(float(x[idx]) - float(target_rt)), -float(positive_metric[idx])),
        ))
    else:
        apex_idx = int(np.argmin(np.abs(x[left_idx:right_idx + 1] - float(target_rt)))) + left_idx

    apex_height = float(positive_metric[apex_idx])
    if apex_height <= 0:
        return None

    local_noise = max(_robust_sigma(y_corrected_raw[left_idx:right_idx + 1]), 1.0)
    threshold = max(apex_height * float(threshold_fraction), local_noise * float(threshold_sigma))

    start_idx = apex_idx
    while start_idx > left_idx and positive_metric[start_idx] > threshold:
        start_idx -= 1
    end_idx = apex_idx
    while end_idx < right_idx and positive_metric[end_idx] > threshold:
        end_idx += 1

    if start_idx > left_idx:
        refine_left = slice(max(left_idx, start_idx - 2), min(apex_idx + 1, start_idx + 3))
        start_idx = int(refine_left.start + np.argmin(positive_metric[refine_left]))
    if end_idx < right_idx:
        refine_right = slice(max(apex_idx, end_idx - 2), min(right_idx + 1, end_idx + 3))
        end_idx = int(refine_right.start + np.argmin(positive_metric[refine_right]))

    if end_idx <= start_idx:
        return None

    area = float(np.trapezoid(np.clip(y_corrected_raw[start_idx:end_idx + 1], 0.0, None), x[start_idx:end_idx + 1]))
    return {
        "start_idx": int(start_idx),
        "apex_idx": int(apex_idx),
        "end_idx": int(end_idx),
        "start_x": float(x[start_idx]),
        "apex_x": float(x[apex_idx]),
        "end_x": float(x[end_idx]),
        "area": area,
    }


def refine_small_peak_integrations(
    df: pd.DataFrame,
    matched_targets_df: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets_df.copy()
    if df is None or df.empty or out is None or out.empty:
        return out

    for code, spec in SMALL_PEAK_SHARP_SPECS.items():
        row_idx_list = out.index[out["code"] == code].tolist()
        if not row_idx_list:
            continue
        row_idx = int(row_idx_list[0])
        found_rt = pd.to_numeric(pd.Series([out.at[row_idx, "found_rt"]]), errors="coerce").iloc[0]
        current_area = pd.to_numeric(pd.Series([out.at[row_idx, "area"]]), errors="coerce").iloc[0]
        current_percent = pd.to_numeric(pd.Series([out.at[row_idx, "percent_area"]]), errors="coerce").iloc[0]
        current_start = pd.to_numeric(pd.Series([out.at[row_idx, "integration_start_x"]]), errors="coerce").iloc[0]
        current_end = pd.to_numeric(pd.Series([out.at[row_idx, "integration_end_x"]]), errors="coerce").iloc[0]
        if not (np.isfinite(found_rt) and np.isfinite(current_area) and np.isfinite(current_start) and np.isfinite(current_end)):
            continue
        if np.isfinite(current_percent) and current_percent > float(spec["max_percent"]):
            continue

        if str(spec.get("mode", "isolated")) == "bounded":
            left_width = float(found_rt - current_start)
            right_width = float(current_end - found_rt)
            current_asymmetry = right_width / max(left_width, 1e-9) if left_width > 0 else np.inf
            if current_asymmetry < float(spec.get("min_asymmetry", 0.0)):
                continue
            geom = _extract_sharp_peak_geometry_within_bounds(
                df=df,
                target_rt=float(found_rt),
                start_x=float(current_start),
                end_x=float(current_end),
                threshold_fraction=float(spec["threshold_fraction"]),
                threshold_sigma=float(spec["threshold_sigma"]),
            )
        else:
            geom = _extract_sharp_isolated_peak_geometry(df, float(found_rt))
        if geom is None:
            continue

        current_width = float(current_end - current_start)
        new_width = float(geom["end_x"] - geom["start_x"])
        if not (np.isfinite(current_width) and current_width > 0 and np.isfinite(new_width) and new_width > 0):
            continue

        area_ratio = float(geom["area"] / current_area) if current_area > 0 else np.nan
        if not np.isfinite(area_ratio):
            continue
        if area_ratio < float(spec["min_area_ratio"]) or area_ratio > float(spec["max_area_ratio"]):
            continue
        if new_width >= current_width * float(spec["max_width_ratio"]):
            continue

        out.at[row_idx, "found_rt"] = float(geom["apex_x"])
        out.at[row_idx, "area"] = float(geom["area"])
        out.at[row_idx, "integration_start_x"] = float(geom["start_x"])
        out.at[row_idx, "integration_end_x"] = float(geom["end_x"])
        out.at[row_idx, "status"] = _append_status_suffix(out.at[row_idx, "status"], "sharp")

    return _recompute_matched_percent_area(out)


def tighten_overwide_c22_cluster_tails(
    df: pd.DataFrame,
    matched_targets_df: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets_df.copy()
    if (
        not ENABLE_C22_TAIL_TIGHTENING
        or df is None or df.empty
        or out is None or out.empty
    ):
        return out

    c22_codes = ["C22:6", "C22:5", "C22:4"]
    cluster = out[out["code"].isin(c22_codes)].copy()
    if len(cluster) != len(c22_codes):
        return out

    cluster["found_rt"] = pd.to_numeric(cluster["found_rt"], errors="coerce")
    cluster["integration_start_x"] = pd.to_numeric(cluster["integration_start_x"], errors="coerce")
    cluster["integration_end_x"] = pd.to_numeric(cluster["integration_end_x"], errors="coerce")
    cluster["area"] = pd.to_numeric(cluster["area"], errors="coerce")
    if cluster[["found_rt", "integration_start_x", "integration_end_x", "area"]].isna().any().any():
        return out

    ordered = cluster.set_index("code").loc[c22_codes].reset_index()
    current_widths = (ordered["integration_end_x"] - ordered["integration_start_x"]).to_numpy(dtype=float)
    if not np.all(np.isfinite(current_widths)):
        return out
    dpa_area = float(ordered.loc[ordered["code"] == "C22:5", "area"].iloc[0])
    c22_4_area = float(ordered.loc[ordered["code"] == "C22:4", "area"].iloc[0])
    dha_width = float(ordered.loc[ordered["code"] == "C22:6", "integration_end_x"].iloc[0] - ordered.loc[ordered["code"] == "C22:6", "integration_start_x"].iloc[0])
    mean_width = float(np.mean(current_widths))
    ratio_trigger = (
        c22_4_area > 0
        and dpa_area / c22_4_area > C22_TAIL_TIGHTENING_DPA_RATIO_TRIGGER
        and dha_width > C22_TAIL_TIGHTENING_DHA_WIDTH_TRIGGER
    )
    if mean_width <= C22_TAIL_TIGHTENING_MEAN_WIDTH and not ratio_trigger:
        return out

    x_col = _get_x_column_name(df)
    x = df[x_col].to_numpy(dtype=float)
    y_corrected_raw = df["y_corrected"].to_numpy(dtype=float)
    y_smooth = df["y_smooth"].to_numpy(dtype=float)
    left_x = float(np.min(ordered["integration_start_x"]) - 0.006)
    right_x = float(np.max(ordered["integration_end_x"]) + 0.006)
    left_idx = int(np.searchsorted(x, left_x, side="left"))
    right_idx = int(np.searchsorted(x, right_x, side="right") - 1)
    left_idx = max(0, min(left_idx, len(x) - 2))
    right_idx = max(left_idx + 1, min(right_idx, len(x) - 1))

    local_metric_pack = _build_cluster_local_metric(
        x=x,
        y_corrected=y_corrected_raw,
        y_smooth=y_smooth,
        start_idx=left_idx,
        end_idx=right_idx,
    )
    if local_metric_pack is None:
        return out
    _, _, _, boundary_metric = local_metric_pack
    y_corrected = np.clip(y_corrected_raw, 0.0, None)

    tightened = []
    previous_end_idx = None
    for _, row in ordered.iterrows():
        current_start_idx = int(np.searchsorted(x, float(row["integration_start_x"]), side="left"))
        current_end_idx = int(np.searchsorted(x, float(row["integration_end_x"]), side="right") - 1)
        apex_idx = int(np.argmin(np.abs(x - float(row["found_rt"]))))
        current_start_idx = max(left_idx, min(current_start_idx, apex_idx - 1))
        current_end_idx = min(right_idx, max(current_end_idx, apex_idx + 1))
        left_half = max(apex_idx - current_start_idx, 1)
        right_half = max(current_end_idx - apex_idx, 1)
        target_start_idx = max(left_idx, apex_idx - int(round(left_half * C22_TAIL_TIGHTENING_WIDTH_SCALE)))
        target_end_idx = min(right_idx, apex_idx + int(round(right_half * C22_TAIL_TIGHTENING_WIDTH_SCALE)))
        new_start_idx = _find_preferred_minimum_index(
            boundary_metric,
            current_start_idx,
            apex_idx,
            target_idx=target_start_idx,
        )
        new_end_idx = _find_preferred_minimum_index(
            boundary_metric,
            apex_idx,
            current_end_idx,
            target_idx=target_end_idx,
        )
        if previous_end_idx is not None and new_start_idx <= previous_end_idx:
            new_start_idx = previous_end_idx + 1
        if new_end_idx <= new_start_idx:
            return out
        previous_end_idx = new_end_idx
        new_area = float(np.trapezoid(y_corrected[new_start_idx:new_end_idx + 1], x[new_start_idx:new_end_idx + 1]))
        tightened.append({
            "code": row["code"],
            "start_idx": int(new_start_idx),
            "end_idx": int(new_end_idx),
            "area": new_area,
        })

    current_cluster_area = float(ordered["area"].sum())
    new_cluster_area = float(sum(item["area"] for item in tightened))
    if current_cluster_area <= 0:
        return out
    area_ratio = new_cluster_area / current_cluster_area
    if not (C22_TAIL_TIGHTENING_AREA_RATIO_MIN <= area_ratio <= C22_TAIL_TIGHTENING_AREA_RATIO_MAX):
        return out

    for item in tightened:
        row_idx = out.index[out["code"] == item["code"]][0]
        out.at[row_idx, "area"] = float(item["area"])
        out.at[row_idx, "integration_start_x"] = float(x[item["start_idx"]])
        out.at[row_idx, "integration_end_x"] = float(x[item["end_idx"]])
        out.at[row_idx, "status"] = _append_status_suffix(out.at[row_idx, "status"], "tailtight")

    return _recompute_matched_percent_area(out)


def refine_overwide_c22_cluster_with_pvfit(
    df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    matched_targets_df: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets_df.copy()
    if (
        not ENABLE_OVERWIDE_C22_PVFIT_REFINEMENT
        or df is None or df.empty
        or peaks_df is None or peaks_df.empty
        or out is None or out.empty
    ):
        return out

    c22_codes = ["C22:6", "C22:5", "C22:4"]
    cluster = out[out["code"].isin(c22_codes)].copy()
    if len(cluster) != len(c22_codes):
        return out

    cluster["integration_start_x"] = pd.to_numeric(cluster["integration_start_x"], errors="coerce")
    cluster["integration_end_x"] = pd.to_numeric(cluster["integration_end_x"], errors="coerce")
    cluster["area"] = pd.to_numeric(cluster["area"], errors="coerce")
    if cluster[["integration_start_x", "integration_end_x", "area"]].isna().any().any():
        return out

    status_text = " ".join(cluster["status"].fillna("").astype(str).tolist())
    if "split" in status_text or "tailtight" not in status_text:
        return out

    ordered = cluster.set_index("code").loc[c22_codes].reset_index()
    widths = (ordered["integration_end_x"] - ordered["integration_start_x"]).to_numpy(dtype=float)
    if not np.all(np.isfinite(widths)):
        return out
    mean_width = float(np.mean(widths))
    dha_width = float(widths[0])
    c22_4_width = float(widths[2])
    if (
        mean_width <= C22_PVFIT_OVERWIDE_MEAN_WIDTH_MIN
        or dha_width <= C22_PVFIT_OVERWIDE_DHA_WIDTH_MIN
        or c22_4_width <= C22_PVFIT_OVERWIDE_C22_4_WIDTH_MIN
    ):
        return out

    previous_total = float(ordered["area"].sum())
    if previous_total <= 0:
        return out

    previous_flag = ENABLE_SPLIT_PSEUDOVOIGT_CLUSTER_FIT
    try:
        globals()["ENABLE_SPLIT_PSEUDOVOIGT_CLUSTER_FIT"] = True
        fit_out, delta_area = _refine_cluster_with_deconvolution(
            df=df,
            peaks_df=peaks_df,
            matched_targets_df=out,
            cluster_codes=c22_codes,
            default_centers=[9.247, 9.280, 9.310],
            window_left=9.22,
            window_right=9.33,
            center_tolerances=[0.010, 0.010, 0.010],
            status="matched_c22_pvfit_tail",
        )
    finally:
        globals()["ENABLE_SPLIT_PSEUDOVOIGT_CLUSTER_FIT"] = previous_flag

    if delta_area == 0.0:
        return out

    fitted_cluster = fit_out[fit_out["code"].isin(c22_codes)].copy()
    fitted_cluster["area"] = pd.to_numeric(fitted_cluster["area"], errors="coerce")
    fitted_total = float(fitted_cluster["area"].fillna(0.0).sum())
    ratio = fitted_total / previous_total if previous_total > 0 else np.nan
    if (
        not np.isfinite(ratio)
        or ratio < C22_PVFIT_AREA_RATIO_MIN
        or ratio > C22_PVFIT_AREA_RATIO_MAX
        or fitted_total >= previous_total
    ):
        return out

    return _recompute_matched_percent_area(fit_out)


def _reintegrate_cluster_by_local_minima(
    df: pd.DataFrame,
    matched_targets_df: pd.DataFrame,
    cluster_codes,
    window_left: float,
    window_right: float,
    status_suffix: str,
    force: bool = False,
) -> pd.DataFrame:
    out = matched_targets_df.copy()
    if not force and not (
        _cluster_has_integration_overlap(out, cluster_codes)
        or _cluster_has_duplicate_peak_ids(out, cluster_codes)
    ):
        return out

    cluster = out[out["code"].isin(cluster_codes)].copy()
    if df is None or df.empty or cluster.empty or len(cluster) != len(cluster_codes):
        return out

    cluster["found_rt"] = pd.to_numeric(cluster["found_rt"], errors="coerce")
    cluster = cluster.dropna(subset=["found_rt"]).sort_values("found_rt").reset_index(drop=False)
    if len(cluster) != len(cluster_codes):
        return out

    x_col = _get_x_column_name(df)
    x = df[x_col].to_numpy(dtype=float)
    y_smooth = df["y_smooth"].to_numpy(dtype=float)
    y_corrected_raw = df["y_corrected"].to_numpy(dtype=float)
    if len(x) < 3:
        return out

    left_idx = int(np.searchsorted(x, float(window_left), side="left"))
    right_idx = int(np.searchsorted(x, float(window_right), side="right") - 1)
    left_idx = max(0, min(left_idx, len(x) - 2))
    right_idx = max(left_idx + 1, min(right_idx, len(x) - 1))

    local_metric_pack = _build_cluster_local_metric(
        x=x,
        y_corrected=y_corrected_raw,
        y_smooth=y_smooth,
        start_idx=left_idx,
        end_idx=right_idx,
    )
    if local_metric_pack is None:
        return out
    _, corrected_local, smooth_local, boundary_metric = local_metric_pack
    y_corrected = np.clip(y_corrected_raw, 0.0, None)
    support_signal = np.maximum(np.clip(corrected_local, 0.0, None), np.clip(smooth_local, 0.0, None))
    cluster_noise = max(_robust_sigma(y_corrected_raw[left_idx:right_idx + 1]), 1.0)

    apex_indices = [int(np.argmin(np.abs(x - float(rt)))) for rt in cluster["found_rt"]]
    if any(apex_idx <= left_idx or apex_idx >= right_idx for apex_idx in apex_indices):
        return out
    if any(apex_indices[i] >= apex_indices[i + 1] for i in range(len(apex_indices) - 1)):
        return out

    left_target = left_idx + 0.25 * max(apex_indices[0] - left_idx, 1)
    left_boundary = _find_preferred_minimum_index(
        boundary_metric,
        left_idx,
        apex_indices[0],
        target_idx=left_target,
    )
    right_target = apex_indices[-1] + 0.75 * max(right_idx - apex_indices[-1], 1)
    right_boundary = _find_preferred_minimum_index(
        boundary_metric,
        apex_indices[-1],
        right_idx,
        target_idx=right_target,
    )
    left_support_threshold = max(
        cluster_noise * PEAK_SUPPORT_THRESHOLD_SIGMA,
        float(max(support_signal[apex_indices[0]], 0.0)) * PEAK_SUPPORT_THRESHOLD_FRACTION,
    )
    right_support_threshold = max(
        cluster_noise * PEAK_SUPPORT_THRESHOLD_SIGMA,
        float(max(support_signal[apex_indices[-1]], 0.0)) * PEAK_SUPPORT_THRESHOLD_FRACTION,
    )
    left_boundary = _extend_boundary_to_support(
        signal=support_signal,
        start_idx=int(left_boundary),
        limit_idx=int(left_idx),
        direction=-1,
        threshold=float(left_support_threshold),
        consecutive_points=PEAK_SUPPORT_CONSECUTIVE_POINTS,
    )
    right_boundary = _extend_boundary_to_support(
        signal=support_signal,
        start_idx=int(right_boundary),
        limit_idx=int(right_idx),
        direction=1,
        threshold=float(right_support_threshold),
        consecutive_points=PEAK_SUPPORT_CONSECUTIVE_POINTS,
    )
    boundaries = [left_boundary]
    for i in range(len(apex_indices) - 1):
        split_start = apex_indices[i]
        split_end = apex_indices[i + 1]
        if split_end <= split_start:
            return out
        split_idx = _find_preferred_minimum_index(
            boundary_metric,
            split_start,
            split_end,
            target_idx=0.5 * (split_start + split_end),
        )
        if split_idx <= boundaries[-1]:
            split_idx = max(boundaries[-1] + 1, int(round(0.5 * (apex_indices[i] + apex_indices[i + 1]))))
        boundaries.append(split_idx)
    boundaries.append(right_boundary)

    if any(boundaries[i] >= boundaries[i + 1] for i in range(len(boundaries) - 1)):
        return out

    for cluster_pos, (_, row) in enumerate(cluster.iterrows()):
        row_idx = int(row["index"])
        start_idx = int(boundaries[cluster_pos])
        end_idx = int(boundaries[cluster_pos + 1])
        if end_idx <= start_idx:
            continue
        area = float(np.trapezoid(y_corrected[start_idx:end_idx + 1], x[start_idx:end_idx + 1]))
        out.at[row_idx, "area"] = area
        out.at[row_idx, "integration_start_x"] = float(x[start_idx])
        out.at[row_idx, "integration_end_x"] = float(x[end_idx])
        out.at[row_idx, "status"] = _append_status_suffix(out.at[row_idx, "status"], status_suffix)

    return out


def _should_force_c18_valley_split(matched_targets_df: pd.DataFrame) -> bool:
    cluster = matched_targets_df[matched_targets_df["code"].isin(["C18:2N6C", "C18:1N9C", "C18:3N3", "C18:0"])].copy()
    if len(cluster) != 4:
        return False

    cluster["area"] = pd.to_numeric(cluster["area"], errors="coerce")
    if cluster["area"].isna().any():
        return False

    area_by_code = cluster.set_index("code")["area"]
    c18_2 = float(area_by_code.get("C18:2N6C", np.nan))
    c18_1 = float(area_by_code.get("C18:1N9C", np.nan))
    c18_3 = float(area_by_code.get("C18:3N3", np.nan))
    if not (np.isfinite(c18_2) and np.isfinite(c18_1) and np.isfinite(c18_3)):
        return False

    return bool(
        c18_2 > max(c18_1 * 2.2, 4500.0)
        and c18_3 < max(350.0, c18_1 * 0.18)
    )


def refine_cluster_areas_by_local_valleys(
    df: pd.DataFrame,
    matched_targets_df: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets_df.copy()
    if df is None or df.empty or out is None or out.empty:
        return out

    out = _reintegrate_cluster_by_local_minima(
        df=df,
        matched_targets_df=out,
        cluster_codes=["C18:2N6C", "C18:1N9C", "C18:3N3", "C18:0"],
        window_left=7.56,
        window_right=7.79,
        status_suffix="valley",
        force=_should_force_c18_valley_split(out),
    )
    out = _reintegrate_cluster_by_local_minima(
        df=df,
        matched_targets_df=out,
        cluster_codes=["C20:4N6", "C20:5", "C20:3N8"],
        window_left=8.34,
        window_right=8.50,
        status_suffix="valley",
    )
    out = _reintegrate_cluster_by_local_minima(
        df=df,
        matched_targets_df=out,
        cluster_codes=["C22:6", "C22:5", "C22:4"],
        window_left=9.22,
        window_right=9.34,
        status_suffix="valley",
    )
    return _recompute_matched_percent_area(out)


def refine_c18_c20_cluster_matches(
    df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    matched_targets_df: pd.DataFrame,
):
    out = matched_targets_df.copy()
    peaks_out = peaks_df.copy()
    if df is None or df.empty or out is None or out.empty:
        return peaks_out, out

    c18_codes = ["C18:2N6C", "C18:1N9C", "C18:3N3", "C18:0"]
    c18_choice = None
    if _cluster_has_integration_overlap(out, c18_codes):
        c18_candidates = _collect_local_cluster_peak_geometries(
            df,
            window_left=7.56,
            window_right=7.79,
            min_prominence=100.0,
            min_area=10.0,
        )
        c18_choice = _select_ordered_cluster_peaks(
            c18_candidates,
            target_apexes=[7.593, 7.623, 7.650, 7.750],
            max_distances=[0.022, 0.022, 0.025, 0.028],
        )
    if c18_choice is not None:
            peaks_lookup = peaks_out.copy()
            for code, (_, geom) in zip(c18_codes, c18_choice.iterrows()):
                row_idx = out.index[out["code"] == code]
                if len(row_idx):
                    row_idx_int = int(row_idx[0])
                    _assign_local_geometry_bounds_to_row(out, row_idx_int, geom, "matched_c18_local_bounds")
                    peak_match = peaks_lookup[(peaks_lookup["apex_x"] - float(geom["apex_x"])).abs() <= 0.006]
                    if not peak_match.empty:
                        out.at[row_idx_int, "matched_peak_id"] = int(peak_match.sort_values("area", ascending=False).iloc[0]["peak_id"])

    c20_codes = ["C20:4N6", "C20:5", "C20:3N8"]
    c20_candidates = _collect_local_cluster_peak_geometries(
        df,
        window_left=8.34,
        window_right=8.50,
        min_prominence=10.0,
        min_area=2.0,
    )
    c20_choice = _select_ordered_cluster_peaks(
        c20_candidates,
        target_apexes=[8.381, 8.410, 8.467],
        max_distances=[0.025, 0.022, 0.025],
    )
    if c20_choice is not None:
        current_cluster = out[out["code"].isin(c20_codes)].copy()
        current_cluster["area"] = pd.to_numeric(current_cluster["area"], errors="coerce")
        current_cluster["integration_start_x"] = pd.to_numeric(current_cluster["integration_start_x"], errors="coerce")
        current_cluster["integration_end_x"] = pd.to_numeric(current_cluster["integration_end_x"], errors="coerce")

        current_area = float(current_cluster["area"].fillna(0.0).sum())
        local_area = float(c20_choice["area"].sum())
        has_missing = current_cluster["area"].isna().any()
        has_duplicate = _cluster_has_duplicate_peak_ids(out, c20_codes)

        extends_boundaries = False
        ordered_current = current_cluster.set_index("code").reindex(c20_codes)
        for code, (_, geom) in zip(c20_codes, c20_choice.iterrows()):
            current_row = ordered_current.loc[code]
            current_start = current_row.get("integration_start_x")
            current_end = current_row.get("integration_end_x")
            if np.isfinite(current_start) and float(geom["start_x"]) < float(current_start) - C20_LOCAL_BOUNDARY_EXTENSION:
                extends_boundaries = True
            if np.isfinite(current_end) and float(geom["end_x"]) > float(current_end) + C20_LOCAL_BOUNDARY_EXTENSION:
                extends_boundaries = True

        if has_missing or has_duplicate or (
            current_area > 0
            and local_area >= current_area * C20_LOCAL_AREA_RATIO_TRIGGER
            and extends_boundaries
        ):
            peaks_out = _attach_local_peak_records(peaks_out, c20_choice)
            peaks_lookup = peaks_out.copy()
            for code, (_, geom) in zip(c20_codes, c20_choice.iterrows()):
                row_idx = out.index[out["code"] == code]
                if not len(row_idx):
                    continue
                _assign_local_geometry_to_row(out, int(row_idx[0]), geom, "matched_c20_local")
                peak_match = peaks_lookup[(peaks_lookup["apex_x"] - float(geom["apex_x"])).abs() <= 0.006]
                if not peak_match.empty:
                    out.at[int(row_idx[0]), "matched_peak_id"] = int(peak_match.sort_values("area", ascending=False).iloc[0]["peak_id"])

    out = _recompute_matched_percent_area(out)
    return peaks_out, out


def augment_targeted_cluster_peaks(df: pd.DataFrame, peaks_df: pd.DataFrame) -> pd.DataFrame:
    if df is None or df.empty:
        return peaks_df

    noise = max(_robust_sigma(df["y_corrected"].to_numpy(dtype=float)), 1.0)
    specs = [
        {"target_x": 7.623, "search_radius": 0.018, "min_prominence": max(8.0 * noise, 1200.0), "min_area": 120.0},
        {"target_x": 7.650, "search_radius": 0.020, "min_prominence": max(6.0 * noise, 900.0), "min_area": 120.0},
        {"target_x": 7.750, "search_radius": 0.020, "min_prominence": max(8.0 * noise, 1200.0), "min_area": 120.0},
        {"target_x": 8.381, "search_radius": 0.018, "min_prominence": max(8.0 * noise, 1200.0), "min_area": 120.0},
        {"target_x": 8.410, "search_radius": 0.018, "min_prominence": max(5.0 * noise, 700.0), "min_area": 80.0},
        {"target_x": 8.467, "search_radius": 0.020, "min_prominence": max(5.0 * noise, 700.0), "min_area": 80.0},
        {"target_x": 9.252, "search_radius": 0.018, "min_prominence": max(2.0 * noise, 250.0), "min_area": 20.0},
        {"target_x": 9.285, "search_radius": 0.018, "min_prominence": max(2.0 * noise, 220.0), "min_area": 20.0},
        {"target_x": 9.316, "search_radius": 0.018, "min_prominence": max(2.0 * noise, 220.0), "min_area": 20.0},
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
    global _PYOPENMS_PEAK_PICKER
    if _PYOPENMS_PEAK_PICKER is not None:
        return _PYOPENMS_PEAK_PICKER
    picker = oms.PeakPickerChromatogram()
    params = picker.getParameters()
    params.setValue(b"gauss_width", float(PYOPENMS_GAUSS_WIDTH_SECONDS))
    params.setValue(b"signal_to_noise", float(PYOPENMS_SIGNAL_TO_NOISE))
    params.setValue(b"sn_win_len", float(PYOPENMS_SN_WIN_LEN))
    params.setValue(b"use_gauss", b"true")
    params.setValue(b"remove_overlapping_peaks", b"false")
    picker.setParameters(params)
    _PYOPENMS_PEAK_PICKER = picker
    return picker


def detect_peaks_with_pyopenms(df: pd.DataFrame) -> pd.DataFrame:
    if not ENABLE_PYOPENMS_PEAK_ASSIST or oms is None or df is None or df.empty:
        return pd.DataFrame(columns=PEAK_RECORD_COLUMNS)

    try:
        x_col = _get_x_column_name(df)
        x = df[x_col].to_numpy(dtype=float)
        y_corrected = np.clip(df["y_corrected"].to_numpy(dtype=float), 0.0, None)
        if x.size < 20 or not np.any(y_corrected > 0):
            return pd.DataFrame(columns=PEAK_RECORD_COLUMNS)

        chromatogram = oms.MSChromatogram()
        chromatogram.set_peaks((x * 60.0, y_corrected))
        picked = oms.MSChromatogram()
        picker = _get_pyopenms_peak_picker()
        picker.pickChromatogram(chromatogram, picked)
        peak_rts, _ = picked.get_peaks()
    except Exception:
        return pd.DataFrame(columns=PEAK_RECORD_COLUMNS)

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


def _fit_chebyshev_baseline(x: np.ndarray, y: np.ndarray, degree: int, n_bins: int, lower_quantile: float, n_iter: int, sigma_threshold: float):
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


def add_baseline_to_dataframe(
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
        raise ValueError("Недостаточно точек для baseline correction.")

    resolved_degree = 6 if degree is None else int(degree)
    baseline, coeffs = _fit_chebyshev_baseline(x, y, resolved_degree, n_bins, lower_quantile, n_iter, sigma_threshold)
    out["baseline"] = baseline
    out["y_corrected"] = y - baseline

    if WRITE_CHEBYSHEV_COEFFICIENTS:
        coeff_path = get_runtime_app_dir() / "chebyshev_coefficients.csv"
        coeff_df = pd.DataFrame({
            "coefficient_index": np.arange(len(coeffs), dtype=int),
            "coefficient": coeffs,
        })
        coeff_df.to_csv(coeff_path, index=False)
    return out


def add_arpls_baseline_to_dataframe(
    df: pd.DataFrame,
    lam: float = ARPLS_BASELINE_LAM,
) -> pd.DataFrame:
    if Baseline is None:
        return add_baseline_to_dataframe(df, **BASELINE_KWARGS)

    out = df.copy()
    y = out["y"].to_numpy(dtype=float)
    if y.size < 8:
        raise ValueError("Недостаточно точек для baseline correction.")

    baseline, _ = Baseline().arpls(y, lam=float(lam))
    baseline = np.asarray(baseline, dtype=float)
    out["baseline"] = baseline
    out["y_corrected"] = y - baseline
    return out


def add_savgol_and_derivatives_to_dataframe(df: pd.DataFrame, polyorder: int = 3, candidate_windows=None):
    out = df.copy()
    y = out["y_corrected"].to_numpy(dtype=float)
    x = out["x_corrected"].to_numpy(dtype=float)
    if candidate_windows is None:
        candidate_windows = [11, 15, 21, 31, 41, 51, 61, 81, 101, 151]

    if y.size <= polyorder + 2:
        raise ValueError("Недостаточно точек для Savitzky-Golay smoothing.")

    valid_windows = sorted({
        int(w) for w in candidate_windows
        if int(w) % 2 == 1 and int(w) > polyorder and int(w) <= (len(y) if len(y) % 2 == 1 else len(y) - 1)
    })
    if SAVGOL_MAX_SELECTED_WINDOW is not None:
        capped_windows = [int(w) for w in valid_windows if int(w) <= int(SAVGOL_MAX_SELECTED_WINDOW)]
        if capped_windows:
            valid_windows = capped_windows
    if not valid_windows:
        fallback = len(y) if len(y) % 2 == 1 else len(y) - 1
        fallback = max(polyorder + 2 + ((polyorder + 2) % 2 == 0), min(fallback, 11))
        valid_windows = [fallback]

    best_window = valid_windows[0]
    best_score = math.inf
    best_smooth = None
    raw_scale = max(_robust_sigma(y), 1e-9)
    y_p99 = float(np.percentile(y, 99))
    y_p99_abs = max(abs(y_p99), 1e-9)

    for window in valid_windows:
        smooth = savgol_filter(y, window_length=window, polyorder=polyorder, mode="interp")
        residual = y - smooth
        noise_score = _robust_sigma(residual) / raw_scale
        curvature_score = _robust_sigma(np.diff(smooth, n=2))
        peak_loss = abs(np.percentile(smooth, 99) - y_p99) / y_p99_abs
        score = noise_score + 0.15 * curvature_score + 2.0 * peak_loss
        if score < best_score:
            best_score = float(score)
            best_window = int(window)
            best_smooth = smooth

    if best_smooth is None:
        best_smooth = savgol_filter(y, window_length=best_window, polyorder=polyorder, mode="interp")

    dy = np.gradient(best_smooth, x)
    out["y_smooth"] = best_smooth
    out["dy"] = dy
    out["d2y"] = np.gradient(dy, x)
    return out, best_window


def detect_peaks_from_derivatives(
    df: pd.DataFrame,
    best_window=None,
    height_sigma: float = 1.5,
    prominence_sigma: float = 2.0,
    rel_height: float = PEAK_INTEGRATION_REL_HEIGHT,
) -> pd.DataFrame:
    if df is None or df.empty:
        return pd.DataFrame(columns=["peak_id", "start_x", "apex_x", "end_x", "area", "percent_area"])

    x_col = _get_x_column_name(df)
    x = df[x_col].to_numpy(dtype=float)
    y_corrected = df["y_corrected"].to_numpy(dtype=float)
    y_smooth = df["y_smooth"].to_numpy(dtype=float)
    y_smooth_positive = np.clip(y_smooth, 0.0, None)
    dx = float(np.median(np.diff(x))) if len(x) > 1 else 1.0
    noise = max(_robust_sigma(y_corrected), 1e-9)

    height_floor = max(np.median(y_smooth) + height_sigma * noise, np.quantile(y_smooth, 0.60))
    prominence_floor = max(prominence_sigma * noise, np.quantile(y_smooth_positive, 0.75) * 0.05)
    min_distance = max(1, int(round(0.03 / max(dx, 1e-9))))
    min_width = max(1, int(round(0.009 / max(dx, 1e-9))))

    peaks, props = find_peaks(
        y_smooth,
        height=height_floor,
        prominence=prominence_floor,
        distance=min_distance,
        width=min_width,
    )
    if peaks.size == 0:
        return pd.DataFrame(columns=["peak_id", "start_x", "apex_x", "end_x", "area", "percent_area"])

    widths = peak_widths(y_smooth, peaks, rel_height=rel_height)
    left_ips = widths[2]
    right_ips = widths[3]

    records = []
    for order, peak_idx in enumerate(peaks, start=1):
        start_idx = max(0, int(np.floor(left_ips[order - 1])))
        end_idx = min(len(x) - 1, int(np.ceil(right_ips[order - 1])))
        if end_idx <= start_idx:
            continue
        x_seg = x[start_idx:end_idx + 1]
        y_seg = np.clip(y_corrected[start_idx:end_idx + 1], 0.0, None)
        area = float(np.trapezoid(y_seg, x_seg))
        records.append({
            "peak_id": order,
            "start_idx": start_idx,
            "apex_idx": int(peak_idx),
            "end_idx": end_idx,
            "start_x": float(x[start_idx]),
            "apex_x": float(x[peak_idx]),
            "end_x": float(x[end_idx]),
            "height": float(props["peak_heights"][order - 1]),
            "prominence": float(props["prominences"][order - 1]),
            "width_points": float(props["widths"][order - 1]),
            "area": area,
        })

    peaks_df = pd.DataFrame(records).sort_values("apex_x").reset_index(drop=True)
    if peaks_df.empty:
        return pd.DataFrame(columns=["peak_id", "start_x", "apex_x", "end_x", "area", "percent_area"])

    peaks_df["peak_id"] = np.arange(1, len(peaks_df) + 1)
    total_area = float(peaks_df["area"].sum())
    peaks_df["percent_area"] = 100.0 * peaks_df["area"] / total_area if total_area > 0 else np.nan
    peaks_df = augment_targeted_cluster_peaks(df, peaks_df)
    pyopenms_peaks_df = detect_peaks_with_pyopenms(df)
    if not pyopenms_peaks_df.empty:
        peaks_df = _merge_peak_records(peaks_df, pyopenms_peaks_df.to_dict("records"))
    return peaks_df


def _pseudo_voigt_unit_area(x: np.ndarray, center: float, fwhm: float, eta: float) -> np.ndarray:
    width = max(float(fwhm), 1e-6)
    mixing = float(np.clip(eta, 0.0, 1.0))
    dx = np.asarray(x, dtype=float) - float(center)
    scaled = dx / width
    gaussian = math.sqrt(4.0 * math.log(2.0) / math.pi) / width * np.exp(-4.0 * math.log(2.0) * scaled * scaled)
    lorentzian = (2.0 / (math.pi * width)) / (1.0 + 4.0 * scaled * scaled)
    return mixing * lorentzian + (1.0 - mixing) * gaussian


def _derive_pseudo_voigt_boundaries(components, x_left: float, x_right: float):
    if not components:
        return []

    ordered = sorted(components, key=lambda item: float(item["center"]))
    crossings = []
    for left, right in zip(ordered[:-1], ordered[1:]):
        dense_x = np.linspace(float(left["center"]), float(right["center"]), 360)
        left_curve = float(left["area"]) * _pseudo_voigt_unit_area(dense_x, left["center"], left["fwhm"], left["eta"])
        right_curve = float(right["area"]) * _pseudo_voigt_unit_area(dense_x, right["center"], right["fwhm"], right["eta"])
        diff = left_curve - right_curve
        crossing = None
        for idx in range(1, len(diff)):
            if diff[idx - 1] == 0.0:
                crossing = float(dense_x[idx - 1])
                break
            if diff[idx] == 0.0 or np.sign(diff[idx]) != np.sign(diff[idx - 1]):
                crossing = float(dense_x[idx])
                break
        if crossing is None:
            crossing = 0.5 * (float(left["center"]) + float(right["center"]))
        crossings.append(crossing)

    resolved = []
    for idx, component in enumerate(ordered):
        center = float(component["center"])
        fwhm = float(component["fwhm"])
        soft_left = max(float(x_left), center - 2.6 * fwhm)
        soft_right = min(float(x_right), center + 2.6 * fwhm)
        start_x = soft_left if idx == 0 else max(soft_left, float(crossings[idx - 1]))
        end_x = soft_right if idx == len(ordered) - 1 else min(soft_right, float(crossings[idx]))
        if end_x <= start_x:
            half_width = max(0.5 * fwhm, 1e-4)
            start_x = max(float(x_left), center - half_width)
            end_x = min(float(x_right), center + half_width)
        resolved.append((float(start_x), float(end_x)))
    return resolved


def extract_sample_name_from_header(file_path: Path) -> str:
    pattern = re.compile(r"\bO\d+_[A-Za-z0-9._-]+\b")
    try:
        with file_path.open("r", encoding="utf-8", errors="ignore") as f:
            header_lines = [next(f, "") for _ in range(3)]
        match = pattern.search(" ".join(header_lines))
        if match:
            return match.group(0)
    except Exception:
        pass
    return file_path.stem


def extract_sample_name_from_text(text: str, fallback: str) -> str:
    pattern = re.compile(r"\bO\d+_[A-Za-z0-9._-]+\b")
    match = pattern.search(text or "")
    return match.group(0) if match else fallback


def omega_sample_sort_key(name: str):
    text = name or ""
    match = re.search(r"\bO(\d+)\b", text)
    if match:
        return (0, int(match.group(1)), text)
    match = re.search(r"\bO(\d+)", text)
    if match:
        return (0, int(match.group(1)), text)
    return (1, text)


def finalize_chromatogram_dataframe(df: pd.DataFrame, cutoff_minutes: float = 4.0) -> pd.DataFrame:
    df = df.copy()
    df["x"] = pd.to_numeric(df["x"], errors="coerce")
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    df = df.dropna(subset=["x", "y"]).reset_index(drop=True)

    unique_x = np.sort(df["x"].unique())
    if len(unique_x) < 2:
        raise ValueError("Недостаточно уникальных x для определения шага.")

    step = np.median(np.diff(unique_x))
    corrected_x = []
    i = 0
    n = len(df)
    while i < n:
        current_x = df.loc[i, "x"]
        j = i + 1
        while j < n and df.loc[j, "x"] == current_x:
            j += 1
        group_size = j - i
        offsets = np.linspace(0, step * (group_size - 1) / max(group_size, 1), group_size)
        corrected_x.extend(current_x + offsets)
        i = j

    df["x_corrected"] = corrected_x
    return df[df["x_corrected"] >= cutoff_minutes].reset_index(drop=True)


def load_chromatogram_csv(file_path: Path, cutoff_minutes: float = 4.0) -> pd.DataFrame:
    df = pd.read_csv(file_path, skiprows=3, header=None, names=["x", "y"])
    return finalize_chromatogram_dataframe(df, cutoff_minutes=cutoff_minutes)


def _is_chromtab_file(file_path: Path) -> bool:
    try:
        with file_path.open("r", encoding="utf-8", errors="ignore") as f:
            first_line = next(f, "").strip()
            second_line = next(f, "").strip()
            third_line = next(f, "").strip()
    except Exception:
        return False
    return (
        first_line.startswith('"Path","File","Date Acquired"')
        and third_line.startswith('"Signal: ')
    )


def load_chromtab_batches(file_path: Path, cutoff_minutes: float = 4.0):
    batches = []
    current_meta = None
    current_rows = []

    def flush_current():
        nonlocal current_meta, current_rows
        if current_meta is None or not current_rows:
            current_rows = []
            return
        raw_df = pd.DataFrame(current_rows, columns=["x", "y"])
        batch_df = finalize_chromatogram_dataframe(raw_df, cutoff_minutes=cutoff_minutes)
        file_name = current_meta.get("file_name", "")
        sample_name = extract_sample_name_from_text(
            " ".join(filter(None, [current_meta.get("signal_name", ""), file_name])),
            fallback=Path(file_name).stem if file_name else f"batch_{len(batches)+1}",
        )
        batches.append({
            "file_name": file_name,
            "signal_name": current_meta.get("signal_name", ""),
            "acquired_at": current_meta.get("acquired_at", ""),
            "source_path": current_meta.get("source_path", ""),
            "sample_name": sample_name,
            "dataframe": batch_df,
        })
        current_meta = None
        current_rows = []

    with file_path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith('"Path","File","Date Acquired"'):
                flush_current()
                current_meta = {"source_path": "", "file_name": "", "acquired_at": "", "signal_name": ""}
                continue

            if current_meta is not None and not current_meta["file_name"] and line.startswith('"'):
                parsed = next(csv.reader([line]))
                current_meta["source_path"] = parsed[0] if len(parsed) > 0 else ""
                current_meta["file_name"] = parsed[1] if len(parsed) > 1 else ""
                current_meta["acquired_at"] = parsed[2] if len(parsed) > 2 else ""
                continue

            if current_meta is not None and line.startswith('"Signal: '):
                current_meta["signal_name"] = line.strip('"')
                continue

            if current_meta is not None:
                parsed = next(csv.reader([line]), [])
                if len(parsed) < 2:
                    continue
                try:
                    current_rows.append((float(parsed[0]), float(parsed[1])))
                except ValueError:
                    continue

    flush_current()
    if not batches:
        raise ValueError("Не удалось разобрать серии из CHROMTAB.CSV.")
    return sorted(batches, key=lambda batch: omega_sample_sort_key(batch.get("sample_name", batch.get("file_name", ""))))


def compute_clean_omega_metrics(matched_targets_df: pd.DataFrame) -> dict:
    result = {
        "omega3_trio": np.nan,
        "omega3_trio_strict": np.nan,
        "omega3_trio_corrected": np.nan,
        "total_area": np.nan,
        "effective_total_area": np.nan,
        "epa_area": 0.0,
        "dha_area": 0.0,
        "dpa_area": 0.0,
        "epa_effective_area": 0.0,
        "epa_neighbor_area": 0.0,
        "epa_overlap_credit_area": 0.0,
        "epa_overlap_fraction": 0.0,
        "epa_overlap_model_applied": False,
        "epa_overlap_extra_scale": 1.0,
        "c22_overlap_source_area": 0.0,
        "c22_overlap_credit_area": 0.0,
        "c22_overlap_fraction": 0.0,
        "c22_overlap_legacy_fraction": 0.0,
        "c22_overlap_model_fraction": np.nan,
        "c22_overlap_model_applied": False,
        "c22_reference_ratio": np.nan,
        "c22_width_scale": 1.0,
        "c18_denominator_scale": 1.0,
    }
    if matched_targets_df is None or matched_targets_df.empty:
        return result
    valid = matched_targets_df.copy()
    valid["area"] = pd.to_numeric(valid.get("area"), errors="coerce")
    valid = valid.dropna(subset=["area"])
    if valid.empty:
        return result
    total_area = float(valid["area"].sum())
    if total_area <= 0 or not np.isfinite(total_area):
        return result

    def area_of(code: str) -> float:
        row = valid[valid["code"] == code]
        return float(row["area"].iloc[0]) if not row.empty else 0.0

    epa = area_of("C20:5")
    dha = area_of("C22:6")
    dpa = area_of("C22:5")
    omega_value = 100.0 * (epa + dha + dpa) / total_area
    result.update({
        "omega3_trio": omega_value,
        "omega3_trio_strict": omega_value,
        "omega3_trio_corrected": omega_value,
        "total_area": total_area,
        "effective_total_area": total_area,
        "epa_area": epa,
        "dha_area": dha,
        "dpa_area": dpa,
        "epa_effective_area": epa,
    })
    return result


def _add_derivatives_for_gui(processed_df: pd.DataFrame) -> pd.DataFrame:
    out = processed_df.copy()
    x_col = _get_x_column_name(out)
    x = out[x_col].to_numpy(dtype=float)
    y_smooth = out["y_smooth"].to_numpy(dtype=float)
    if "dy" not in out.columns:
        out["dy"] = np.gradient(y_smooth, x) if len(x) > 2 else 0.0
    if "ddy" not in out.columns:
        dy = out["dy"].to_numpy(dtype=float)
        out["ddy"] = np.gradient(dy, x) if len(x) > 2 else 0.0
    return out


def _integration_row_to_peak_record(processed_df: pd.DataFrame, row: pd.Series) -> dict | None:
    x_col = _get_x_column_name(processed_df)
    x = processed_df[x_col].to_numpy(dtype=float)
    y = np.clip(processed_df["y_corrected"].to_numpy(dtype=float), 0.0, None)
    y_smooth = processed_df["y_smooth"].to_numpy(dtype=float)
    start_x = pd.to_numeric(pd.Series([row.get("integration_start_x")]), errors="coerce").iloc[0]
    end_x = pd.to_numeric(pd.Series([row.get("integration_end_x")]), errors="coerce").iloc[0]
    apex_x = pd.to_numeric(pd.Series([row.get("found_rt")]), errors="coerce").iloc[0]
    if not all(np.isfinite(value) for value in [start_x, end_x, apex_x]):
        return None
    start_idx = int(np.argmin(np.abs(x - float(start_x))))
    end_idx = int(np.argmin(np.abs(x - float(end_x))))
    apex_idx = int(np.argmin(np.abs(x - float(apex_x))))
    if end_idx <= start_idx:
        return None
    area = float(np.trapezoid(y[start_idx:end_idx + 1], x[start_idx:end_idx + 1]))
    prominence = float(max(0.0, y_smooth[apex_idx] - max(y_smooth[start_idx], y_smooth[end_idx])))
    return {
        "start_idx": start_idx,
        "apex_idx": apex_idx,
        "end_idx": end_idx,
        "start_x": float(x[start_idx]),
        "apex_x": float(x[apex_idx]),
        "end_x": float(x[end_idx]),
        "height": float(y_smooth[apex_idx]),
        "prominence": prominence,
        "width_points": float(end_idx - start_idx),
        "area": max(area, 0.0),
    }


def _build_clean_peaks_df(processed_df: pd.DataFrame, clean_result: dict) -> pd.DataFrame:
    records = []
    used_apex_indices: set[int] = set()
    for peak_idx in np.asarray(clean_result.get("peaks", []), dtype=int):
        try:
            geometry = _extract_peak_geometry(processed_df, int(peak_idx))
        except Exception:
            geometry = None
        if geometry is None:
            continue
        records.append(geometry)
        used_apex_indices.add(int(geometry["apex_idx"]))

    matched = clean_result.get("matched_targets_df")
    if isinstance(matched, pd.DataFrame) and not matched.empty:
        for _, row in matched.iterrows():
            synthetic = _integration_row_to_peak_record(processed_df, row)
            if synthetic is None:
                continue
            apex_idx = int(synthetic["apex_idx"])
            if any(abs(apex_idx - used_idx) <= 2 for used_idx in used_apex_indices):
                continue
            records.append(synthetic)
            used_apex_indices.add(apex_idx)

    if not records:
        return pd.DataFrame(columns=PEAK_RECORD_COLUMNS)
    peaks_df = pd.DataFrame(records).sort_values("apex_x").reset_index(drop=True)
    peaks_df["peak_id"] = np.arange(1, len(peaks_df) + 1)
    total_area = float(pd.to_numeric(peaks_df["area"], errors="coerce").fillna(0.0).sum())
    peaks_df["percent_area"] = 100.0 * peaks_df["area"] / total_area if total_area > 0 else np.nan
    return peaks_df.reindex(columns=PEAK_RECORD_COLUMNS)


def _attach_clean_matched_peak_ids(matched_targets_df: pd.DataFrame, peaks_df: pd.DataFrame) -> pd.DataFrame:
    out = matched_targets_df.copy()
    if "matched_peak_id" not in out.columns:
        out["matched_peak_id"] = np.nan
    if "match_score" not in out.columns:
        out["match_score"] = np.nan
    if "corrected_target_rt" not in out.columns:
        out["corrected_target_rt"] = out.get("target_rt", out.get("expected_rt", np.nan))
    if peaks_df is None or peaks_df.empty:
        return out
    peak_rt = peaks_df["apex_x"].to_numpy(dtype=float)
    peak_ids = peaks_df["peak_id"].to_numpy(dtype=int)
    used: set[int] = set()
    for row_idx, row in out.iterrows():
        found_rt = pd.to_numeric(pd.Series([row.get("found_rt")]), errors="coerce").iloc[0]
        if not np.isfinite(found_rt):
            continue
        distances = np.abs(peak_rt - float(found_rt))
        order = np.argsort(distances)
        chosen_pos = None
        for pos in order:
            if int(peak_ids[pos]) not in used:
                chosen_pos = int(pos)
                break
        if chosen_pos is None:
            chosen_pos = int(order[0])
        if float(distances[chosen_pos]) <= 0.030:
            out.at[row_idx, "matched_peak_id"] = int(peak_ids[chosen_pos])
            out.at[row_idx, "match_score"] = float(distances[chosen_pos])
            used.add(int(peak_ids[chosen_pos]))
    return out


def process_chromatogram_batch(dataframe: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    # The GUI must use the same engine that is exercised by omega_regression.py.
    # The former clean-only route bypassed the validated C18/C20/C22 matching,
    # cluster deconvolution, boundary judge and metric safeguards.  As a result a
    # visually plausible chromatogram could still assign a shoulder to the wrong
    # fatty acid and produce a large field-batch error.
    result = dict(omega_core.process_batch(dataframe, reference_targets))
    result["engine"] = "omega_core"
    result.setdefault("total_area", result.get("omega", {}).get("total_area", np.nan))
    return result


def _compute_cluster_quality_score(matched_targets_df: pd.DataFrame) -> float:
    if matched_targets_df is None or matched_targets_df.empty:
        return -np.inf

    cluster_groups = [
        ["C18:2N6C", "C18:1N9C", "C18:3N3", "C18:0"],
        ["C20:4N6", "C20:5", "C20:3N8"],
        ["C22:6", "C22:5", "C22:4"],
    ]
    score = 0.0
    for cluster_codes in cluster_groups:
        cluster = matched_targets_df.loc[matched_targets_df["code"].isin(cluster_codes), ["found_rt", "area", "matched_peak_id"]]
        found_rt = pd.to_numeric(cluster["found_rt"], errors="coerce")
        area = pd.to_numeric(cluster["area"], errors="coerce")
        matched_peak_id = pd.to_numeric(cluster["matched_peak_id"], errors="coerce")
        score += 5.0 * float(found_rt.notna().sum())
        score -= 8.0 * float(area.isna().sum())

        peak_ids = matched_peak_id.dropna().astype(int)
        score -= 10.0 * float(peak_ids.duplicated().sum())

        rt_values = found_rt.dropna().to_numpy(dtype=float)
        if len(rt_values) >= 2:
            score -= 1000.0 * float(np.sum(np.maximum(0.0, 0.004 - np.diff(rt_values))))
    return float(score)


def _annotate_processing_result(result: dict, baseline_mode: str) -> dict:
    annotated = dict(result)
    annotated["baseline_mode"] = baseline_mode
    annotated["cluster_quality_score"] = _compute_cluster_quality_score(annotated["matched_targets_df"])
    annotated["confidence"] = build_confidence_assessment(
        matched_targets_df=annotated["matched_targets_df"],
        peaks_df=annotated["peaks_df"],
        omega=annotated["omega"],
        baseline_mode=annotated["baseline_mode"],
        cluster_quality_score=annotated["cluster_quality_score"],
    )
    return annotated


def _process_chromatogram_from_baseline(processed_df: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    processed_df, best_window = add_savgol_and_derivatives_to_dataframe(
        processed_df,
        polyorder=SAVGOL_POLYORDER,
        candidate_windows=SAVGOL_CANDIDATE_WINDOWS,
    )
    raw_peaks_df = detect_peaks_from_derivatives(
        processed_df,
        best_window=best_window,
        height_sigma=PEAK_DETECTION_HEIGHT_SIGMA,
        prominence_sigma=PEAK_DETECTION_PROMINENCE_SIGMA,
    )
    peaks_df = raw_peaks_df
    matched_targets_df, rt_shift = dynamic_match_targets_to_peaks(reference_targets, peaks_df)
    peaks_df, matched_targets_df = refine_c18_c20_cluster_matches(processed_df, peaks_df, matched_targets_df)
    matched_targets_df = refine_overlapped_c22_cluster_areas(processed_df, peaks_df, matched_targets_df)
    matched_targets_df = refine_cluster_areas_by_local_valleys(processed_df, matched_targets_df)
    matched_targets_df = recover_missing_c22_components_with_fit(processed_df, peaks_df, matched_targets_df)
    matched_targets_df = recover_underintegrated_c20_components_with_fit(processed_df, peaks_df, matched_targets_df)
    matched_targets_df = recover_overlapped_c18_components_with_fit(processed_df, peaks_df, matched_targets_df)
    matched_targets_df = tighten_overwide_c22_cluster_tails(processed_df, matched_targets_df)
    matched_targets_df = refine_overwide_c22_cluster_with_pvfit(processed_df, peaks_df, matched_targets_df)
    matched_targets_df = refine_small_peak_integrations(processed_df, matched_targets_df)
    omega = compute_omega_metrics(matched_targets_df)
    return {
        "processed_df": processed_df,
        "best_window": best_window,
        "peaks_df": peaks_df,
        "matched_targets_df": matched_targets_df,
        "rt_shift": rt_shift,
        "omega": omega,
        "omega_report": omega["omega3_trio"],
    }


def estimate_rt_shift(expected_rts: np.ndarray, observed_rts: np.ndarray, max_shift: float = 0.15, step: float = 0.001, sigma: float = 0.025) -> float:
    expected_rts = np.asarray(expected_rts, dtype=float)
    observed_rts = np.asarray(observed_rts, dtype=float)
    if expected_rts.size == 0 or observed_rts.size == 0:
        return 0.0
    shifts = np.arange(-max_shift, max_shift + step, step)
    best_shift, best_score = 0.0, -np.inf
    for shift in shifts:
        score = 0.0
        for rt in expected_rts + shift:
            d = np.min(np.abs(observed_rts - rt))
            score += np.exp(-0.5 * (d / sigma) ** 2)
        if score > best_score:
            best_score, best_shift = score, float(shift)
    return best_shift


def _clear_match(out: pd.DataFrame, row_idx: int):
    out.at[row_idx, "found_rt"] = np.nan
    out.at[row_idx, "area"] = np.nan
    out.at[row_idx, "percent_area"] = np.nan
    out.at[row_idx, "matched_peak_id"] = np.nan
    out.at[row_idx, "match_score"] = np.nan
    out.at[row_idx, "status"] = "not_found"


def _assign_peak(out: pd.DataFrame, row_idx: int, peak_row: pd.Series, status: str, match_score: float):
    out.at[row_idx, "found_rt"] = float(peak_row["apex_x"])
    out.at[row_idx, "area"] = float(peak_row["area"])
    out.at[row_idx, "percent_area"] = float(peak_row["percent_area"])
    out.at[row_idx, "matched_peak_id"] = int(peak_row["peak_id"])
    out.at[row_idx, "match_score"] = float(match_score)
    out.at[row_idx, "status"] = status
    out.at[row_idx, "integration_start_x"] = float(peak_row["start_x"])
    out.at[row_idx, "integration_end_x"] = float(peak_row["end_x"])


def _get_order_bounds(out: pd.DataFrame, row_idx: int):
    lower = -np.inf
    upper = np.inf

    for j in range(row_idx - 1, -1, -1):
        if pd.notna(out.at[j, "found_rt"]):
            lower = float(out.at[j, "found_rt"])
            break

    for j in range(row_idx + 1, len(out)):
        if pd.notna(out.at[j, "found_rt"]):
            upper = float(out.at[j, "found_rt"])
            break

    return lower, upper


def _select_reliable_rt_peak(candidates: pd.DataFrame) -> int:
    nearest_idx = int(candidates["distance"].idxmin())
    dominant_idx = int(candidates["area"].idxmax())
    nearest_area = float(candidates.loc[nearest_idx, "area"])
    dominant_area = float(candidates.loc[dominant_idx, "area"])
    dominant_distance = float(candidates.loc[dominant_idx, "distance"])
    dominant_is_clear = dominant_area > max(
        nearest_area * RELIABLE_RT_DOMINANT_AREA_MULTIPLIER,
        nearest_area + RELIABLE_RT_DOMINANT_AREA_MIN_DELTA,
    )
    if dominant_is_clear and dominant_distance <= RELIABLE_RT_DOMINANT_DISTANCE_MAX:
        return dominant_idx
    return nearest_idx


def _select_best_peak_position(
    candidate_positions: np.ndarray,
    distances: np.ndarray,
    prominences: np.ndarray,
    areas: np.ndarray,
) -> int:
    return int(min(
        candidate_positions.tolist(),
        key=lambda pos: (float(distances[pos]), -float(prominences[pos]), -float(areas[pos]), int(pos)),
    ))


def _select_reliable_rt_peak_position(
    candidate_positions: np.ndarray,
    distances: np.ndarray,
    areas: np.ndarray,
) -> int:
    nearest_pos = int(candidate_positions[np.argmin(distances[candidate_positions])])
    dominant_pos = int(candidate_positions[np.argmax(areas[candidate_positions])])
    nearest_area = float(areas[nearest_pos])
    dominant_area = float(areas[dominant_pos])
    dominant_distance = float(distances[dominant_pos])
    dominant_is_clear = dominant_area > max(
        nearest_area * RELIABLE_RT_DOMINANT_AREA_MULTIPLIER,
        nearest_area + RELIABLE_RT_DOMINANT_AREA_MIN_DELTA,
    )
    if dominant_is_clear and dominant_distance <= RELIABLE_RT_DOMINANT_DISTANCE_MAX:
        return dominant_pos
    return nearest_pos


def _apply_target_cluster_override(
    matched_targets_df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    cluster_codes,
    target_apexes,
    max_distance: float,
    status: str,
    rt_shift: float,
) -> pd.DataFrame:
    out = matched_targets_df.copy()
    if out.empty or peaks_df is None or peaks_df.empty:
        return out

    peaks_apex = peaks_df["apex_x"].to_numpy(dtype=float)
    peaks_prominence = peaks_df["prominence"].to_numpy(dtype=float)
    peaks_area = peaks_df["area"].to_numpy(dtype=float)
    peak_ids = peaks_df["peak_id"].to_numpy(dtype=int)
    used_mask = np.zeros(len(peaks_df), dtype=bool)
    chosen = {}
    for code, target_apex in zip(cluster_codes, target_apexes):
        adjusted_target = float(target_apex + rt_shift)
        distances = np.abs(peaks_apex - adjusted_target)
        candidate_positions = np.flatnonzero((~used_mask) & (distances <= max_distance))
        if candidate_positions.size == 0:
            continue
        best_pos = _select_best_peak_position(candidate_positions, distances, peaks_prominence, peaks_area)
        peak_row = peaks_df.iloc[best_pos]
        chosen[code] = peak_row
        used_mask[best_pos] = True

    if not chosen:
        return out

    selected_peak_ids = {int(peak["peak_id"]) for peak in chosen.values()}
    conflict_mask = out["matched_peak_id"].isin(selected_peak_ids) & ~out["code"].isin(cluster_codes)
    for row_idx in out.index[conflict_mask]:
        _clear_match(out, int(row_idx))

    for row_idx, row in out.iterrows():
        peak_row = chosen.get(row["code"])
        if peak_row is None:
            continue
        adjusted_target = float(target_apexes[list(cluster_codes).index(row["code"])] + rt_shift)
        _assign_peak(out, row_idx, peak_row, status, abs(float(peak_row["apex_x"]) - adjusted_target))

    return out


def apply_c22_cluster_override(
    matched_targets_df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    rt_shift: float = 0.0,
) -> pd.DataFrame:
    out = matched_targets_df.copy()
    if out.empty or peaks_df is None or peaks_df.empty:
        return out

    peaks_apex = peaks_df["apex_x"].to_numpy(dtype=float)
    peaks_prominence = peaks_df["prominence"].to_numpy(dtype=float)
    peaks_area = peaks_df["area"].to_numpy(dtype=float)
    cluster_codes = ["C22:6", "C22:5", "C22:4"]
    target_apexes = [9.252, 9.285, 9.316]
    max_distances = [0.020, 0.020, 0.020]

    chosen = {}
    used_mask = np.zeros(len(peaks_df), dtype=bool)
    for code, target_apex, max_distance in zip(cluster_codes, target_apexes, max_distances):
        adjusted_target = float(target_apex + rt_shift)
        distances = np.abs(peaks_apex - adjusted_target)
        candidate_positions = np.flatnonzero((~used_mask) & (distances <= max_distance))
        if candidate_positions.size == 0:
            continue
        best_pos = _select_best_peak_position(candidate_positions, distances, peaks_prominence, peaks_area)
        peak_row = peaks_df.iloc[best_pos]
        chosen[code] = peak_row
        used_mask[best_pos] = True

    if len(chosen) < 2:
        return out

    for code in cluster_codes:
        row_idx = out.index[out["code"] == code]
        if len(row_idx):
            _clear_match(out, int(row_idx[0]))

    selected_peak_ids = {int(peak["peak_id"]) for peak in chosen.values()}
    conflict_mask = out["matched_peak_id"].isin(selected_peak_ids) & ~out["code"].isin(cluster_codes)
    for row_idx in out.index[conflict_mask]:
        _clear_match(out, int(row_idx))

    for code, target_apex in zip(cluster_codes, target_apexes):
        peak_row = chosen.get(code)
        if peak_row is None:
            continue
        row_idx = out.index[out["code"] == code][0]
        adjusted_target = float(target_apex + rt_shift)
        _assign_peak(out, row_idx, peak_row, "matched_c22_rule", abs(float(peak_row["apex_x"]) - adjusted_target))

    return out


def dynamic_match_targets_to_peaks(targets_df: pd.DataFrame, peaks_df: pd.DataFrame):
    targets = targets_df.sort_values("order_index").reset_index(drop=True)
    if peaks_df is None or peaks_df.empty:
        out = targets.copy()
        for col in ["corrected_target_rt", "found_rt", "area", "percent_area", "matched_peak_id", "match_score"]:
            out[col] = np.nan
        out["status"] = "not_found"
        return out, 0.0

    peaks = peaks_df.sort_values("apex_x").reset_index(drop=True)
    peak_apex = peaks["apex_x"].to_numpy(dtype=float)
    peak_area = peaks["area"].to_numpy(dtype=float)
    peak_prominence = peaks["prominence"].to_numpy(dtype=float)
    refs = targets[targets["rt_reliable"] & targets["expected_rt"].notna()]
    shift = estimate_rt_shift(refs["expected_rt"].to_numpy(dtype=float), peak_apex) if not refs.empty else 0.0
    targets["corrected_target_rt"] = targets["expected_rt"] + shift

    out = targets.copy()
    out["found_rt"] = np.nan
    out["area"] = np.nan
    out["percent_area"] = np.nan
    out["matched_peak_id"] = np.nan
    out["match_score"] = np.nan
    out["integration_start_x"] = np.nan
    out["integration_end_x"] = np.nan
    out["status"] = "not_found"

    used_mask = np.zeros(len(peaks), dtype=bool)
    reliable_rows = out.index[out["rt_reliable"] & out["corrected_target_rt"].notna()].tolist()
    for i in reliable_rows:
        expected = float(out.at[i, "corrected_target_rt"])
        distances = np.abs(peak_apex - expected)
        candidate_positions = np.flatnonzero(~used_mask)
        if candidate_positions.size == 0:
            continue
        local_positions = candidate_positions[distances[candidate_positions] <= RELIABLE_RT_WINDOW]
        if local_positions.size == 0:
            best_pos = int(candidate_positions[np.argmin(distances[candidate_positions])])
            best_distance = float(distances[best_pos])
        else:
            best_pos = _select_reliable_rt_peak_position(local_positions, distances, peak_area)
            best_distance = float(distances[best_pos])
        if best_distance > 0.3:
            continue
        used_mask[best_pos] = True
        _assign_peak(out, i, peaks.iloc[best_pos], "matched_rt", best_distance)

    soft_rows = out.index[(~out["rt_reliable"]) & out["corrected_target_rt"].notna()].tolist()
    for i in soft_rows:
        lower, upper = _get_order_bounds(out, i)
        expected = float(out.at[i, "corrected_target_rt"])
        distances = np.abs(peak_apex - expected)
        candidate_positions = np.flatnonzero((~used_mask) & (peak_apex > lower) & (peak_apex < upper))
        if candidate_positions.size == 0:
            continue
        best_pos = int(min(candidate_positions.tolist(), key=lambda pos: (float(distances[pos]), int(pos))))
        best_distance = float(distances[best_pos])
        if best_distance > 0.2:
            continue
        used_mask[best_pos] = True
        _assign_peak(out, i, peaks.iloc[best_pos], "matched_soft_rt", best_distance)

    order_rows = out.index[out["corrected_target_rt"].isna()].tolist()
    for i in order_rows:
        lower, upper = _get_order_bounds(out, i)
        candidate_positions = np.flatnonzero((~used_mask) & (peak_apex > lower) & (peak_apex < upper))
        if candidate_positions.size == 0:
            continue
        best_pos = int(candidate_positions[0])
        used_mask[best_pos] = True
        _assign_peak(out, i, peaks.iloc[best_pos], "matched_order", 0.0)

    out = _apply_target_cluster_override(
        out,
        peaks_df,
        cluster_codes=["C18:1N9C", "C18:3N3", "C18:0"],
        target_apexes=[7.623, 7.650, 7.750],
        max_distance=0.025,
        status="matched_c18_rule",
        rt_shift=shift,
    )
    out = _apply_target_cluster_override(
        out,
        peaks_df,
        cluster_codes=["C20:4N6", "C20:5", "C20:3N8"],
        target_apexes=[8.381, 8.410, 8.467],
        max_distance=0.025,
        status="matched_c20_rule",
        rt_shift=shift,
    )
    out = apply_c22_cluster_override(out, peaks_df, rt_shift=shift)
    return out, shift


def _gaussian_component(x: np.ndarray, amplitude: float, center: float, sigma: float) -> np.ndarray:
    sigma = max(float(sigma), 1e-6)
    return float(amplitude) * np.exp(-0.5 * ((x - float(center)) / sigma) ** 2)


def _split_pseudo_voigt_unit_area(
    x: np.ndarray,
    center: float,
    fwhm_left: float,
    fwhm_right: float,
    eta: float,
) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    center = float(center)
    fwhm_left = max(float(fwhm_left), 1e-6)
    fwhm_right = max(float(fwhm_right), 1e-6)
    eta = float(np.clip(eta, 0.0, 1.0))
    profile = np.zeros_like(x, dtype=float)

    left_mask = x <= center
    right_mask = ~left_mask
    left_scale = 2.0 * fwhm_left / (fwhm_left + fwhm_right)
    right_scale = 2.0 * fwhm_right / (fwhm_left + fwhm_right)
    if np.any(left_mask):
        profile[left_mask] = left_scale * _pseudo_voigt_unit_area(x[left_mask], center=center, fwhm=fwhm_left, eta=eta)
    if np.any(right_mask):
        profile[right_mask] = right_scale * _pseudo_voigt_unit_area(x[right_mask], center=center, fwhm=fwhm_right, eta=eta)
    return profile


def _build_cluster_fit_weights(
    x: np.ndarray,
    y: np.ndarray,
    centers: np.ndarray,
    noise: float,
) -> np.ndarray:
    positive_y = np.clip(np.asarray(y, dtype=float), 0.0, None)
    noise = max(float(noise), 1.0)
    weights = 1.0 / np.maximum(np.sqrt(positive_y + noise * noise), noise)
    low_level = max(float(np.quantile(positive_y, 0.50)), noise)
    foot_fraction = 1.0 - np.clip(positive_y / max(low_level * 2.0, noise * 2.0), 0.0, 1.0)
    weights *= 1.0 + SPLIT_PSEUDOVOIGT_FOOT_WEIGHT_BOOST * foot_fraction

    if len(centers) > 1:
        for midpoint in 0.5 * (centers[:-1] + centers[1:]):
            weights *= 1.0 + SPLIT_PSEUDOVOIGT_VALLEY_WEIGHT_BOOST * np.exp(-0.5 * ((x - float(midpoint)) / 0.010) ** 2)
    weights *= 1.0 + SPLIT_PSEUDOVOIGT_EDGE_WEIGHT_BOOST * np.exp(-0.5 * ((x - float(centers[0])) / 0.018) ** 2)
    weights *= 1.0 + SPLIT_PSEUDOVOIGT_EDGE_WEIGHT_BOOST * np.exp(-0.5 * ((x - float(centers[-1])) / 0.018) ** 2)
    return weights


def _snap_boundary_to_change_point(
    x: np.ndarray,
    boundary_metric: np.ndarray,
    d2y: np.ndarray,
    target_x: float,
    left_limit_idx: int,
    right_limit_idx: int,
    search_half_window: float = SPLIT_PSEUDOVOIGT_BOUNDARY_SNAP_WINDOW,
) -> int:
    left_limit_idx = int(max(0, left_limit_idx))
    right_limit_idx = int(min(len(x) - 1, right_limit_idx))
    if right_limit_idx <= left_limit_idx:
        return left_limit_idx

    search_left = max(left_limit_idx, int(np.searchsorted(x, float(target_x - search_half_window), side="left")))
    search_right = min(right_limit_idx, int(np.searchsorted(x, float(target_x + search_half_window), side="right") - 1))
    if search_right <= search_left:
        return _find_preferred_minimum_index(
            boundary_metric,
            left_limit_idx,
            right_limit_idx,
            target_idx=float(np.argmin(np.abs(x[left_limit_idx:right_limit_idx + 1] - float(target_x))) + left_limit_idx),
        )

    metric_scale = max(float(np.nanmax(boundary_metric[search_left:search_right + 1])), 1.0)
    candidate_indices = []
    for idx in range(max(search_left + 1, 1), min(search_right, len(x) - 2)):
        local_min = boundary_metric[idx - 1] >= boundary_metric[idx] <= boundary_metric[idx + 1]
        d2_cross = d2y[idx] == 0.0 or np.sign(d2y[idx]) != np.sign(d2y[idx - 1])
        if local_min or d2_cross:
            candidate_indices.append((idx, local_min, d2_cross))
    if not candidate_indices:
        candidate_indices = [
            (idx, False, False)
            for idx in range(search_left, search_right + 1)
        ]

    span = max(float(search_half_window), 1e-6)

    def score(item):
        idx, local_min, d2_cross = item
        distance = abs(float(x[idx]) - float(target_x)) / span
        metric_term = float(boundary_metric[idx]) / metric_scale
        feature_penalty = 0.0 if (local_min and d2_cross) else (0.04 if (local_min or d2_cross) else 0.10)
        return metric_term + 0.22 * distance + feature_penalty

    best_idx = int(min(candidate_indices, key=score)[0])
    return int(min(max(best_idx, left_limit_idx), right_limit_idx))


def _derive_split_pseudovoigt_boundaries(
    df: pd.DataFrame,
    fitted_components,
    window_left: float,
    window_right: float,
):
    if not fitted_components:
        return []

    x_col = _get_x_column_name(df)
    sub = df[(df[x_col] >= window_left) & (df[x_col] <= window_right)].copy()
    if len(sub) < 6:
        return []

    x = sub[x_col].to_numpy(dtype=float)
    y_corrected_raw = sub["y_corrected"].to_numpy(dtype=float)
    y_smooth = sub["y_smooth"].to_numpy(dtype=float)
    d2y = sub["d2y"].to_numpy(dtype=float) if "d2y" in sub.columns else np.zeros(len(sub), dtype=float)
    baseline = _estimate_local_linear_baseline(x, y_corrected_raw, 0, len(x) - 1)
    corrected_local = np.clip(y_corrected_raw - baseline, 0.0, None)
    smooth_local = np.clip(y_smooth - baseline, 0.0, None)
    boundary_metric = 0.70 * smooth_local + 0.30 * corrected_local
    cluster_noise = max(_robust_sigma(y_corrected_raw), 1.0)

    ordered = sorted(fitted_components, key=lambda item: float(item["center"]))
    total_model = np.zeros_like(x, dtype=float)
    for component in ordered:
        total_model += float(component["area"]) * _split_pseudo_voigt_unit_area(
            x,
            center=component["center"],
            fwhm_left=component["fwhm_left"],
            fwhm_right=component["fwhm_right"],
            eta=component["eta"],
        )
    support_signal = np.maximum(boundary_metric, total_model)

    left_component = ordered[0]
    right_component = ordered[-1]
    left_seed_x = max(float(window_left), float(left_component["center"] - SPLIT_PSEUDOVOIGT_OUTER_SUPPORT_WIDTH_FACTOR * left_component["fwhm_left"]))
    right_seed_x = min(float(window_right), float(right_component["center"] + SPLIT_PSEUDOVOIGT_OUTER_SUPPORT_WIDTH_FACTOR * right_component["fwhm_right"]))
    left_seed_idx = int(np.argmin(np.abs(x - left_seed_x)))
    right_seed_idx = int(np.argmin(np.abs(x - right_seed_x)))

    left_peak_height = float(np.max(total_model[: max(int(np.argmin(np.abs(x - left_component["center"]))) + 1, 1)]))
    right_peak_height = float(np.max(total_model[min(int(np.argmin(np.abs(x - right_component["center"]))), len(x) - 1) :]))
    left_threshold = max(cluster_noise * PEAK_SUPPORT_THRESHOLD_SIGMA, left_peak_height * PEAK_SUPPORT_THRESHOLD_FRACTION)
    right_threshold = max(cluster_noise * PEAK_SUPPORT_THRESHOLD_SIGMA, right_peak_height * PEAK_SUPPORT_THRESHOLD_FRACTION)

    left_boundary = _extend_boundary_to_support(
        signal=support_signal,
        start_idx=left_seed_idx,
        limit_idx=0,
        direction=-1,
        threshold=float(left_threshold),
        consecutive_points=PEAK_SUPPORT_CONSECUTIVE_POINTS,
    )
    right_boundary = _extend_boundary_to_support(
        signal=support_signal,
        start_idx=right_seed_idx,
        limit_idx=len(x) - 1,
        direction=1,
        threshold=float(right_threshold),
        consecutive_points=PEAK_SUPPORT_CONSECUTIVE_POINTS,
    )

    boundaries = [int(left_boundary)]
    for left_component, right_component in zip(ordered[:-1], ordered[1:]):
        dense_x = np.linspace(float(left_component["center"]), float(right_component["center"]), 360)
        left_curve = float(left_component["area"]) * _split_pseudo_voigt_unit_area(
            dense_x,
            center=left_component["center"],
            fwhm_left=left_component["fwhm_left"],
            fwhm_right=left_component["fwhm_right"],
            eta=left_component["eta"],
        )
        right_curve = float(right_component["area"]) * _split_pseudo_voigt_unit_area(
            dense_x,
            center=right_component["center"],
            fwhm_left=right_component["fwhm_left"],
            fwhm_right=right_component["fwhm_right"],
            eta=right_component["eta"],
        )
        diff = left_curve - right_curve
        crossing_x = None
        for idx in range(1, len(diff)):
            if diff[idx - 1] == 0.0:
                crossing_x = float(dense_x[idx - 1])
                break
            if diff[idx] == 0.0 or np.sign(diff[idx]) != np.sign(diff[idx - 1]):
                crossing_x = float(dense_x[idx])
                break
        if crossing_x is None:
            crossing_x = 0.5 * (float(left_component["center"]) + float(right_component["center"]))

        left_limit_idx = max(boundaries[-1] + 1, int(np.argmin(np.abs(x - float(left_component["center"])))))
        right_limit_idx = max(left_limit_idx + 1, int(np.argmin(np.abs(x - float(right_component["center"])))))
        split_idx = _snap_boundary_to_change_point(
            x=x,
            boundary_metric=boundary_metric,
            d2y=d2y,
            target_x=float(crossing_x),
            left_limit_idx=left_limit_idx,
            right_limit_idx=right_limit_idx,
        )
        if split_idx <= boundaries[-1]:
            split_idx = max(boundaries[-1] + 1, int(round(0.5 * (left_limit_idx + right_limit_idx))))
        boundaries.append(int(split_idx))
    boundaries.append(int(right_boundary))

    if any(boundaries[i] >= boundaries[i + 1] for i in range(len(boundaries) - 1)):
        return []

    resolved = []
    for idx, component in enumerate(ordered):
        start_idx = int(boundaries[idx])
        end_idx = int(boundaries[idx + 1])
        if end_idx <= start_idx:
            return []
        resolved.append((float(x[start_idx]), float(x[end_idx])))
    return resolved


def _derive_fit_component_boundaries(
    fitted_components,
    window_left: float,
    window_right: float,
):
    if not fitted_components:
        return []

    components = sorted(fitted_components, key=lambda item: float(item["center"]))
    boundaries = [float(window_left)]
    for i in range(len(components) - 1):
        left = components[i]
        right = components[i + 1]
        dense_x = np.linspace(float(left["center"]), float(right["center"]), 240)
        left_curve = _gaussian_component(dense_x, left["amplitude"], left["center"], left["sigma"])
        right_curve = _gaussian_component(dense_x, right["amplitude"], right["center"], right["sigma"])
        diff = left_curve - right_curve
        sign = np.sign(diff)
        crossing_idx = None
        for j in range(1, len(sign)):
            if sign[j - 1] == 0:
                crossing_idx = j - 1
                break
            if sign[j] == 0 or sign[j] != sign[j - 1]:
                crossing_idx = j
                break
        if crossing_idx is None:
            boundary = 0.5 * (float(left["center"]) + float(right["center"]))
        else:
            boundary = float(dense_x[crossing_idx])
        boundaries.append(boundary)
    boundaries.append(float(window_right))

    resolved = []
    for i in range(len(components)):
        start_x = float(boundaries[i])
        end_x = float(boundaries[i + 1])
        if end_x <= start_x:
            center_left = float(components[i]["center"])
            center_right = float(components[i + 1]["center"]) if i + 1 < len(components) else float(window_right)
            end_x = max(start_x + 1e-4, 0.5 * (center_left + center_right))
        resolved.append((start_x, end_x))
    return resolved


def _fit_cluster_components_split_pseudovoigt(
    df: pd.DataFrame,
    initial_centers,
    window_left: float,
    window_right: float,
    center_tolerances,
    initial_areas=None,
    sigma_bounds=(0.003, 0.025),
):
    if not ENABLE_SPLIT_PSEUDOVOIGT_CLUSTER_FIT:
        return None, {}

    x_col = _get_x_column_name(df)
    sub = df[(df[x_col] >= window_left) & (df[x_col] <= window_right)].copy()
    if len(sub) < 20:
        return None, {}

    x = sub[x_col].to_numpy(dtype=float)
    y = np.clip(sub["y_corrected"].to_numpy(dtype=float), 0.0, None)
    if not np.any(y > 0):
        return None, {}

    initial_centers = np.asarray(initial_centers, dtype=float)
    center_tolerances = np.asarray(center_tolerances, dtype=float)
    if initial_areas is None:
        initial_areas = [np.nan] * len(initial_centers)

    spacing = np.diff(initial_centers)
    base_sigma = 0.0075 if len(spacing) == 0 else float(np.clip(np.min(spacing) * 0.28, 0.0045, 0.010))
    base_fwhm = float(np.clip(2.0 * base_sigma, 0.0075, 0.024))
    fwhm_min = max(2.0 * float(sigma_bounds[0]), base_fwhm * 0.65)
    fwhm_max = min(2.0 * float(sigma_bounds[1]), max(base_fwhm * 1.55, fwhm_min + 1e-4))
    noise = max(_robust_sigma(sub["y_corrected"].to_numpy(dtype=float)), 1.0)
    baseline_init = _estimate_local_linear_baseline(x, y, 0, len(x) - 1)
    if baseline_init.size == 0:
        return None, {}
    mid = float(np.mean(x))
    baseline_intercept = float(np.median(baseline_init))
    baseline_slope = 0.0 if len(x) < 2 else float((baseline_init[-1] - baseline_init[0]) / max(x[-1] - x[0], 1e-9))
    weights = _build_cluster_fit_weights(x, y, initial_centers, noise=noise)
    min_dx = float(np.median(np.diff(x))) if len(x) > 1 else 0.001

    amplitude_init = []
    for center, prior_area in zip(initial_centers, initial_areas):
        nearest_idx = int(np.argmin(np.abs(x - center)))
        peak_height = max(float(y[nearest_idx] - baseline_init[nearest_idx]), max(float(np.max(y)) * 0.01, 1.0))
        rough_area = peak_height * base_fwhm
        amplitude_guess = float(prior_area) if np.isfinite(prior_area) and float(prior_area) > 0 else rough_area
        amplitude_init.append(max(amplitude_guess, rough_area, 1.0))

    n = len(initial_centers)
    p0 = np.array(
        amplitude_init
        + list(initial_centers)
        + [base_fwhm] * n
        + [base_fwhm] * n
        + [0.35]
        + [baseline_intercept, baseline_slope],
        dtype=float,
    )
    lower_bounds = np.array(
        [0.0] * n
        + list(initial_centers - center_tolerances)
        + [fwhm_min] * n
        + [fwhm_min] * n
        + [0.0]
        + [0.0, -max(float(np.max(y)) * 4.0, 200.0)],
        dtype=float,
    )
    upper_bounds = np.array(
        [max(float(np.max(y)) * 60.0, max(amplitude_init) * 3.5, 10.0)] * n
        + list(initial_centers + center_tolerances)
        + [fwhm_max] * n
        + [fwhm_max] * n
        + [1.0]
        + [max(float(np.max(y)) * 0.20, 80.0), max(float(np.max(y)) * 4.0, 200.0)],
        dtype=float,
    )

    def unpack(params):
        amplitudes = params[:n]
        centers = params[n : 2 * n]
        fwhm_left = params[2 * n : 3 * n]
        fwhm_right = params[3 * n : 4 * n]
        eta = float(params[4 * n])
        baseline_0 = float(params[-2])
        baseline_1 = float(params[-1])
        return amplitudes, centers, fwhm_left, fwhm_right, eta, baseline_0, baseline_1

    def residuals(params):
        amplitudes, centers, fwhm_left, fwhm_right, eta, baseline_0, baseline_1 = unpack(params)
        if np.any(np.diff(centers) <= min_dx * 0.5):
            return np.full(n + len(x), 1e6, dtype=float)

        baseline = baseline_0 + baseline_1 * (x - mid)
        prediction = baseline.copy()
        penalties = np.zeros(n, dtype=float)
        asymmetry_cap = max(abs(math.log(SPLIT_PSEUDOVOIGT_ASYMMETRY_MIN)), abs(math.log(SPLIT_PSEUDOVOIGT_ASYMMETRY_MAX)))
        for amplitude, center, left_width, right_width in zip(amplitudes, centers, fwhm_left, fwhm_right):
            prediction += float(amplitude) * _split_pseudo_voigt_unit_area(
                x,
                center=float(center),
                fwhm_left=float(left_width),
                fwhm_right=float(right_width),
                eta=float(eta),
            )
        for idx, (left_width, right_width) in enumerate(zip(fwhm_left, fwhm_right)):
            asymmetry = left_width / max(right_width, 1e-9)
            asymmetry_log = abs(math.log(max(asymmetry, 1e-9)))
            if asymmetry_log > asymmetry_cap:
                penalties[idx] = (asymmetry_log - asymmetry_cap) * 80.0
        residual = (prediction - y) * weights
        return np.concatenate([residual, penalties])

    result = least_squares(
        residuals,
        p0,
        bounds=(lower_bounds, upper_bounds),
        max_nfev=18000,
    )
    if not result.success:
        return None, {}

    amplitudes, centers, fwhm_left, fwhm_right, eta, baseline_0, baseline_1 = unpack(result.x)
    baseline = baseline_0 + baseline_1 * (x - mid)
    prediction = baseline.copy()
    fitted_components = []
    for amplitude, center, left_width, right_width in zip(amplitudes, centers, fwhm_left, fwhm_right):
        component_curve = float(amplitude) * _split_pseudo_voigt_unit_area(
            x,
            center=float(center),
            fwhm_left=float(left_width),
            fwhm_right=float(right_width),
            eta=float(eta),
        )
        prediction += component_curve
        fitted_components.append({
            "center": float(center),
            "area": float(np.trapezoid(component_curve, x)),
            "amplitude": float(amplitude),
            "fwhm_left": float(left_width),
            "fwhm_right": float(right_width),
            "fwhm": float(0.5 * (left_width + right_width)),
            "eta": float(eta),
        })

    residual = y - prediction
    ss_res = float(np.sum(residual * residual))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    if not np.isfinite(r2) or r2 < SPLIT_PSEUDOVOIGT_MIN_R2:
        return None, {"r2": r2}

    return fitted_components, {
        "r2": float(r2),
        "baseline_intercept": float(baseline_0),
        "baseline_slope": float(baseline_1),
        "split_pv": True,
    }


def _fit_cluster_components_gaussian(
    df: pd.DataFrame,
    initial_centers,
    window_left: float,
    window_right: float,
    center_tolerances,
    sigma_bounds=(0.003, 0.025),
):
    x_col = _get_x_column_name(df)
    sub = df[(df[x_col] >= window_left) & (df[x_col] <= window_right)].copy()
    if len(sub) < 20:
        return None

    x = sub[x_col].to_numpy(dtype=float)
    y = np.clip(sub["y_smooth"].to_numpy(dtype=float), 0.0, None)
    if not np.any(y > 0):
        return None

    initial_centers = np.asarray(initial_centers, dtype=float)
    center_tolerances = np.asarray(center_tolerances, dtype=float)
    floor = float(np.quantile(y, 0.05))
    mid = float(np.mean(x))
    min_dx = float(np.median(np.diff(x))) if len(x) > 1 else 0.001
    weights = 1.0 / np.sqrt(np.maximum(y, np.quantile(y, 0.35)) + 1.0)

    amplitude_init = []
    for center in initial_centers:
        idx = int(np.argmin(np.abs(x - center)))
        amplitude_init.append(max(float(y[idx] - floor), max(y) * 0.01, 1.0))

    spacing = np.diff(initial_centers)
    base_sigma = 0.0075 if len(spacing) == 0 else float(np.clip(np.min(spacing) * 0.32, 0.005, 0.012))
    sigma_init = [base_sigma] * len(initial_centers)

    p0 = np.array(amplitude_init + list(initial_centers) + sigma_init + [floor, 0.0], dtype=float)
    lower_bounds = np.array(
        [0.0] * len(initial_centers)
        + list(initial_centers - center_tolerances)
        + [sigma_bounds[0]] * len(initial_centers)
        + [0.0, -max(y) * 20.0],
        dtype=float,
    )
    upper_bounds = np.array(
        [max(y) * 20.0] * len(initial_centers)
        + list(initial_centers + center_tolerances)
        + [sigma_bounds[1]] * len(initial_centers)
        + [max(y) * 2.0, max(y) * 20.0],
        dtype=float,
    )

    def residuals(params):
        n = len(initial_centers)
        amplitudes = params[:n]
        centers = params[n : 2 * n]
        sigmas = params[2 * n : 3 * n]
        baseline_0 = params[-2]
        baseline_1 = params[-1]

        if np.any(np.diff(centers) <= min_dx * 0.5):
            return np.full_like(x, 1e6)

        baseline = baseline_0 + baseline_1 * (x - mid)
        prediction = baseline.copy()
        for amplitude, center, sigma in zip(amplitudes, centers, sigmas):
            prediction += _gaussian_component(x, amplitude, center, sigma)
        return (prediction - y) * weights

    result = least_squares(residuals, p0, bounds=(lower_bounds, upper_bounds), max_nfev=12000)
    if not result.success:
        return None

    params = result.x
    n = len(initial_centers)
    amplitudes = params[:n]
    centers = params[n : 2 * n]
    sigmas = params[2 * n : 3 * n]
    fitted_components = []
    for amplitude, center, sigma in zip(amplitudes, centers, sigmas):
        area = float(amplitude * sigma * math.sqrt(2.0 * math.pi))
        fitted_components.append({
            "center": float(center),
            "area": area,
            "amplitude": float(amplitude),
            "sigma": float(sigma),
        })

    return fitted_components


def _fit_cluster_components_lmfit_pseudovoigt(
    df: pd.DataFrame,
    initial_centers,
    window_left: float,
    window_right: float,
    center_tolerances,
    initial_areas=None,
    sigma_bounds=(0.003, 0.025),
):
    if not ENABLE_LMFIT_LOCAL_PSEUDOVOIGT or PseudoVoigtModel is None or LinearModel is None:
        return None, {}

    x_col = _get_x_column_name(df)
    sub = df[(df[x_col] >= window_left) & (df[x_col] <= window_right)].copy()
    if len(sub) < 20:
        return None, {}

    x = sub[x_col].to_numpy(dtype=float)
    y_raw = np.clip(sub["y_corrected"].to_numpy(dtype=float), 0.0, None)
    if not np.any(y_raw > 0):
        return None, {}

    initial_centers = np.asarray(initial_centers, dtype=float)
    center_tolerances = np.asarray(center_tolerances, dtype=float)
    if initial_areas is None:
        initial_areas = [np.nan] * len(initial_centers)

    baseline = _estimate_local_linear_baseline(x, y_raw, 0, len(x) - 1)
    y = np.clip(y_raw - baseline, 0.0, None)
    if not np.any(y > 0):
        return None, {}

    floor = 0.0
    spacing = np.diff(initial_centers)
    base_sigma = 0.0075 if len(spacing) == 0 else float(np.clip(np.min(spacing) * 0.28, 0.0045, 0.010))
    sigma_min = max(float(sigma_bounds[0]), base_sigma * 0.78)
    sigma_max = min(float(sigma_bounds[1]), max(base_sigma * 1.30, sigma_min + 1e-4))
    noise = max(_robust_sigma(sub["y_corrected"].to_numpy(dtype=float)), 1.0)
    weights = 1.0 / np.maximum(np.sqrt(y + noise * noise), noise)
    if len(initial_centers) > 1:
        for midpoint in 0.5 * (initial_centers[:-1] + initial_centers[1:]):
            weights *= 1.0 + 0.45 * np.exp(-0.5 * ((x - float(midpoint)) / 0.010) ** 2)
    weights *= 1.0 + 0.25 * np.exp(-0.5 * ((x - float(initial_centers[0])) / 0.018) ** 2)
    weights *= 1.0 + 0.25 * np.exp(-0.5 * ((x - float(initial_centers[-1])) / 0.018) ** 2)

    model = LinearModel(prefix="b_")
    params = model.make_params(intercept=floor, slope=0.0)
    params["b_intercept"].set(value=floor, min=0.0, max=max(float(np.max(y)) * 0.15, 50.0))
    params["b_slope"].set(value=0.0, min=-max(float(np.max(y)) * 4.0, 200.0), max=max(float(np.max(y)) * 4.0, 200.0))

    for idx, (center, tolerance, prior_area) in enumerate(zip(initial_centers, center_tolerances, initial_areas)):
        prefix = f"c{idx}_"
        component_model = PseudoVoigtModel(prefix=prefix)
        model += component_model

        nearest_idx = int(np.argmin(np.abs(x - center)))
        peak_height = max(float(y[nearest_idx] - floor), max(float(np.max(y)) * 0.01, 1.0))
        rough_area = peak_height * base_sigma * math.sqrt(2.0 * math.pi)
        amplitude_init = float(prior_area) if np.isfinite(prior_area) and float(prior_area) > 0 else rough_area
        amplitude_init = max(amplitude_init, rough_area, 1.0)

        params.update(component_model.make_params())
        params[f"{prefix}center"].set(value=float(center), min=float(center - tolerance), max=float(center + tolerance))
        params[f"{prefix}amplitude"].set(value=float(amplitude_init), min=0.0, max=max(float(np.max(y)) * 40.0, amplitude_init * 6.0, 10.0))
        if idx == 0:
            params[f"{prefix}sigma"].set(value=float(base_sigma), min=float(sigma_min), max=float(sigma_max))
            params[f"{prefix}fraction"].set(value=0.35, min=0.0, max=1.0)
        else:
            params[f"{prefix}sigma"].set(expr="c0_sigma")
            params[f"{prefix}fraction"].set(expr="c0_fraction")

    try:
        result = model.fit(y, params, x=x, weights=weights, nan_policy="omit")
    except Exception:
        return None, {}
    if not getattr(result, "success", False):
        return None, {}

    y_fit = np.asarray(result.best_fit, dtype=float)
    residual = y - y_fit
    ss_res = float(np.sum(residual * residual))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    if not np.isfinite(r2) or r2 < LMFIT_LOCAL_PSEUDOVOIGT_MIN_R2:
        return None, {"r2": r2}

    fitted_components = []
    for idx, _center in enumerate(initial_centers):
        prefix = f"c{idx}_"
        center = float(result.params[f"{prefix}center"].value)
        area = float(max(result.params[f"{prefix}amplitude"].value, 0.0))
        sigma = float(max(result.params[f"{prefix}sigma"].value, 1e-6))
        eta = float(np.clip(result.params[f"{prefix}fraction"].value, 0.0, 1.0))
        fitted_components.append({
            "center": center,
            "area": area,
            "amplitude": area,
            "sigma": sigma,
            "fwhm": 2.0 * sigma,
            "eta": eta,
        })

    return fitted_components, {
        "r2": float(r2),
        "redchi": float(getattr(result, "redchi", np.nan)),
        "weighted": True,
        "shared_sigma": True,
    }


def _fit_cluster_components(
    df: pd.DataFrame,
    initial_centers,
    window_left: float,
    window_right: float,
    center_tolerances,
    initial_areas=None,
    sigma_bounds=(0.003, 0.025),
):
    fit, meta = _fit_cluster_components_split_pseudovoigt(
        df=df,
        initial_centers=initial_centers,
        window_left=window_left,
        window_right=window_right,
        center_tolerances=center_tolerances,
        initial_areas=initial_areas,
        sigma_bounds=sigma_bounds,
    )
    if fit is not None:
        return fit, meta

    fit, meta = _fit_cluster_components_lmfit_pseudovoigt(
        df=df,
        initial_centers=initial_centers,
        window_left=window_left,
        window_right=window_right,
        center_tolerances=center_tolerances,
        initial_areas=initial_areas,
        sigma_bounds=sigma_bounds,
    )
    if fit is not None:
        return fit, meta

    gaussian_fit = _fit_cluster_components_gaussian(
        df=df,
        initial_centers=initial_centers,
        window_left=window_left,
        window_right=window_right,
        center_tolerances=center_tolerances,
        sigma_bounds=sigma_bounds,
    )
    return gaussian_fit, {"r2": np.nan, "redchi": np.nan}


def _refine_cluster_with_deconvolution(
    df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    matched_targets_df: pd.DataFrame,
    cluster_codes,
    default_centers,
    window_left: float,
    window_right: float,
    center_tolerances,
    status: str,
):
    out = matched_targets_df.copy()
    current_rows = out[out["code"].isin(cluster_codes)].copy()
    if len(current_rows) != len(cluster_codes):
        return out, 0.0

    initial_centers = []
    initial_areas = []
    for code, default_center in zip(cluster_codes, default_centers):
        row = current_rows[current_rows["code"] == code]
        found_rt = pd.to_numeric(row["found_rt"], errors="coerce").iloc[0]
        area = pd.to_numeric(row["area"], errors="coerce").iloc[0]
        initial_centers.append(float(found_rt) if np.isfinite(found_rt) else float(default_center))
        initial_areas.append(float(area) if np.isfinite(area) and float(area) > 0 else np.nan)

    fit, fit_meta = _fit_cluster_components(
        df,
        initial_centers=initial_centers,
        window_left=window_left,
        window_right=window_right,
        center_tolerances=center_tolerances,
        initial_areas=initial_areas,
    )
    if fit is None:
        return out, 0.0

    previous_area_sum = float(current_rows["area"].fillna(0.0).sum())
    use_split_pv_boundaries = all("fwhm_left" in component and "fwhm_right" in component and "eta" in component for component in fit)
    if use_split_pv_boundaries:
        boundaries = _derive_split_pseudovoigt_boundaries(
            df=df,
            fitted_components=fit,
            window_left=window_left,
            window_right=window_right,
        )
    elif all("fwhm" in component and "eta" in component for component in fit):
        boundaries = _derive_pseudo_voigt_boundaries(fit, x_left=window_left, x_right=window_right)
    else:
        boundaries = _derive_fit_component_boundaries(fit, window_left=window_left, window_right=window_right)
    if len(boundaries) != len(cluster_codes):
        return out, 0.0

    x_col = _get_x_column_name(df)
    x_all = df[x_col].to_numpy(dtype=float)
    fitted_area_sum = 0.0
    fitted_rows = []
    for component, (start_x, end_x) in zip(fit, boundaries):
        if end_x <= start_x:
            return out, 0.0
        start_idx = int(np.searchsorted(x_all, float(start_x), side="left"))
        end_idx = int(np.searchsorted(x_all, float(end_x), side="right") - 1)
        start_idx = max(0, min(start_idx, len(x_all) - 2))
        end_idx = max(start_idx + 1, min(end_idx, len(x_all) - 1))
        local_x = x_all[start_idx:end_idx + 1]
        if "fwhm_left" in component and "fwhm_right" in component and "eta" in component:
            local_curve = float(component["area"]) * _split_pseudo_voigt_unit_area(
                local_x,
                center=component["center"],
                fwhm_left=component["fwhm_left"],
                fwhm_right=component["fwhm_right"],
                eta=component["eta"],
            )
            assigned_area = float(np.trapezoid(local_curve, local_x))
        elif "fwhm" in component and "eta" in component:
            local_curve = float(component["area"]) * _pseudo_voigt_unit_area(
                local_x,
                center=component["center"],
                fwhm=component["fwhm"],
                eta=component["eta"],
            )
            assigned_area = float(np.trapezoid(local_curve, local_x))
        else:
            assigned_area = float(component["area"])
        fitted_area_sum += assigned_area
        fitted_rows.append((float(component["center"]), float(assigned_area), float(start_x), float(end_x)))

    if previous_area_sum > 0:
        area_ratio = fitted_area_sum / previous_area_sum
        area_ratio_min = SPLIT_PSEUDOVOIGT_AREA_RATIO_MIN if use_split_pv_boundaries else LMFIT_LOCAL_AREA_RATIO_MIN
        area_ratio_max = SPLIT_PSEUDOVOIGT_AREA_RATIO_MAX if use_split_pv_boundaries else LMFIT_LOCAL_AREA_RATIO_MAX
        if area_ratio < area_ratio_min or area_ratio > area_ratio_max:
            return out, 0.0

    for component_idx, (code, default_center, component) in enumerate(zip(cluster_codes, default_centers, fit)):
        row_idx = out.index[out["code"] == code][0]
        component_center, component_area, start_x, end_x = fitted_rows[component_idx]
        out.at[row_idx, "found_rt"] = component_center
        out.at[row_idx, "area"] = component_area
        out.at[row_idx, "status"] = status
        out.at[row_idx, "match_score"] = abs(component_center - float(default_center))
        if np.isfinite(fit_meta.get("r2", np.nan)):
            out.at[row_idx, "fit_r2"] = float(fit_meta["r2"])
        out.at[row_idx, "integration_start_x"] = start_x
        out.at[row_idx, "integration_end_x"] = end_x

    return out, fitted_area_sum - previous_area_sum


def compute_omega_metrics(matched_targets_df: pd.DataFrame):
    result = {
        "omega3_trio": np.nan,
        "omega3_trio_strict": np.nan,
        "omega3_trio_corrected": np.nan,
        "total_area": np.nan,
        "effective_total_area": np.nan,
        "epa_area": 0.0,
        "dha_area": 0.0,
        "dpa_area": 0.0,
        "epa_neighbor_area": 0.0,
        "epa_overlap_credit_area": 0.0,
        "epa_effective_area": 0.0,
        "epa_overlap_fraction": 0.0,
        "epa_overlap_model_applied": False,
        "epa_overlap_extra_scale": 1.0,
        "c22_overlap_source_area": 0.0,
        "c22_overlap_credit_area": 0.0,
        "c22_overlap_fraction": 0.0,
        "c22_overlap_legacy_fraction": 0.0,
        "c22_overlap_model_fraction": np.nan,
        "c22_overlap_model_applied": False,
        "c22_reference_ratio": np.nan,
        "c22_width_scale": 1.0,
        "c18_denominator_scale": 1.0,
    }
    if matched_targets_df is None or matched_targets_df.empty:
        return result
    valid = matched_targets_df[pd.notna(matched_targets_df["area"])].copy()
    if valid.empty:
        return result
    total_area = float(valid["area"].sum())
    if total_area <= 0 or not np.isfinite(total_area):
        return result

    def area_of(code):
        row = valid[valid["code"] == code]
        return float(row["area"].iloc[0]) if not row.empty else 0.0

    def width_of(code):
        row = valid[valid["code"] == code]
        if row.empty:
            return np.nan
        start_x = pd.to_numeric(row["integration_start_x"], errors="coerce").iloc[0] if "integration_start_x" in row else np.nan
        end_x = pd.to_numeric(row["integration_end_x"], errors="coerce").iloc[0] if "integration_end_x" in row else np.nan
        if not (np.isfinite(start_x) and np.isfinite(end_x)):
            return np.nan
        return float(end_x - start_x)

    def status_of(code):
        row = valid[valid["code"] == code]
        if row.empty or "status" not in row:
            return ""
        value = row["status"].iloc[0]
        return "" if pd.isna(value) else str(value)

    epa, dha, dpa = area_of("C20:5"), area_of("C22:6"), area_of("C22:5")
    c20_3 = area_of("C20:3N8")
    c20_4 = area_of("C20:4N6")
    c22_4 = area_of("C22:4")
    c18_2 = area_of("C18:2N6C")
    c18_1 = area_of("C18:1N9C")
    c18_3 = area_of("C18:3N3")
    c18_denominator_scale = 1.0
    effective_total_area = total_area
    if (
        c18_1 > 0
        and c18_2 > 0
        and c18_1 > c18_2 * C18_DENOMINATOR_DOMINANCE_RATIO
        and c18_3 < c18_1 * C18_DENOMINATOR_SMALL_N3_FRACTION
    ):
        c18_denominator_scale = C18_DENOMINATOR_AREA_SCALE
        effective_total_area = total_area - c18_1 * (1.0 - c18_denominator_scale)
    effective_total_area = max(effective_total_area, 1e-9)
    strict_value = 100.0 * (epa + dha + dpa) / effective_total_area
    c18_width = width_of("C18:1N9C")
    if (
        c18_denominator_scale < 0.999
        and c18_2 > 0
        and c18_1 > c18_2 * C18_DENOMINATOR_EXTREME_RATIO
        and c18_3 < c18_1 * C18_DENOMINATOR_EXTREME_SMALL_N3_FRACTION
        and np.isfinite(c18_width)
        and c18_width > C18_DENOMINATOR_EXTREME_WIDTH_MIN
        and strict_value < C18_DENOMINATOR_EXTREME_STRICT_MAX
    ):
        c18_denominator_scale = C18_DENOMINATOR_EXTREME_AREA_SCALE
        effective_total_area = total_area - c18_1 * (1.0 - c18_denominator_scale)
        effective_total_area = max(effective_total_area, 1e-9)
        strict_value = 100.0 * (epa + dha + dpa) / effective_total_area
    epa_credit_area = 0.0
    epa_overlap_fraction = 0.0
    epa_model_applied = False
    epa_extra_scale = 1.0
    epa_to_c20_3_ratio = epa / c20_3 if c20_3 > 0 else np.nan
    if ENABLE_DATA_DRIVEN_C20_EPA_MODEL and c20_3 > 0 and np.isfinite(epa_to_c20_3_ratio) and epa_to_c20_3_ratio < C20_EPA_MODEL_GATE_RATIO_MAX:
        epa_percent = 100.0 * epa / total_area
        c20_3_percent = 100.0 * c20_3 / total_area
        c20_features = np.asarray([
            strict_value,
            math.log(max(epa_to_c20_3_ratio, 1e-6)),
            epa_percent,
            c20_3_percent,
            width_of("C20:5"),
            width_of("C20:3N8"),
            width_of("C20:4N6"),
        ], dtype=float)
        if np.all(np.isfinite(c20_features)):
            z_value = float(C20_EPA_MODEL_PARAMS[0] + np.dot(c20_features / C20_EPA_MODEL_SCALES, C20_EPA_MODEL_PARAMS[1:]))
            epa_overlap_fraction = float(C20_EPA_MODEL_BLEND / (1.0 + math.exp(-z_value)))
            epa_credit_area = float(np.clip(c20_3 * epa_overlap_fraction, 0.0, c20_3))
            epa_model_applied = True
    c20_width_epa = width_of("C20:5")
    c20_width_neighbor = width_of("C20:3N8")
    if (
        epa_credit_area > 0
        and np.isfinite(epa_to_c20_3_ratio)
        and epa_to_c20_3_ratio < 0.50
        and np.isfinite(c20_width_epa)
        and np.isfinite(c20_width_neighbor)
        and c20_width_neighbor > c20_width_epa * C20_EPA_OVERLAP_WIDE_NEIGHBOR_RATIO
    ):
        epa_extra_scale = C20_EPA_OVERLAP_EXTRA_SCALE
        epa_credit_area *= epa_extra_scale
    if (
        epa_credit_area > C20_EPA_UNDERFIT_CREDIT_MIN
        and "matched_c20_fit" in status_of("C20:5")
        and np.isfinite(epa_to_c20_3_ratio)
        and epa_to_c20_3_ratio < C20_EPA_UNDERFIT_RATIO_MAX
        and np.isfinite(c20_width_epa)
        and np.isfinite(c20_width_neighbor)
        and c20_width_neighbor > c20_width_epa * C20_EPA_UNDERFIT_WIDTH_RATIO
        and strict_value < C20_EPA_UNDERFIT_STRICT_MAX
    ):
        epa_extra_scale *= C20_EPA_UNDERFIT_EXTRA_SCALE
        epa_credit_area *= C20_EPA_UNDERFIT_EXTRA_SCALE
    c22_ratio = dpa / c22_4 if c22_4 > 0 else np.nan
    c22_width_values = np.asarray([width_of("C22:6"), width_of("C22:5"), width_of("C22:4")], dtype=float)
    legacy_fraction = 0.0
    if c22_4 > 0 and np.isfinite(c22_ratio) and strict_value >= C22_OVERLAP_TRIGGER_OMEGA_MIN:
        legacy_fraction = float(np.clip(
            C22_OVERLAP_RATIO_OFFSET - C22_OVERLAP_RATIO_SLOPE * c22_ratio,
            0.0,
            C22_OVERLAP_FRACTION_CAP,
        ))
    c22_fraction = legacy_fraction
    model_fraction = np.nan
    model_applied = False
    if ENABLE_DATA_DRIVEN_C22_OVERLAP_MODEL and c22_4 > 0 and legacy_fraction >= C22_OVERLAP_MODEL_APPLY_FRACTION_MIN:
        c22_rows = valid[valid["code"].isin(["C22:6", "C22:5", "C22:4"])].copy()
        c22_status_text = " ".join(c22_rows["status"].fillna("").astype(str).tolist()) if "status" in c22_rows.columns else ""
        if "matched_c22_fit" not in c22_status_text and dpa > 0:
            c22_rows["integration_start_x"] = pd.to_numeric(c22_rows.get("integration_start_x"), errors="coerce")
            c22_rows["integration_end_x"] = pd.to_numeric(c22_rows.get("integration_end_x"), errors="coerce")
            width_by_code = (
                c22_rows.assign(width=lambda frame: frame["integration_end_x"] - frame["integration_start_x"])
                .set_index("code")["width"]
                .to_dict()
            )
            dpa_percent = 100.0 * dpa / total_area
            c22_4_percent = 100.0 * c22_4 / total_area
            dha_dpa_ratio = dha / dpa if dpa > 0 else np.nan
            features = np.asarray([
                strict_value,
                math.log(max(c22_ratio, 1e-6)),
                math.log(max(dha_dpa_ratio, 1e-6)) if np.isfinite(dha_dpa_ratio) else 0.0,
                float(width_by_code.get("C22:6", np.nan)),
                float(width_by_code.get("C22:5", np.nan)),
                float(width_by_code.get("C22:4", np.nan)),
                c22_4_percent,
                dpa_percent,
            ], dtype=float)
            if np.all(np.isfinite(features)):
                scaled = features / C22_OVERLAP_MODEL_SCALES
                z_value = float(C22_OVERLAP_MODEL_PARAMS[0] + np.dot(scaled, C22_OVERLAP_MODEL_PARAMS[1:]))
                model_fraction = float(C22_OVERLAP_FRACTION_CAP / (1.0 + math.exp(-z_value)))
                c22_fraction = float(np.clip(
                    C22_OVERLAP_MODEL_BLEND * model_fraction + (1.0 - C22_OVERLAP_MODEL_BLEND) * legacy_fraction,
                    0.0,
                    C22_OVERLAP_FRACTION_CAP,
                ))
                model_applied = True
    c22_width_scale = 1.0
    finite_c22_widths = c22_width_values[np.isfinite(c22_width_values)]
    if finite_c22_widths.size == 3 and float(np.mean(finite_c22_widths)) > C22_OVERLAP_WIDE_CLUSTER_MEAN_WIDTH:
        c22_width_scale = C22_OVERLAP_WIDE_CLUSTER_SCALE
        c22_fraction *= c22_width_scale
    c22_credit_area = c22_4 * c22_fraction
    corrected_value = 100.0 * (epa + dha + dpa + epa_credit_area + c22_credit_area) / effective_total_area
    result.update({
        "omega3_trio": corrected_value,
        "omega3_trio_strict": strict_value,
        "omega3_trio_corrected": corrected_value,
        "total_area": total_area,
        "effective_total_area": effective_total_area,
        "epa_area": epa,
        "dha_area": dha,
        "dpa_area": dpa,
        "epa_neighbor_area": c20_3,
        "epa_overlap_credit_area": epa_credit_area,
        "epa_effective_area": epa + epa_credit_area,
        "epa_overlap_fraction": epa_overlap_fraction,
        "epa_overlap_model_applied": epa_model_applied,
        "epa_overlap_extra_scale": epa_extra_scale,
        "c22_overlap_source_area": c22_4,
        "c22_overlap_credit_area": c22_credit_area,
        "c22_overlap_fraction": c22_fraction,
        "c22_overlap_legacy_fraction": legacy_fraction,
        "c22_overlap_model_fraction": model_fraction,
        "c22_overlap_model_applied": model_applied,
        "c22_reference_ratio": c22_ratio,
        "c22_width_scale": c22_width_scale,
        "c18_denominator_scale": c18_denominator_scale,
    })
    return result


def build_confidence_assessment(
    matched_targets_df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    omega: dict,
    baseline_mode: str = "chebyshev",
    cluster_quality_score: float | None = None,
) -> dict:
    result = {
        "score": np.nan,
        "level": "unknown",
        "label": "—",
        "button_text": "Уверенность: —",
        "reasons": [],
        "metrics": [],
    }
    if matched_targets_df is None or matched_targets_df.empty or not isinstance(omega, dict):
        return result

    valid = matched_targets_df.copy()
    valid["area"] = pd.to_numeric(valid.get("area"), errors="coerce")
    valid["found_rt"] = pd.to_numeric(valid.get("found_rt"), errors="coerce")
    valid["matched_peak_id"] = pd.to_numeric(valid.get("matched_peak_id"), errors="coerce")
    valid["integration_start_x"] = pd.to_numeric(valid.get("integration_start_x"), errors="coerce")
    valid["integration_end_x"] = pd.to_numeric(valid.get("integration_end_x"), errors="coerce")

    score = 100.0
    reason_items: list[tuple[float, str]] = []

    def penalize(points: float, reason: str):
        nonlocal score
        score -= float(points)
        reason_items.append((float(points), reason))

    def area_of(code: str) -> float:
        row = valid[valid["code"] == code]
        if row.empty:
            return 0.0
        value = row["area"].iloc[0]
        return float(value) if pd.notna(value) else 0.0

    def width_of(code: str) -> float:
        row = valid[valid["code"] == code]
        if row.empty:
            return np.nan
        start_x = row["integration_start_x"].iloc[0]
        end_x = row["integration_end_x"].iloc[0]
        if not (np.isfinite(start_x) and np.isfinite(end_x)):
            return np.nan
        return float(end_x - start_x)

    def status_of(code: str) -> str:
        row = valid[valid["code"] == code]
        if row.empty or "status" not in row:
            return ""
        value = row["status"].iloc[0]
        return "" if pd.isna(value) else str(value)

    omega_value = float(omega.get("omega3_trio", np.nan))
    strict_value = float(omega.get("omega3_trio_strict", np.nan))
    spread = abs(omega_value - strict_value) if np.isfinite(omega_value) and np.isfinite(strict_value) else np.nan

    if baseline_mode != "chebyshev":
        penalize(8.0, "Потребовался fallback baseline")

    if cluster_quality_score is not None and np.isfinite(cluster_quality_score) and cluster_quality_score < CLUSTER_QUALITY_COMPLETE_SCORE:
        gap = float(CLUSTER_QUALITY_COMPLETE_SCORE - cluster_quality_score)
        penalize(min(18.0, 6.0 + gap * 0.45), f"Качество кластеров ниже целевого ({cluster_quality_score:.1f})")

    matched_count = int(valid["matched_peak_id"].notna().sum())
    missing_count = int(max(0, len(valid) - matched_count))
    if missing_count > 0:
        penalize(min(18.0, 4.0 * missing_count), f"Есть неполные матчинг-пики ({missing_count})")

    trio_codes = ["C20:5", "C22:6", "C22:5"]
    trio_missing = [code for code in trio_codes if area_of(code) <= 0]
    if trio_missing:
        penalize(30.0, f"Отсутствуют ключевые omega-3 пики: {', '.join(trio_missing)}")

    if np.isfinite(spread):
        if spread > 0.45:
            penalize(18.0, f"Сильный разброс strict/final omega ({spread:.2f})")
        elif spread > 0.25:
            penalize(10.0, f"Заметный разброс strict/final omega ({spread:.2f})")
        elif spread > 0.12:
            penalize(5.0, f"Небольшой разброс strict/final omega ({spread:.2f})")

    c18_scale = float(omega.get("c18_denominator_scale", 1.0))
    c18_1 = area_of("C18:1N9C")
    c18_2 = area_of("C18:2N6C")
    c18_3 = area_of("C18:3N3")
    c18_width = width_of("C18:1N9C")
    c18_ratio = c18_1 / c18_2 if c18_2 > 0 else np.nan
    c18_n3_fraction = c18_3 / c18_1 if c18_1 > 0 else np.nan
    c18_status = status_of("C18:1N9C")
    if c18_scale < 0.999:
        penalize(8.0 if c18_scale >= 0.80 else 12.0, f"Сработала коррекция denominator для C18 ({c18_scale:.2f})")
    if (
        np.isfinite(c18_ratio)
        and np.isfinite(c18_n3_fraction)
        and np.isfinite(c18_width)
        and c18_ratio > C18_DENOMINATOR_EXTREME_RATIO
        and c18_n3_fraction < C18_DENOMINATOR_EXTREME_SMALL_N3_FRACTION
        and c18_width > C18_DENOMINATOR_EXTREME_WIDTH_MIN
    ):
        penalize(10.0, "C18:1N9C выглядит перерастянутым относительно C18-кластера")
    if "matched_c18_local_bounds" in c18_status or "matched_c18_pvfit" in c18_status:
        penalize(6.0, "C18-кластер потребовал локальную переинтеграцию")

    epa = area_of("C20:5")
    c20_3 = area_of("C20:3N8")
    epa_status = status_of("C20:5")
    epa_credit = float(omega.get("epa_overlap_credit_area", 0.0))
    epa_ratio = epa / c20_3 if c20_3 > 0 else np.nan
    w_epa = width_of("C20:5")
    w_c20_3 = width_of("C20:3N8")
    if "matched_c20_fit" in epa_status:
        penalize(10.0, "EPA восстанавливался через C20 fit")
    elif "matched_c20_local" in epa_status:
        penalize(5.0, "EPA потребовал локальную C20-коррекцию")
    if epa_credit > C20_EPA_UNDERFIT_CREDIT_MIN:
        penalize(5.0, f"Существенный overlap-credit у EPA ({epa_credit:.1f})")
    if (
        np.isfinite(epa_ratio)
        and np.isfinite(w_epa)
        and np.isfinite(w_c20_3)
        and epa_ratio < C20_EPA_UNDERFIT_RATIO_MAX
        and w_c20_3 > w_epa * C20_EPA_UNDERFIT_WIDTH_RATIO
    ):
        penalize(10.0, "EPA выглядит недобранным относительно C20:3N8")

    c22_statuses = [status_of(code) for code in ["C22:6", "C22:5", "C22:4"]]
    c22_status_text = " ".join(status for status in c22_statuses if status)
    c22_credit = float(omega.get("c22_overlap_credit_area", 0.0))
    c22_widths = np.asarray([width_of("C22:6"), width_of("C22:5"), width_of("C22:4")], dtype=float)
    c22_mean_width = float(np.nanmean(c22_widths)) if np.isfinite(c22_widths).any() else np.nan
    if "matched_c22_pvfit" in c22_status_text:
        penalize(12.0, "C22-кластер потребовал pvfit refinement")
    elif "matched_c22_fit" in c22_status_text:
        penalize(10.0, "C22-кластер потребовал fit-восстановление")
    elif "tailtight" in c22_status_text:
        penalize(6.0, "C22-кластер потребовал tail tightening")
    if c22_credit > 0:
        penalize(4.0 if c22_credit < 80 else 8.0, f"C22 overlap-credit участвует в расчёте ({c22_credit:.1f})")
    if np.isfinite(c22_mean_width):
        if c22_mean_width > 0.036:
            penalize(12.0, f"C22-кластер всё ещё широкий ({c22_mean_width:.3f} min)")
        elif c22_mean_width > 0.032:
            penalize(6.0, f"C22-кластер умеренно широкий ({c22_mean_width:.3f} min)")

    peak_count = len(peaks_df) if peaks_df is not None else 0
    if peak_count >= 65:
        penalize(5.0, f"Необычно много детектированных пиков ({peak_count})")

    score = max(0.0, min(100.0, score))
    if score >= 85.0:
        level = "Авто OK"
    elif score >= 70.0:
        level = "Быстро проверить"
    elif score >= 55.0:
        level = "Проверить руками"
    else:
        level = "Ручная проверка"

    reason_items.sort(key=lambda item: item[0], reverse=True)
    reasons = [f"-{int(round(points))}: {text}" for points, text in reason_items]
    metrics = []
    if np.isfinite(strict_value) and np.isfinite(omega_value):
        metrics.append(f"Omega final / strict: {omega_value:.2f}% / {strict_value:.2f}%")
    if np.isfinite(spread):
        metrics.append(f"Разброс strict/final: {spread:.2f}")
    metrics.append(f"Baseline: {baseline_mode}")
    if cluster_quality_score is not None and np.isfinite(cluster_quality_score):
        metrics.append(f"Cluster quality: {cluster_quality_score:.1f}")
    metrics.append(f"Matched peaks: {matched_count}/{len(valid)}")

    result["score"] = score
    result["level"] = level
    result["label"] = level
    result["button_text"] = f"Уверенность: {int(round(score))}"
    result["reasons"] = reasons
    result["metrics"] = metrics
    return result


def refine_overlapped_c22_cluster_areas(
    df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    matched_targets_df: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets_df.copy()
    if df is None or df.empty or peaks_df is None or peaks_df.empty or out is None or out.empty:
        return out

    c22_codes = ["C22:6", "C22:5", "C22:4"]
    cluster_rows = out[out["code"].isin(c22_codes)].copy()
    if len(cluster_rows) != len(c22_codes):
        return out
    if cluster_rows["matched_peak_id"].isna().any() or cluster_rows["found_rt"].isna().any():
        return out

    ordered_rows = cluster_rows.set_index("code").loc[c22_codes].reset_index()
    peak_rows = []
    for peak_id in ordered_rows["matched_peak_id"]:
        matched_peak = peaks_df[peaks_df["peak_id"] == int(peak_id)]
        if matched_peak.empty:
            return out
        peak_rows.append(matched_peak.iloc[0])

    first_peak = peak_rows[0]
    third_peak = peak_rows[2]
    first_peak_width = float(first_peak["end_x"] - first_peak["start_x"])
    if float(first_peak["end_x"]) <= float(third_peak["apex_x"]) or first_peak_width < 0.06:
        return out

    x_col = _get_x_column_name(df)
    x = df[x_col].to_numpy(dtype=float)
    y_corrected_raw = df["y_corrected"].to_numpy(dtype=float)
    y_smooth = df["y_smooth"].to_numpy(dtype=float)
    apex_indices = [int(np.argmin(np.abs(x - float(rt)))) for rt in ordered_rows["found_rt"]]

    left_limit = max(0, apex_indices[0] - 80)
    right_limit = min(len(x) - 1, apex_indices[2] + 80)
    local_metric_pack = _build_cluster_local_metric(
        x=x,
        y_corrected=y_corrected_raw,
        y_smooth=y_smooth,
        start_idx=left_limit,
        end_idx=right_limit,
    )
    if local_metric_pack is None:
        return out
    _, _, _, boundary_metric = local_metric_pack
    y_corrected = np.clip(y_corrected_raw, 0.0, None)

    left_boundary = _find_preferred_minimum_index(
        boundary_metric,
        left_limit,
        apex_indices[0],
        target_idx=left_limit + 0.25 * max(apex_indices[0] - left_limit, 1),
    )
    split_1 = _find_preferred_minimum_index(
        boundary_metric,
        apex_indices[0],
        apex_indices[1],
        target_idx=0.5 * (apex_indices[0] + apex_indices[1]),
    )
    split_2 = _find_preferred_minimum_index(
        boundary_metric,
        apex_indices[1],
        apex_indices[2],
        target_idx=0.5 * (apex_indices[1] + apex_indices[2]),
    )
    right_boundary = _find_preferred_minimum_index(
        boundary_metric,
        apex_indices[2],
        right_limit,
        target_idx=apex_indices[2] + 0.75 * max(right_limit - apex_indices[2], 1),
    )
    boundaries = [left_boundary, split_1, split_2, right_boundary]

    if any(boundaries[i] >= boundaries[i + 1] for i in range(len(boundaries) - 1)):
        return out

    for i, code in enumerate(c22_codes):
        start_idx = boundaries[i]
        end_idx = boundaries[i + 1]
        area = float(np.trapezoid(y_corrected[start_idx:end_idx + 1], x[start_idx:end_idx + 1]))
        row_idx = out.index[out["code"] == code][0]
        out.at[row_idx, "area"] = area
        out.at[row_idx, "percent_area"] = np.nan
        out.at[row_idx, "status"] = f"{out.at[row_idx, 'status']}_split"
        out.at[row_idx, "integration_start_x"] = float(x[start_idx])
        out.at[row_idx, "integration_end_x"] = float(x[end_idx])

    total_area = float(pd.to_numeric(out["area"], errors="coerce").fillna(0.0).sum())
    if total_area > 0:
        out["percent_area"] = 100.0 * pd.to_numeric(out["area"], errors="coerce") / total_area
    return out


def recover_missing_c22_components_with_fit(
    df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    matched_targets_df: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets_df.copy()
    if df is None or df.empty or peaks_df is None or peaks_df.empty or out is None or out.empty:
        return out

    c22_codes = ["C22:6", "C22:5", "C22:4"]
    cluster = out[out["code"].isin(c22_codes)].copy()
    if len(cluster) != len(c22_codes):
        return out

    area_by_code = cluster.set_index("code")["area"].apply(lambda value: pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0])
    found_count = int(area_by_code.notna().sum())
    dpa_missing = bool(pd.isna(area_by_code.get("C22:5")))
    if not dpa_missing or found_count > 2:
        return out
    current_cluster_total = float(area_by_code.fillna(0.0).sum())

    fit_out, _ = _refine_cluster_with_deconvolution(
        df=df,
        peaks_df=peaks_df,
        matched_targets_df=out,
        cluster_codes=c22_codes,
        default_centers=[9.247, 9.280, 9.310],
        window_left=9.22,
        window_right=9.33,
        center_tolerances=[0.010, 0.010, 0.010],
        status="matched_c22_fit",
    )
    fitted_cluster = fit_out[fit_out["code"].isin(c22_codes)].copy()
    fitted_cluster_total = float(pd.to_numeric(fitted_cluster["area"], errors="coerce").fillna(0.0).sum())
    if current_cluster_total > 0 and fitted_cluster_total < current_cluster_total * 0.995:
        return out
    return _recompute_matched_percent_area(fit_out)


def recover_underintegrated_c20_components_with_fit(
    df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    matched_targets_df: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets_df.copy()
    if df is None or df.empty or peaks_df is None or peaks_df.empty or out is None or out.empty:
        return out

    epa_row = out[out["code"] == "C20:5"]
    if epa_row.empty:
        return out

    epa_area = pd.to_numeric(epa_row["area"], errors="coerce").iloc[0]
    matched_peak_id = pd.to_numeric(epa_row["matched_peak_id"], errors="coerce").iloc[0]
    if not np.isfinite(epa_area) or not np.isfinite(matched_peak_id):
        return out
    if float(epa_area) > C20_FIT_EPA_AREA_MAX:
        return out

    peak_row = peaks_df[peaks_df["peak_id"] == int(matched_peak_id)]
    if peak_row.empty:
        return out
    epa_prominence = float(peak_row.iloc[0].get("raw_prominence", peak_row.iloc[0]["prominence"]))
    if epa_prominence > C20_FIT_EPA_PROMINENCE_MAX:
        return out

    fit_out, _ = _refine_cluster_with_deconvolution(
        df=df,
        peaks_df=peaks_df,
        matched_targets_df=out,
        cluster_codes=["C20:4N6", "C20:5", "C20:3N8"],
        default_centers=[8.382, 8.410, 8.467],
        window_left=8.35,
        window_right=8.50,
        center_tolerances=[0.010, 0.010, 0.015],
        status="matched_c20_fit",
    )
    return _recompute_matched_percent_area(fit_out)


def recover_overlapped_c18_components_with_fit(
    df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    matched_targets_df: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets_df.copy()
    if (
        not ENABLE_LMFIT_C18_RECOVERY
        or df is None or df.empty
        or peaks_df is None or peaks_df.empty
        or out is None or out.empty
    ):
        return out

    c18_codes = ["C18:1N9C", "C18:3N3", "C18:0"]
    if not (
        _should_force_c18_valley_split(out)
        or _cluster_has_duplicate_peak_ids(out, c18_codes)
        or _cluster_has_integration_overlap(out, c18_codes)
    ):
        return out

    fit_out, _ = _refine_cluster_with_deconvolution(
        df=df,
        peaks_df=peaks_df,
        matched_targets_df=out,
        cluster_codes=c18_codes,
        default_centers=[7.623, 7.650, 7.750],
        window_left=7.57,
        window_right=7.79,
        center_tolerances=[0.015, 0.015, 0.018],
        status="matched_c18_pvfit",
    )
    return _recompute_matched_percent_area(fit_out)


class ChromatogramApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Chromatogram Peak Detector (Integrated)")
        self.root.geometry("1780x1040")

        self.reference_json_path = ensure_runtime_file("reference_targets_reverted_c22fixed.json")
        self.reference_targets = omega_chromatopy_clean.load_reference_targets(self.reference_json_path)

        self.current_file = None
        self.current_sample_name = ""
        self.loaded_batches = []
        self.current_batch_index = 0
        self.batch_tree = None
        self.batch_results_window = None
        self.batch_results_tree = None
        self._batch_tree_syncing = False
        self.df_processed = None
        self.best_window = None
        self.peaks_df = pd.DataFrame()
        self.matched_targets_df = pd.DataFrame()
        self.current_rt_shift = 0.0
        self.selected_target_code = None
        self.manual_start_var = tk.StringVar(value="")
        self.manual_end_var = tk.StringVar(value="")
        self._manual_drag_active_boundary = None

        self.status_var = tk.StringVar(value="Выбери CSV-файл.")
        self.file_var = tk.StringVar(value="Файл не выбран")
        self.omega_var = tk.StringVar(value="Omega-3: —")
        self.integration_var = tk.StringVar(value="Integration: —")
        self.gamma_var = tk.StringVar(value="γ-Linolenic: —")
        self.batch_var = tk.StringVar(value="Series: —")
        self.confidence_var = tk.StringVar(value="Уверенность: —")
        self.current_confidence = None

        self._build_ui()

    def _build_ui(self):
        controls = ttk.Frame(self.root, padding=10)
        controls.pack(fill="x")
        self.confidence_button = ttk.Button(
            controls,
            textvariable=self.confidence_var,
            command=self.show_confidence_details,
            width=22,
        )
        self.confidence_button.pack(side="right")
        self.confidence_button.state(["disabled"])
        ttk.Button(controls, text="Открыть CSV", command=self.open_file).pack(side="left", padx=(0, 10))
        self.prev_button = ttk.Button(controls, text="←", width=4, command=self.prev_batch)
        self.prev_button.pack(side="left", padx=(0, 4))
        self.next_button = ttk.Button(controls, text="→", width=4, command=self.next_batch)
        self.next_button.pack(side="left", padx=(0, 10))
        self.batch_results_button = ttk.Button(controls, text="Результаты Batch", command=self.open_batch_results_window)
        self.batch_results_button.pack(side="left", padx=(0, 10))
        ttk.Label(controls, textvariable=self.batch_var).pack(side="left", padx=(0, 14))
        ttk.Label(controls, textvariable=self.file_var).pack(side="left", padx=(10, 20))
        ttk.Label(controls, textvariable=self.omega_var).pack(side="left")
        ttk.Label(controls, textvariable=self.integration_var).pack(side="left", padx=(20, 0))
        ttk.Label(controls, textvariable=self.gamma_var).pack(side="left", padx=(20, 0))

        content = ttk.Frame(self.root, padding=(10, 0, 10, 0))
        content.pack(fill="both", expand=True)

        plot_frame = ttk.Frame(content)
        plot_frame.pack(side="left", fill="both", expand=True)
        sidebar = ttk.Frame(content, width=560)
        sidebar.pack(side="right", fill="y", padx=(12, 0))
        sidebar.pack_propagate(False)

        self.figure = Figure(figsize=(12, 8), dpi=100)
        self.figure.subplots_adjust(left=0.055, right=0.985, top=0.955, bottom=0.065)
        grid = self.figure.add_gridspec(3, 2, height_ratios=[2.5, 1.0, 1.0], hspace=0.26, wspace=0.18)
        self.ax = self.figure.add_subplot(grid[0, :])
        self.preview_axes = []
        self.preview_specs = PREVIEW_WINDOWS[:]
        for index, spec in enumerate(self.preview_specs):
            row = 1 + index // 2
            col = index % 2
            preview_ax = self.figure.add_subplot(grid[row, col])
            self.preview_axes.append(preview_ax)

        self.canvas = FigureCanvasTkAgg(self.figure, master=plot_frame)
        self.canvas.mpl_connect("button_press_event", self.handle_manual_boundary_press)
        self.canvas.mpl_connect("motion_notify_event", self.handle_manual_boundary_motion)
        self.canvas.mpl_connect("button_release_event", self.handle_manual_boundary_release)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar_host = tk.Frame(plot_frame)
        toolbar_host.pack(fill="x")
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_host, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side="left", fill="x")

        batch_frame = ttk.LabelFrame(sidebar, text="Batch", padding=(8, 8))
        batch_frame.pack(fill="x")
        batch_columns = ("sample_name", "omega_value")
        batch_tree_frame = ttk.Frame(batch_frame)
        batch_tree_frame.pack(fill="both", expand=True)
        self.batch_tree = ttk.Treeview(batch_tree_frame, columns=batch_columns, show="headings", height=14, selectmode="browse")
        self.batch_tree.heading("sample_name", text="Образец")
        self.batch_tree.heading("omega_value", text="Omega-3")
        self.batch_tree.column("sample_name", width=300, anchor="w")
        self.batch_tree.column("omega_value", width=110, anchor="center")
        self.batch_tree.pack(side="left", fill="both", expand=True)
        batch_scroll = ttk.Scrollbar(batch_tree_frame, orient="vertical", command=self.batch_tree.yview)
        batch_scroll.pack(side="right", fill="y")
        self.batch_tree.configure(yscrollcommand=batch_scroll.set)
        self.batch_tree.bind("<<TreeviewSelect>>", self.handle_batch_tree_selection)

        table_frame = ttk.LabelFrame(sidebar, text="Пики", padding=(8, 8))
        table_frame.pack(fill="both", expand=True, pady=(10, 0))
        cols = ["display_name", "code", "expected_rt", "found_rt", "area", "percent_area", "status"]
        peaks_tree_frame = ttk.Frame(table_frame)
        peaks_tree_frame.pack(fill="both", expand=True)
        self.tree = ttk.Treeview(peaks_tree_frame, columns=cols, show="headings", height=18)
        for c in cols:
            self.tree.heading(c, text=c)
            width = 118 if c not in {"display_name", "status"} else 160
            anchor = "w" if c in {"display_name", "status"} else "center"
            self.tree.column(c, width=width, anchor=anchor)
        self.tree.pack(side="left", fill="both", expand=True)
        peaks_scroll = ttk.Scrollbar(peaks_tree_frame, orient="vertical", command=self.tree.yview)
        peaks_scroll.pack(side="right", fill="y")
        self.tree.configure(yscrollcommand=peaks_scroll.set)
        self.tree.bind("<<TreeviewSelect>>", self.handle_target_selection)

        manual_frame = ttk.LabelFrame(sidebar, text="Ручная интеграция", padding=(8, 8))
        manual_frame.pack(fill="x", pady=(10, 0))
        ttk.Label(manual_frame, text="Start RT").grid(row=0, column=0, sticky="w")
        ttk.Entry(manual_frame, textvariable=self.manual_start_var, width=12).grid(row=0, column=1, sticky="ew", padx=(6, 10))
        ttk.Label(manual_frame, text="End RT").grid(row=0, column=2, sticky="w")
        ttk.Entry(manual_frame, textvariable=self.manual_end_var, width=12).grid(row=0, column=3, sticky="ew", padx=(6, 0))
        ttk.Button(manual_frame, text="Взять текущие", command=self.load_selected_integration_bounds).grid(row=1, column=0, columnspan=2, sticky="ew", pady=(8, 0), padx=(0, 6))
        ttk.Button(manual_frame, text="Применить", command=self.apply_manual_integration).grid(row=1, column=2, columnspan=2, sticky="ew", pady=(8, 0))
        manual_frame.columnconfigure(1, weight=1)
        manual_frame.columnconfigure(3, weight=1)

        status = ttk.Label(self.root, textvariable=self.status_var, anchor="w", padding=(10, 4))
        status.pack(fill="x")
        self.update_batch_navigation()

    def update_batch_navigation(self):
        total = len(self.loaded_batches)
        if total <= 0:
            self.batch_var.set("Series: —")
            self.prev_button.state(["disabled"])
            self.next_button.state(["disabled"])
            self.batch_results_button.state(["disabled"])
            self.confidence_var.set("Уверенность: —")
            self.current_confidence = None
            if self.confidence_button is not None:
                self.confidence_button.state(["disabled"])
            if self.batch_tree is not None:
                for item_id in self.batch_tree.get_children():
                    self.batch_tree.delete(item_id)
            return

        self.batch_var.set(f"Series: {self.current_batch_index + 1}/{total}")
        self.batch_results_button.state(["!disabled"])
        if self.current_batch_index > 0:
            self.prev_button.state(["!disabled"])
        else:
            self.prev_button.state(["disabled"])

        if self.current_batch_index < total - 1:
            self.next_button.state(["!disabled"])
        else:
            self.next_button.state(["disabled"])
        self.populate_main_batch_tree()

    def show_confidence_details(self):
        confidence = self.current_confidence
        if not confidence or not np.isfinite(confidence.get("score", np.nan)):
            messagebox.showinfo("Уверенность", "Нет данных для оценки уверенности.", parent=self.root)
            return

        lines = [
            f"Уверенность: {int(round(confidence['score']))}/100",
            f"Статус: {confidence.get('label', '—')}",
            "",
        ]
        reasons = confidence.get("reasons") or []
        if reasons:
            lines.append("Причины снижения уверенности:")
            lines.extend(reasons)
        else:
            lines.append("Сильных причин для ручной проверки не найдено.")

        metrics = confidence.get("metrics") or []
        if metrics:
            lines.append("")
            lines.append("Контекст:")
            lines.extend(metrics)

        messagebox.showinfo("Уверенность", "\n".join(lines), parent=self.root)

    def handle_target_selection(self, event=None):
        selection = self.tree.selection()
        self.selected_target_code = selection[0] if selection else None
        self.load_selected_integration_bounds(silent=True)
        if self.df_processed is not None:
            self.update_plot()
            return

    def load_selected_integration_bounds(self, silent: bool = False):
        if not self.selected_target_code or self.matched_targets_df.empty:
            self.manual_start_var.set("")
            self.manual_end_var.set("")
            if not silent:
                self.status_var.set("Выбери пик в таблице перед ручной интеграцией.")
            return

        row = self.matched_targets_df[self.matched_targets_df["code"] == self.selected_target_code]
        if row.empty:
            self.manual_start_var.set("")
            self.manual_end_var.set("")
            return
        start_x = pd.to_numeric(row["integration_start_x"], errors="coerce").iloc[0]
        end_x = pd.to_numeric(row["integration_end_x"], errors="coerce").iloc[0]
        self.manual_start_var.set("" if not np.isfinite(start_x) else f"{float(start_x):.5f}")
        self.manual_end_var.set("" if not np.isfinite(end_x) else f"{float(end_x):.5f}")
        if not silent:
            self.status_var.set(f"Границы {self.selected_target_code} загружены для ручной правки.")

    def apply_manual_integration(self):
        if self.df_processed is None or self.matched_targets_df.empty:
            messagebox.showwarning("Ручная интеграция", "Сначала открой CSV и выбери образец.", parent=self.root)
            return
        if not self.selected_target_code:
            messagebox.showwarning("Ручная интеграция", "Сначала выбери пик в таблице справа.", parent=self.root)
            return

        try:
            start_x = float(str(self.manual_start_var.get()).replace(",", "."))
            end_x = float(str(self.manual_end_var.get()).replace(",", "."))
        except ValueError:
            messagebox.showerror("Ручная интеграция", "Start RT и End RT должны быть числами.", parent=self.root)
            return
        if not (np.isfinite(start_x) and np.isfinite(end_x) and end_x > start_x):
            messagebox.showerror("Ручная интеграция", "End RT должен быть больше Start RT.", parent=self.root)
            return

        x_col = _get_x_column_name(self.df_processed)
        x = self.df_processed[x_col].to_numpy(dtype=float)
        y = self.df_processed["y_corrected"].to_numpy(dtype=float)
        if start_x < float(np.nanmin(x)) or end_x > float(np.nanmax(x)):
            messagebox.showerror("Ручная интеграция", "Границы вне диапазона текущей хроматограммы.", parent=self.root)
            return

        start_idx = int(np.searchsorted(x, start_x, side="left"))
        end_idx = int(np.searchsorted(x, end_x, side="right") - 1)
        start_idx = max(0, min(start_idx, len(x) - 2))
        end_idx = max(start_idx + 1, min(end_idx, len(x) - 1))
        segment_y = np.clip(y[start_idx:end_idx + 1], 0.0, None)
        segment_x = x[start_idx:end_idx + 1]
        area = float(np.trapezoid(segment_y, segment_x))
        apex_idx = int(start_idx + np.argmax(segment_y))

        row_mask = self.matched_targets_df["code"] == self.selected_target_code
        if not row_mask.any():
            messagebox.showerror("Ручная интеграция", f"Пик {self.selected_target_code} не найден в таблице.", parent=self.root)
            return
        row_idx = self.matched_targets_df.index[row_mask][0]
        old_status = str(self.matched_targets_df.at[row_idx, "status"] or "")
        manual_status = old_status if "manual" in old_status else f"{old_status}_manual" if old_status else "manual"
        self.matched_targets_df.at[row_idx, "integration_start_x"] = float(x[start_idx])
        self.matched_targets_df.at[row_idx, "integration_end_x"] = float(x[end_idx])
        self.matched_targets_df.at[row_idx, "found_rt"] = float(x[apex_idx])
        self.matched_targets_df.at[row_idx, "area"] = area
        self.matched_targets_df.at[row_idx, "matched_peak_id"] = np.nan
        self.matched_targets_df.at[row_idx, "match_score"] = np.nan
        self.matched_targets_df.at[row_idx, "status"] = manual_status
        self.matched_targets_df = _recompute_matched_percent_area(self.matched_targets_df)
        self.selected_target_code = str(self.selected_target_code)
        self.refresh_peaks()
        if self.selected_target_code in self.tree.get_children():
            self.tree.selection_set(self.selected_target_code)
            self.tree.focus(self.selected_target_code)
        self.status_var.set(
            f"Ручная интеграция {self.selected_target_code}: {x[start_idx]:.5f}–{x[end_idx]:.5f}, area {area:.4f}"
        )

    def _manual_bounds_from_vars(self):
        try:
            start_x = float(str(self.manual_start_var.get()).replace(",", "."))
            end_x = float(str(self.manual_end_var.get()).replace(",", "."))
        except ValueError:
            return np.nan, np.nan
        return start_x, end_x

    def _selected_manual_drag_bounds(self):
        if not self.selected_target_code or self.matched_targets_df.empty:
            return np.nan, np.nan
        start_x, end_x = self._manual_bounds_from_vars()
        if np.isfinite(start_x) and np.isfinite(end_x):
            return start_x, end_x
        row = self.matched_targets_df[self.matched_targets_df["code"] == self.selected_target_code]
        if row.empty:
            return np.nan, np.nan
        start_x = pd.to_numeric(row["integration_start_x"], errors="coerce").iloc[0]
        end_x = pd.to_numeric(row["integration_end_x"], errors="coerce").iloc[0]
        return float(start_x), float(end_x)

    def handle_manual_boundary_press(self, event):
        if event.button != 1 or event.inaxes not in [self.ax, *getattr(self, "preview_axes", [])]:
            return
        if event.xdata is None or self.df_processed is None or not self.selected_target_code:
            return
        start_x, end_x = self._selected_manual_drag_bounds()
        if not (np.isfinite(start_x) and np.isfinite(end_x)):
            return

        x_min, x_max = event.inaxes.get_xlim()
        tolerance = max(0.0025, abs(float(x_max) - float(x_min)) * 0.010)
        distances = {"start": abs(float(event.xdata) - start_x), "end": abs(float(event.xdata) - end_x)}
        boundary, distance = min(distances.items(), key=lambda item: item[1])
        if distance > tolerance:
            return
        self._manual_drag_active_boundary = boundary
        self.status_var.set(f"Тяни {'левую' if boundary == 'start' else 'правую'} границу {self.selected_target_code} мышкой…")

    def handle_manual_boundary_motion(self, event):
        if self._manual_drag_active_boundary is None or event.xdata is None or self.df_processed is None:
            return
        x_col = _get_x_column_name(self.df_processed)
        x_values = self.df_processed[x_col].to_numpy(dtype=float)
        x_value = float(np.clip(float(event.xdata), float(np.nanmin(x_values)), float(np.nanmax(x_values))))
        start_x, end_x = self._selected_manual_drag_bounds()
        if self._manual_drag_active_boundary == "start":
            if np.isfinite(end_x):
                x_value = min(x_value, float(end_x) - 1e-5)
            self.manual_start_var.set(f"{x_value:.5f}")
        else:
            if np.isfinite(start_x):
                x_value = max(x_value, float(start_x) + 1e-5)
            self.manual_end_var.set(f"{x_value:.5f}")
        self.update_plot()

    def handle_manual_boundary_release(self, event):
        if self._manual_drag_active_boundary is None:
            return
        self._manual_drag_active_boundary = None
        self.status_var.set(f"Граница {self.selected_target_code} изменена мышью; пересчитываю пик…")
        self.apply_manual_integration()

    def handle_batch_tree_selection(self, event=None):
        if self._batch_tree_syncing or self.batch_tree is None:
            return
        selection = self.batch_tree.selection()
        if not selection:
            return
        target_index = int(selection[0])
        if target_index != self.current_batch_index:
            self.load_batch_at_index(target_index)

    def prev_batch(self):
        if self.current_batch_index <= 0:
            return
        self.load_batch_at_index(self.current_batch_index - 1)

    def next_batch(self):
        if self.current_batch_index >= len(self.loaded_batches) - 1:
            return
        self.load_batch_at_index(self.current_batch_index + 1)

    def process_batch(self, batch: dict):
        if batch.get("processed_df") is not None:
            return batch
        batch.update(process_chromatogram_batch(batch["dataframe"], self.reference_targets))
        return batch

    def load_batch_at_index(self, index: int):
        if index < 0 or index >= len(self.loaded_batches):
            return

        batch = self.loaded_batches[index]
        self.process_batch(batch)
        self.current_batch_index = index
        self.current_sample_name = batch["sample_name"]
        self.file_var.set(
            f"Файл: {self.current_file.name} | Sample: {self.current_sample_name} | Source: {batch.get('file_name', '')}"
        )

        self.df_processed = batch["processed_df"]
        self.best_window = batch["best_window"]
        self.peaks_df = batch["peaks_df"]
        self.matched_targets_df = batch["matched_targets_df"]
        self.current_rt_shift = batch["rt_shift"]
        self.update_batch_navigation()
        self.refresh_peaks()

    def build_batch_results_rows(self, process_all: bool = False):
        rows = []
        for index, batch in enumerate(self.loaded_batches):
            if process_all:
                self.process_batch(batch)
            omega_report = batch.get("omega_report")
            if omega_report is None and isinstance(batch.get("omega"), dict):
                omega_report = batch.get("omega", {}).get("omega3_trio", np.nan)
            value = omega_report if omega_report is not None else np.nan
            value_text = f"{value:.4f}" if np.isfinite(value) else ""
            confidence = batch.get("confidence") if isinstance(batch.get("confidence"), dict) else {}
            confidence_score = confidence.get("score", np.nan)
            confidence_label = confidence.get("label", "")
            confidence_text = (
                f"{int(round(confidence_score))}/100 {confidence_label}".strip()
                if np.isfinite(confidence_score)
                else ""
            )
            rows.append((index, batch.get("sample_name", f"Batch {index + 1}"), value_text, confidence_text))
        return rows

    def _populate_batch_tree_widget(self, tree: ttk.Treeview, process_all: bool = False):
        if tree is None:
            return
        selected_iid = str(self.current_batch_index) if self.loaded_batches else None
        for item_id in tree.get_children():
            tree.delete(item_id)
        show_confidence = "confidence" in set(tree["columns"])
        for index, sample_name, value_text, confidence_text in self.build_batch_results_rows(process_all=process_all):
            values = (sample_name, value_text, confidence_text) if show_confidence else (sample_name, value_text)
            tree.insert("", "end", iid=str(index), values=values)
        if selected_iid is not None and tree.exists(selected_iid):
            self._batch_tree_syncing = True
            tree.selection_set(selected_iid)
            tree.focus(selected_iid)
            tree.see(selected_iid)
            self._batch_tree_syncing = False

    def populate_main_batch_tree(self):
        self._populate_batch_tree_widget(self.batch_tree, process_all=False)

    def copy_batch_results(self, selected_only: bool):
        if self.batch_results_tree is None:
            return

        item_ids = list(self.batch_results_tree.selection()) if selected_only else list(self.batch_results_tree.get_children())
        if not item_ids:
            return

        lines = []
        for item_id in item_ids:
            values = self.batch_results_tree.item(item_id, "values")
            lines.append("\t".join(str(v) for v in values))
        self.root.clipboard_clear()
        self.root.clipboard_append("\n".join(lines))
        self.status_var.set(f"Скопировано строк: {len(lines)}")

    def handle_batch_results_copy(self, event):
        if self.batch_results_tree is None:
            return
        self.copy_batch_results(selected_only=True)

    def jump_to_batch_from_results(self, event=None):
        if self.batch_results_tree is None:
            return
        selection = self.batch_results_tree.selection()
        if not selection:
            return
        target_index = int(selection[0])
        self.load_batch_at_index(target_index)

    def populate_batch_results_tree(self):
        self._populate_batch_tree_widget(self.batch_results_tree, process_all=True)

    def open_batch_results_window(self):
        if not self.loaded_batches:
            return

        if self.batch_results_window is not None and self.batch_results_window.winfo_exists():
            self.populate_batch_results_tree()
            self.batch_results_window.deiconify()
            self.batch_results_window.lift()
            self.batch_results_window.focus_force()
            return

        self.batch_results_window = tk.Toplevel(self.root)
        self.batch_results_window.title("Batch Results")
        self.batch_results_window.geometry("560x520")

        frame = ttk.Frame(self.batch_results_window, padding=10)
        frame.pack(fill="both", expand=True)

        toolbar = ttk.Frame(frame)
        toolbar.pack(fill="x", pady=(0, 8))
        ttk.Button(toolbar, text="Copy Selected", command=lambda: self.copy_batch_results(selected_only=True)).pack(side="left")
        ttk.Button(toolbar, text="Copy All", command=lambda: self.copy_batch_results(selected_only=False)).pack(side="left", padx=(8, 0))
        ttk.Button(toolbar, text="Open Selected", command=self.jump_to_batch_from_results).pack(side="left", padx=(8, 0))

        tree_frame = ttk.Frame(frame)
        tree_frame.pack(fill="both", expand=True)
        columns = ("sample_name", "omega_value", "confidence")
        self.batch_results_tree = ttk.Treeview(tree_frame, columns=columns, show="headings", height=18, selectmode="extended")
        self.batch_results_tree.heading("sample_name", text="Номер образца")
        self.batch_results_tree.heading("omega_value", text="Значение")
        self.batch_results_tree.heading("confidence", text="Уверенность")
        self.batch_results_tree.column("sample_name", width=300, anchor="w")
        self.batch_results_tree.column("omega_value", width=120, anchor="center")
        self.batch_results_tree.column("confidence", width=140, anchor="center")
        self.batch_results_tree.pack(side="left", fill="both", expand=True)

        scrollbar = ttk.Scrollbar(tree_frame, orient="vertical", command=self.batch_results_tree.yview)
        scrollbar.pack(side="right", fill="y")
        self.batch_results_tree.configure(yscrollcommand=scrollbar.set)
        self.batch_results_tree.bind("<Double-1>", self.jump_to_batch_from_results)
        self.batch_results_tree.bind("<Control-c>", self.handle_batch_results_copy)

        help_label = ttk.Label(frame, text="Ctrl+C копирует выделенные строки. Двойной клик открывает выбранный образец.")
        help_label.pack(fill="x", pady=(8, 0))

        def on_close():
            if self.batch_results_window is not None:
                self.batch_results_window.destroy()
            self.batch_results_window = None
            self.batch_results_tree = None

        self.batch_results_window.protocol("WM_DELETE_WINDOW", on_close)
        self.populate_batch_results_tree()

    def open_file(self):
        file_path = filedialog.askopenfilename(title="Выберите CSV", filetypes=[("CSV", "*.csv *.CSV"), ("All", "*.*")])
        if not file_path:
            return
        try:
            self.current_file = Path(file_path)
            self.loaded_batches = omega_chromatopy_clean.load_batches(self.current_file, cutoff_minutes=4.0)
            self.load_batch_at_index(0)
            self.preload_loaded_batches()
        except Exception as e:
            messagebox.showerror("Ошибка", str(e), parent=self.root)

    def preload_loaded_batches(self):
        if not self.loaded_batches:
            return
        total = len(self.loaded_batches)
        for index, batch in enumerate(self.loaded_batches):
            if batch.get("processed_df") is None:
                self.status_var.set(f"Предрасчёт batch: {index + 1}/{total}")
                self.root.update_idletasks()
                self.process_batch(batch)
        self.populate_main_batch_tree()
        if self.batch_results_window is not None and self.batch_results_window.winfo_exists():
            self.populate_batch_results_tree()

    def refresh_peaks(self):
        if self.df_processed is None:
            return
        current_batch = self.loaded_batches[self.current_batch_index] if self.loaded_batches else None
        engine = current_batch.get("engine", "") if current_batch is not None else ""
        if engine == "chromatopy_clean":
            omega = compute_clean_omega_metrics(self.matched_targets_df)
        elif engine == "omega_core":
            omega = core_metrics.compute_omega(self.matched_targets_df)
        else:
            omega = compute_omega_metrics(self.matched_targets_df)
        baseline_mode = current_batch.get("baseline_mode", "chebyshev") if current_batch is not None else "chebyshev"
        cluster_quality_score = current_batch.get("cluster_quality_score", np.nan) if current_batch is not None else np.nan
        if engine == "chromatopy_clean":
            cluster_quality_score = _compute_cluster_quality_score(self.matched_targets_df)
        confidence = build_confidence_assessment(
            matched_targets_df=self.matched_targets_df,
            peaks_df=self.peaks_df,
            omega=omega,
            baseline_mode=baseline_mode,
            cluster_quality_score=cluster_quality_score,
        )
        report_value = omega["omega3_trio"]
        if current_batch is not None:
            current_batch["processed_df"] = self.df_processed
            current_batch["best_window"] = self.best_window
            current_batch["peaks_df"] = self.peaks_df
            current_batch["matched_targets_df"] = self.matched_targets_df
            current_batch["rt_shift"] = self.current_rt_shift
            current_batch["omega"] = omega
            current_batch["omega_report"] = omega["omega3_trio"]
            current_batch["cluster_quality_score"] = cluster_quality_score
            current_batch["confidence"] = confidence
            report_value = current_batch["omega_report"]

        if np.isfinite(report_value):
            self.omega_var.set(
                f"Omega-3: {report_value:.2f}% | strict: {omega['omega3_trio_strict']:.2f}%"
            )
        else:
            self.omega_var.set("Omega-3: —")

        gamma_text = "γ-Linolenic: —"
        if not self.matched_targets_df.empty:
            gamma_match = self.matched_targets_df[self.matched_targets_df["code"] == "C18:3N6"]
            if not gamma_match.empty:
                gamma_area = float(pd.to_numeric(gamma_match["area"], errors="coerce").iloc[0])
                detected_total_area = float(pd.to_numeric(self.peaks_df["area"], errors="coerce").fillna(0.0).sum())
                gamma_percent = (100.0 * gamma_area / detected_total_area) if detected_total_area > 0 else np.nan
                if np.isfinite(gamma_area):
                    gamma_text = f"γ-Linolenic: area {gamma_area:.2f}"
                    if np.isfinite(gamma_percent):
                        gamma_text = f"{gamma_text} | peaks {gamma_percent:.2f}%"
        self.gamma_var.set(gamma_text)
        self.current_confidence = confidence
        self.confidence_var.set(confidence.get("button_text", "Уверенность: —"))
        self.confidence_button.state(["!disabled"])

        self.update_plot()
        self.update_table()

        self.status_var.set(
            f"RT shift: {self.current_rt_shift:+.3f} min | matched {int(self.matched_targets_df['matched_peak_id'].notna().sum())}/{len(self.reference_targets)}"
        )
        self.integration_var.set(
            f"Integration: {len(self.peaks_df)} peaks | SG {self.best_window}"
        )
        self.populate_main_batch_tree()
        if self.batch_results_window is not None and self.batch_results_window.winfo_exists():
            self.populate_batch_results_tree()

    def _resolve_selected_plot_items(self):
        selected_row = None
        selected_peak = None
        if self.selected_target_code and not self.matched_targets_df.empty:
            selected_match = self.matched_targets_df[self.matched_targets_df["code"] == self.selected_target_code]
            if not selected_match.empty:
                selected_row = selected_match.iloc[0]
                matched_peak_id = pd.to_numeric(selected_match["matched_peak_id"], errors="coerce").iloc[0]
                if np.isfinite(matched_peak_id) and not self.peaks_df.empty:
                    selected_peak_match = self.peaks_df[self.peaks_df["peak_id"] == int(matched_peak_id)]
                    if not selected_peak_match.empty:
                        selected_peak = selected_peak_match.iloc[0]
        return selected_row, selected_peak

    def _draw_chromatogram_axis(
        self,
        axis,
        x: np.ndarray,
        y: np.ndarray,
        y_smooth: np.ndarray,
        fill_y: np.ndarray,
        marker_y: np.ndarray,
        selected_row,
        selected_peak,
        title: str,
        x_min=None,
        x_max=None,
        compact: bool = False,
        normalized: bool = False,
    ):
        axis.clear()
        axis.set_facecolor("#fcfcfc")
        y_draw = y
        y_smooth_draw = y_smooth
        fill_y_draw = fill_y
        marker_y_draw = marker_y
        if normalized:
            visible_mask = np.ones(len(x), dtype=bool)
            if x_min is not None:
                visible_mask &= x >= float(x_min)
            if x_max is not None:
                visible_mask &= x <= float(x_max)
            if np.any(visible_mask):
                local_candidates = [
                    np.abs(y_smooth[visible_mask]),
                    np.abs(marker_y[visible_mask]),
                    np.abs(fill_y[visible_mask]),
                ]
                local_scale = max(
                    1e-9,
                    max(float(np.nanmax(values)) for values in local_candidates if values.size > 0),
                )
                y_draw = y / local_scale
                y_smooth_draw = y_smooth / local_scale
                fill_y_draw = fill_y / local_scale
                marker_y_draw = marker_y / local_scale

        axis.axhline(0.0, color="#777777", linewidth=0.8, alpha=0.55)
        axis.grid(color="#d9d9d9", linewidth=0.45, alpha=0.55)
        axis.plot(x, y_draw, linewidth=0.95 if compact else 1.0, color="#2a5b84", alpha=0.50, label="Corrected")
        axis.plot(x, y_smooth_draw, linewidth=1.05 if compact else 1.2, color="#111111", alpha=0.92, label="Smoothed")

        if not self.peaks_df.empty:
            for _, peak in self.peaks_df.iterrows():
                peak_start_x = float(peak["start_x"])
                peak_end_x = float(peak["end_x"])
                peak_apex_x = float(peak["apex_x"])
                if x_min is not None and peak_end_x < x_min:
                    continue
                if x_max is not None and peak_start_x > x_max:
                    continue
                start_idx = int(peak["start_idx"])
                end_idx = int(peak["end_idx"])
                apex_idx = int(peak["apex_idx"])
                axis.axvline(peak_start_x, color="#caa25a", linewidth=0.45 if compact else 0.5, alpha=0.18)
                axis.axvline(peak_end_x, color="#caa25a", linewidth=0.45 if compact else 0.5, alpha=0.18)
                axis.scatter(x[apex_idx], marker_y_draw[apex_idx], s=10 if compact else 16, color="#b84a35", alpha=0.75, zorder=4)

        if not self.matched_targets_df.empty:
            for _, target_row in self.matched_targets_df.iterrows():
                start_x = pd.to_numeric(pd.Series([target_row.get("integration_start_x")]), errors="coerce").iloc[0]
                end_x = pd.to_numeric(pd.Series([target_row.get("integration_end_x")]), errors="coerce").iloc[0]
                if not (np.isfinite(start_x) and np.isfinite(end_x)):
                    continue
                if x_min is not None and end_x < x_min:
                    continue
                if x_max is not None and start_x > x_max:
                    continue
                start_idx = int(np.argmin(np.abs(x - float(start_x))))
                end_idx = int(np.argmin(np.abs(x - float(end_x))))
                if end_idx <= start_idx:
                    continue
                axis.fill_between(
                    x[start_idx:end_idx + 1],
                    0.0,
                    fill_y_draw[start_idx:end_idx + 1],
                    color="#f2b134",
                    alpha=0.22 if compact else 0.24,
                    linewidth=0.0,
                    zorder=2,
                )
                axis.axvline(float(start_x), color="#d6a033", linewidth=0.65 if compact else 0.75, alpha=0.32, zorder=3)
                axis.axvline(float(end_x), color="#d6a033", linewidth=0.65 if compact else 0.75, alpha=0.32, zorder=3)

        if selected_row is not None:
            start_x = pd.to_numeric(pd.Series([selected_row.get("integration_start_x") if selected_row is not None else np.nan]), errors="coerce").iloc[0]
            end_x = pd.to_numeric(pd.Series([selected_row.get("integration_end_x") if selected_row is not None else np.nan]), errors="coerce").iloc[0]
            manual_start_x, manual_end_x = self._manual_bounds_from_vars()
            if np.isfinite(manual_start_x) and np.isfinite(manual_end_x):
                start_x, end_x = manual_start_x, manual_end_x
            apex_x = pd.to_numeric(pd.Series([selected_row.get("found_rt") if selected_row is not None else np.nan]), errors="coerce").iloc[0]
            if (not np.isfinite(start_x) or not np.isfinite(end_x)) and selected_peak is not None:
                start_x = float(selected_peak["start_x"])
                end_x = float(selected_peak["end_x"])
            if not np.isfinite(apex_x) and selected_peak is not None:
                apex_x = float(selected_peak["apex_x"])
            if not np.isfinite(apex_x) and np.isfinite(start_x) and np.isfinite(end_x):
                apex_x = 0.5 * (float(start_x) + float(end_x))
            if np.isfinite(start_x) and np.isfinite(end_x) and np.isfinite(apex_x) and (x_min is None or end_x >= x_min) and (x_max is None or start_x <= x_max):
                start_idx = int(np.argmin(np.abs(x - float(start_x))))
                end_idx = int(np.argmin(np.abs(x - float(end_x))))
                apex_idx = int(np.argmin(np.abs(x - float(apex_x))))
                axis.fill_between(
                    x[start_idx:end_idx + 1],
                    0.0,
                    fill_y_draw[start_idx:end_idx + 1],
                    color="#ff4d6d",
                    alpha=0.36 if compact else 0.40,
                    linewidth=0.0,
                    zorder=3,
                )
                axis.axvline(float(start_x), color="#ff4d6d", linewidth=1.0 if compact else 1.2, alpha=0.85, zorder=5)
                axis.axvline(float(end_x), color="#ff4d6d", linewidth=1.0 if compact else 1.2, alpha=0.85, zorder=5)
                axis.scatter(
                    x[apex_idx],
                    marker_y_draw[apex_idx],
                    s=44 if compact else 70,
                    facecolor="#fff3f6",
                    edgecolor="#ff4d6d",
                    linewidth=1.3 if compact else 1.6,
                    zorder=6,
                )

        visible_codes = []
        if not self.matched_targets_df.empty:
            labeled = self.matched_targets_df[self.matched_targets_df["matched_peak_id"].notna()].copy()
            for _, row in labeled.iterrows():
                found_rt = float(row["found_rt"])
                if x_min is not None and found_rt < x_min:
                    continue
                if x_max is not None and found_rt > x_max:
                    continue
                apex_idx = int(np.argmin(np.abs(x - found_rt)))
                visible_codes.append(str(row["code"]))
                axis.text(
                    found_rt,
                    marker_y_draw[apex_idx] + max(np.nanmax(marker_y_draw) * 0.010, 0.06 if normalized else (50.0 if compact else 80.0)),
                    str(row["code"]),
                    fontsize=6.4 if compact else 7,
                    rotation=90,
                    ha="center",
                    va="bottom",
                    color="#2c3e50",
                    alpha=0.92,
                )

        if selected_row is not None and pd.notna(selected_row.get("found_rt")):
            selected_rt = float(selected_row["found_rt"])
            if (x_min is None or selected_rt >= x_min) and (x_max is None or selected_rt <= x_max):
                apex_idx = int(np.argmin(np.abs(x - selected_rt)))
                label_text = f"{selected_row['code']}  RT {selected_rt:.4f}"
                axis.annotate(
                    label_text,
                    xy=(selected_rt, marker_y_draw[apex_idx]),
                    xytext=(10, 12 if compact else 14),
                    textcoords="offset points",
                    fontsize=7 if compact else 8,
                    color="#7a1028",
                    bbox={"boxstyle": "round,pad=0.25", "facecolor": "#fff3f6", "edgecolor": "#ff4d6d", "alpha": 0.95},
                    arrowprops={"arrowstyle": "->", "color": "#ff4d6d", "lw": 0.9},
                    zorder=7,
                )

        if x_min is not None and x_max is not None:
            axis.set_xlim(float(x_min), float(x_max))
            if visible_codes:
                peak_text = ", ".join(visible_codes)
                axis.text(
                    0.01,
                    0.98,
                    peak_text,
                    transform=axis.transAxes,
                    ha="left",
                    va="top",
                    fontsize=7,
                    color="#304860",
                    bbox={"boxstyle": "round,pad=0.20", "facecolor": "#ffffff", "edgecolor": "#d5dde5", "alpha": 0.85},
                )
            if normalized:
                visible_mask = (x >= float(x_min)) & (x <= float(x_max))
                local_min = float(np.nanmin(y_draw[visible_mask])) if np.any(visible_mask) else 0.0
                local_max = float(np.nanmax(np.maximum.reduce([
                    np.asarray(y_draw[visible_mask]),
                    np.asarray(y_smooth_draw[visible_mask]),
                    np.asarray(fill_y_draw[visible_mask]),
                ]))) if np.any(visible_mask) else 1.0
                axis.set_ylim(min(-0.08, local_min * 1.08), max(1.05, local_max * 1.12))

        axis.set_title(title, fontsize=9 if compact else 11, pad=6)
        axis.tick_params(labelsize=7 if compact else 9)
        axis.set_xlabel("Time, min", fontsize=8 if compact else 10)
        axis.set_ylabel("Norm." if normalized else "Signal", fontsize=8 if compact else 10)
        for spine in axis.spines.values():
            spine.set_color("#b8c2cc")
            spine.set_linewidth(0.8)

    def update_plot(self):
        x_col = _get_x_column_name(self.df_processed)
        x = self.df_processed[x_col].to_numpy(dtype=float)
        y = self.df_processed["y_corrected"].to_numpy(dtype=float)
        y_smooth = self.df_processed["y_smooth"].to_numpy(dtype=float)
        fill_y = np.clip(y, 0.0, None)
        marker_y = y
        selected_row, selected_peak = self._resolve_selected_plot_items()

        self._draw_chromatogram_axis(
            self.ax,
            x=x,
            y=y,
            y_smooth=y_smooth,
            fill_y=fill_y,
            marker_y=marker_y,
            selected_row=selected_row,
            selected_peak=selected_peak,
            title="Общая хроматограмма",
            compact=False,
        )
        self.ax.legend(loc="upper right", fontsize=8)

        for preview_ax, (label, x_min, x_max) in zip(self.preview_axes, self.preview_specs):
            self._draw_chromatogram_axis(
                preview_ax,
                x=x,
                y=y,
                y_smooth=y_smooth,
                fill_y=fill_y,
                marker_y=marker_y,
                selected_row=selected_row,
                selected_peak=selected_peak,
                title=f"Участок {label}",
                x_min=x_min,
                x_max=x_max,
                compact=True,
                normalized=True,
            )
        self.canvas.draw()

    def update_table(self):
        for i in self.tree.get_children():
            self.tree.delete(i)
        if self.matched_targets_df.empty:
            self.selected_target_code = None
            return
        available_codes = set()
        for _, row in self.matched_targets_df.iterrows():
            code = str(row.get("code", ""))
            available_codes.add(code)
            vals = (
                row.get("display_name", ""), row.get("code", ""),
                "" if pd.isna(row.get("expected_rt")) else f"{row['expected_rt']:.4f}",
                "" if pd.isna(row.get("found_rt")) else f"{row['found_rt']:.4f}",
                "" if pd.isna(row.get("area")) else f"{row['area']:.6f}",
                "" if pd.isna(row.get("percent_area")) else f"{row['percent_area']:.2f}",
                row.get("status", ""),
            )
            self.tree.insert("", "end", iid=code, values=vals)
        if self.selected_target_code not in available_codes:
            self.selected_target_code = None
        if self.selected_target_code is not None:
            self.tree.selection_set(self.selected_target_code)
            self.tree.focus(self.selected_target_code)


if __name__ == "__main__":
    root = tk.Tk()
    app = ChromatogramApp(root)
    root.mainloop()
