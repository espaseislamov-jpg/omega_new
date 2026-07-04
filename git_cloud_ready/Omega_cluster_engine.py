from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.signal import find_peaks, savgol_filter
from scipy.stats import median_abs_deviation

import omega_core
from omega_core import clusters as core_clusters
from omega_core import matching as core_matching
from omega_core import metrics as core_metrics
from omega_core import signal as core_signal

try:
    from lmfit.models import LinearModel, PseudoVoigtModel
except Exception:  # pragma: no cover - optional dependency
    LinearModel = None
    PseudoVoigtModel = None


PROJECT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = Path(r"C:\Users\marat\Desktop\CSV_Omega")
DEFAULT_REFERENCE_PATH = PROJECT_DIR / "reference_targets_reverted_c22fixed.json"
HYBRID_CONFIDENCE_TRIGGER = 60
HYBRID_MAX_CLUSTER_OMEGA_SHIFT = 0.08
DETECTION_SAVGOL_WINDOW = 21
ANCHOR_CLUSTER_MAX_OMEGA_SHIFT = 0.12
ANCHOR_CLUSTER_MIN_OMEGA_SHIFT = 0.025
ANCHOR_MIN_POINTS = 4
ANCHOR_MIN_RMS_SHIFT_FOR_REPLACEMENT = 0.005
ANCHOR_COHERENT_SPREAD_MAX = 0.0025
ANCHOR_NOISY_SPREAD_MAX = 0.0060

ANCHOR_REFERENCE_RTS: dict[str, float] = {
    "C16:1N7": 6.600833,
    "C16:0": 6.716000,
    "C24:1N9": 10.301000,
    "C24:0": 10.394000,
}

LEARNED_CLUSTER_CENTERS: dict[str, float] = {
    "C18:3N6": 7.504667,
    "C18:2N6C": 7.592000,
    "C18:1N9C": 7.624167,
    "C18:3N3": 7.653000,
    "C18:0": 7.751333,
    "C20:4N6": 8.383000,
    "C20:5": 8.413333,
    "C20:3N8": 8.469667,
    "C22:6": 9.249167,
    "C22:5": 9.282333,
    "C22:4": 9.311667,
}


@dataclass(frozen=True)
class PeakSpec:
    code: str
    center: float
    center_tol: float = 0.055
    sigma: float = 0.012
    sigma_min: float = 0.0025
    sigma_max: float = 0.055


@dataclass(frozen=True)
class ClusterSpec:
    name: str
    left: float
    right: float
    peaks: tuple[PeakSpec, ...]


@dataclass(frozen=True)
class AnchorModel:
    count: int
    rms_shift: float
    max_abs_shift: float
    median_shift: float
    spread: float
    left_shift: float
    right_shift: float
    gradient_shift: float
    mode: str
    anchor_rows: tuple[dict, ...]

    def shift_at(self, rt: float) -> float:
        if not self.anchor_rows:
            return 0.0
        if self.mode == "coherent":
            return float(self.median_shift)
        anchor_rt = np.asarray([row["reference_rt"] for row in self.anchor_rows], dtype=float)
        shifts = np.asarray([row["shift"] for row in self.anchor_rows], dtype=float)
        if len(anchor_rt) == 1:
            return float(shifts[0])
        order = np.argsort(anchor_rt)
        return float(np.interp(float(rt), anchor_rt[order], shifts[order]))


CLUSTERS: tuple[ClusterSpec, ...] = (
    ClusterSpec(
        "C16",
        6.43,
        6.82,
        (
            PeakSpec("C16:1N7", 6.56, center_tol=0.070, sigma=0.015, sigma_max=0.065),
            PeakSpec("C16:0", 6.67, center_tol=0.070, sigma=0.018, sigma_max=0.075),
        ),
    ),
    ClusterSpec(
        "C18",
        7.42,
        7.84,
        (
            PeakSpec("C18:3N6", 7.50, center_tol=0.060, sigma=0.010, sigma_max=0.045),
            PeakSpec("C18:2N6C", 7.59, center_tol=0.065, sigma=0.012, sigma_max=0.055),
            PeakSpec("C18:1N9C", 7.63, center_tol=0.070, sigma=0.014, sigma_max=0.060),
            PeakSpec("C18:3N3", 7.655, center_tol=0.065, sigma=0.010, sigma_max=0.045),
            PeakSpec("C18:0", 7.755, center_tol=0.075, sigma=0.016, sigma_max=0.065),
        ),
    ),
    ClusterSpec(
        "C20",
        8.31,
        8.52,
        (
            PeakSpec("C20:4N6", 8.385, center_tol=0.055, sigma=0.012, sigma_max=0.050),
            PeakSpec("C20:5", 8.425, center_tol=0.055, sigma=0.010, sigma_max=0.045),
            PeakSpec("C20:3N8", 8.468, center_tol=0.055, sigma=0.010, sigma_max=0.045),
        ),
    ),
    ClusterSpec(
        "C22",
        9.18,
        9.36,
        (
            PeakSpec("C22:6", 9.253, center_tol=0.050, sigma=0.011, sigma_max=0.045),
            PeakSpec("C22:5", 9.286, center_tol=0.050, sigma=0.010, sigma_max=0.042),
            PeakSpec("C22:4", 9.316, center_tol=0.050, sigma=0.010, sigma_max=0.042),
        ),
    ),
    ClusterSpec(
        "C24",
        10.24,
        10.46,
        (
            PeakSpec("C24:1N9", 10.305, center_tol=0.065, sigma=0.012, sigma_max=0.055),
            PeakSpec("C24:0", 10.395, center_tol=0.065, sigma=0.014, sigma_max=0.060),
        ),
    ),
)


def robust_sigma(values: np.ndarray) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return 0.0
    sigma = float(median_abs_deviation(arr, scale="normal", nan_policy="omit"))
    if sigma <= 0:
        sigma = float(np.std(arr))
    return sigma


def _x_column(df: pd.DataFrame) -> str:
    if "x_corrected" in df.columns:
        return "x_corrected"
    return "x"


def _build_detection_profile(
    y: np.ndarray,
    window: int = DETECTION_SAVGOL_WINDOW,
    polyorder: int = 3,
    smooth: bool = True,
) -> np.ndarray:
    if not smooth:
        return y.copy()
    if y.size <= polyorder + 2:
        return y.copy()
    window = min(int(window), y.size if y.size % 2 else y.size - 1)
    window = max(polyorder + 2 + ((polyorder + 2) % 2 == 0), window)
    if window % 2 == 0:
        window -= 1
    return savgol_filter(y, window_length=window, polyorder=polyorder, mode="interp")


def prepare_raw_signal(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    y = out["y"].to_numpy(dtype=float)
    if y.size < 8:
        raise ValueError("Not enough points for raw signal processing.")
    out["baseline"] = np.zeros_like(y)
    out["y_corrected"] = y.copy()
    out["y_smooth"] = y.copy()
    return out


def _window_view(df: pd.DataFrame, left: float, right: float):
    x = df[_x_column(df)].to_numpy(dtype=float)
    y_raw = df["y_corrected"].to_numpy(dtype=float)
    mask = (x >= float(left)) & (x <= float(right))
    if mask.sum() < 8:
        return None
    return x[mask], y_raw[mask]


def _estimate_local_background(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    n = len(x)
    edge = max(3, min(12, n // 6))
    left_y = float(np.median(y[:edge]))
    right_y = float(np.median(y[-edge:]))
    return np.interp(x, [x[0], x[-1]], [left_y, right_y])


def _initial_centers_from_signal(
    x: np.ndarray,
    y: np.ndarray,
    specs: Iterable[PeakSpec],
) -> dict[str, float]:
    specs = list(specs)
    noise = max(robust_sigma(y), 1.0)
    prominence = max(noise * 1.2, float(np.nanmax(y)) * 0.015)
    distance = max(1, int(round(0.012 / max(float(np.median(np.diff(x))), 1e-6))))
    peaks, props = find_peaks(y, prominence=prominence, distance=distance)
    peak_x = x[peaks] if len(peaks) else np.asarray([], dtype=float)
    peak_prom = props.get("prominences", np.ones_like(peak_x))
    centers: dict[str, float] = {}
    used: set[int] = set()
    for spec in specs:
        if peak_x.size == 0:
            centers[spec.code] = spec.center
            continue
        candidates = np.where(np.abs(peak_x - spec.center) <= spec.center_tol)[0]
        if candidates.size == 0:
            centers[spec.code] = spec.center
            continue
        candidates = [int(i) for i in candidates if int(i) not in used] or [int(i) for i in candidates]
        best = min(candidates, key=lambda i: (abs(float(peak_x[i]) - spec.center), -float(peak_prom[i])))
        used.add(int(best))
        centers[spec.code] = float(peak_x[best])
    ordered = []
    last = -np.inf
    for spec in specs:
        center = max(float(centers[spec.code]), last + 0.004)
        center = min(center, spec.center + spec.center_tol)
        ordered.append(center)
        last = center
    return {spec.code: center for spec, center in zip(specs, ordered)}


def _fit_cluster(df: pd.DataFrame, cluster: ClusterSpec, smooth_detection: bool = True) -> list[dict]:
    view = _window_view(df, cluster.left, cluster.right)
    if view is None:
        return []
    x, y_raw = view
    background_seed = _estimate_local_background(x, y_raw)
    y_signal = np.clip(y_raw - background_seed, 0.0, None)
    if not np.any(y_signal > 0):
        return []
    detection_profile = _build_detection_profile(y_signal, smooth=smooth_detection)
    centers = _initial_centers_from_signal(x, detection_profile, cluster.peaks)

    if PseudoVoigtModel is None or LinearModel is None:
        return _fallback_valley_cluster(x, y_signal, cluster, centers, detection_profile)

    model = LinearModel(prefix="bkg_")
    params = model.make_params()
    slope0, intercept0 = np.polyfit(x, background_seed, deg=1)
    signal_span = max(float(np.ptp(y_raw)), robust_sigma(y_raw), 1.0)
    x_span = max(float(np.ptp(x)), 1e-6)
    params["bkg_slope"].set(
        value=float(slope0),
        min=float(slope0 - 4.0 * signal_span / x_span),
        max=float(slope0 + 4.0 * signal_span / x_span),
    )
    params["bkg_intercept"].set(
        value=float(intercept0),
        min=float(intercept0 - 4.0 * signal_span),
        max=float(intercept0 + 4.0 * signal_span),
    )

    max_y = max(float(np.nanmax(y_signal)), 1.0)
    for pos, spec in enumerate(cluster.peaks):
        prefix = f"p{pos}_"
        component = PseudoVoigtModel(prefix=prefix)
        model += component
        params.update(component.make_params())
        center0 = float(centers.get(spec.code, spec.center))
        sigma0 = float(spec.sigma)
        height0 = max(float(np.interp(center0, x, detection_profile)), max_y * 0.02)
        amp0 = max(height0 * sigma0 * math.sqrt(2 * math.pi), 1.0)
        min_center = max(cluster.left, spec.center - spec.center_tol)
        max_center = min(cluster.right, spec.center + spec.center_tol)
        if pos > 0:
            min_center = max(min_center, centers[cluster.peaks[pos - 1].code] + 0.003)
        params[f"{prefix}amplitude"].set(value=amp0, min=0.0, max=max(amp0 * 50.0, max_y * 0.5))
        params[f"{prefix}center"].set(value=center0, min=min_center, max=max_center)
        params[f"{prefix}sigma"].set(value=sigma0, min=spec.sigma_min, max=spec.sigma_max)
        params[f"{prefix}fraction"].set(value=0.45, min=0.0, max=1.0)

    noise = max(robust_sigma(y_signal - detection_profile), 1.0)
    weights = 1.0 / np.sqrt(np.clip(y_signal, 0.0, None) + noise)
    try:
        fit = model.fit(y_raw, params, x=x, weights=weights, nan_policy="omit", max_nfev=5000)
    except Exception:
        return _fallback_valley_cluster(x, y_signal, cluster, centers, detection_profile)

    components = fit.eval_components(x=x)
    rows: list[dict] = []
    for pos, spec in enumerate(cluster.peaks):
        curve = np.clip(np.asarray(components.get(f"p{pos}_", np.zeros_like(x)), dtype=float), 0.0, None)
        area = float(np.trapezoid(curve, x))
        center = float(fit.params[f"p{pos}_center"].value)
        threshold = max(float(np.nanmax(curve)) * 0.02, noise * 0.25)
        support = np.where(curve > threshold)[0]
        if support.size:
            start_x = float(x[int(support[0])])
            end_x = float(x[int(support[-1])])
        else:
            start_x = center
            end_x = center
        rows.append({
            "code": spec.code,
            "found_rt": center,
            "area": area,
            "integration_start_x": start_x,
            "integration_end_x": end_x,
            "status": f"cluster_pv_{cluster.name}",
            "matched_peak_id": pos,
        })
    return rows


def _fallback_valley_cluster(
    x: np.ndarray,
    y: np.ndarray,
    cluster: ClusterSpec,
    centers: dict[str, float],
    detection_profile: np.ndarray,
) -> list[dict]:
    apex_indices = [int(np.argmin(np.abs(x - centers[spec.code]))) for spec in cluster.peaks]

    boundaries = [0]
    for left, right in zip(apex_indices[:-1], apex_indices[1:]):
        if right <= left:
            boundaries.append(max(boundaries[-1] + 1, left + 1))
            continue
        split = int(left + np.argmin(detection_profile[left:right + 1]))
        boundaries.append(split)
    boundaries.append(len(x) - 1)
    return _cluster_rows_from_boundaries(x, y, cluster, apex_indices, boundaries, "raw_valley")


def _cluster_rows_from_boundaries(
    x: np.ndarray,
    y: np.ndarray,
    cluster: ClusterSpec,
    apex_indices: list[int],
    boundaries: list[int],
    status: str,
) -> list[dict]:
    rows = []
    for pos, spec in enumerate(cluster.peaks):
        start = int(boundaries[pos])
        end = int(boundaries[pos + 1])
        if end <= start:
            area = 0.0
        else:
            area = float(np.trapezoid(y[start:end + 1], x[start:end + 1]))
        rows.append({
            "code": spec.code,
            "found_rt": float(x[apex_indices[pos]]),
            "area": area,
            "integration_start_x": float(x[start]),
            "integration_end_x": float(x[end]),
            "status": f"{status}_{cluster.name}",
            "matched_peak_id": pos,
        })
    return rows


def estimate_anchor_model(matched_targets: pd.DataFrame) -> AnchorModel:
    if matched_targets is None or matched_targets.empty:
        return AnchorModel(
            0,
            float("nan"),
            float("nan"),
            float("nan"),
            float("nan"),
            float("nan"),
            float("nan"),
            float("nan"),
            "missing",
            tuple(),
        )

    rows: list[dict] = []
    table = matched_targets.copy()
    table["found_rt"] = pd.to_numeric(table.get("found_rt"), errors="coerce")
    table["area"] = pd.to_numeric(table.get("area"), errors="coerce")
    table["integration_start_x"] = pd.to_numeric(table.get("integration_start_x"), errors="coerce")
    table["integration_end_x"] = pd.to_numeric(table.get("integration_end_x"), errors="coerce")
    for code, reference_rt in ANCHOR_REFERENCE_RTS.items():
        hit = table[table["code"] == code]
        if hit.empty:
            continue
        row = hit.iloc[0]
        found_rt = float(row["found_rt"])
        area = float(row["area"])
        width = float(row["integration_end_x"] - row["integration_start_x"])
        if not (np.isfinite(found_rt) and np.isfinite(area) and np.isfinite(width)):
            continue
        if area <= 0.0 or width <= 0.004 or width >= 0.120:
            continue
        shift = found_rt - float(reference_rt)
        if abs(shift) > 0.060:
            continue
        rows.append({
            "code": code,
            "reference_rt": float(reference_rt),
            "found_rt": found_rt,
            "shift": float(shift),
            "area": area,
            "width": width,
        })

    if rows:
        shifts = np.asarray([row["shift"] for row in rows], dtype=float)
        rms_shift = float(np.sqrt(np.mean(np.square(shifts))))
        max_abs_shift = float(np.max(np.abs(shifts)))
        median_shift = float(np.median(shifts))
        spread = float(np.max(shifts) - np.min(shifts))
        by_code = {row["code"]: row["shift"] for row in rows}
        left_values = [by_code[code] for code in ("C16:1N7", "C16:0") if code in by_code]
        right_values = [by_code[code] for code in ("C24:1N9", "C24:0") if code in by_code]
        left_shift = float(np.mean(left_values)) if left_values else float("nan")
        right_shift = float(np.mean(right_values)) if right_values else float("nan")
        gradient_shift = (
            float(right_shift - left_shift)
            if np.isfinite(left_shift) and np.isfinite(right_shift)
            else float("nan")
        )
    else:
        rms_shift = float("nan")
        max_abs_shift = float("nan")
        median_shift = float("nan")
        spread = float("nan")
        left_shift = float("nan")
        right_shift = float("nan")
        gradient_shift = float("nan")

    if len(rows) < ANCHOR_MIN_POINTS:
        mode = "missing"
    elif not np.isfinite(rms_shift) or rms_shift < ANCHOR_MIN_RMS_SHIFT_FOR_REPLACEMENT:
        mode = "quiet"
    elif np.isfinite(spread) and spread <= ANCHOR_COHERENT_SPREAD_MAX:
        mode = "coherent"
    elif np.isfinite(spread) and spread <= ANCHOR_NOISY_SPREAD_MAX:
        mode = "gradient"
    else:
        mode = "noisy"
    return AnchorModel(
        len(rows),
        rms_shift,
        max_abs_shift,
        median_shift,
        spread,
        left_shift,
        right_shift,
        gradient_shift,
        mode,
        tuple(rows),
    )


def _anchored_center(spec: PeakSpec, anchor_model: AnchorModel) -> float:
    reference_center = float(LEARNED_CLUSTER_CENTERS.get(spec.code, spec.center))
    return reference_center + anchor_model.shift_at(reference_center)


def _cluster_reference_center(cluster: ClusterSpec) -> float:
    values = [float(LEARNED_CLUSTER_CENTERS.get(spec.code, spec.center)) for spec in cluster.peaks]
    return float(np.mean(values))


def _anchor_allows_cluster(anchor_model: AnchorModel, cluster: ClusterSpec) -> bool:
    if anchor_model.count < ANCHOR_MIN_POINTS:
        return False
    if anchor_model.mode not in {"coherent", "gradient"}:
        return False
    local_shift = abs(anchor_model.shift_at(_cluster_reference_center(cluster)))
    return bool(local_shift >= ANCHOR_MIN_RMS_SHIFT_FOR_REPLACEMENT)


def _gaussian_curve(x: np.ndarray, height: float, center: float, sigma: float) -> np.ndarray:
    sigma = max(float(sigma), 1e-6)
    return float(height) * np.exp(-0.5 * np.square((x - float(center)) / sigma))


def _fit_cluster_gaussian_anchored(
    df: pd.DataFrame,
    cluster: ClusterSpec,
    anchor_model: AnchorModel,
    smooth_detection: bool = True,
) -> list[dict]:
    view = _window_view(df, cluster.left + anchor_model.shift_at(cluster.left), cluster.right + anchor_model.shift_at(cluster.right))
    if view is None:
        return []

    x, y_raw = view
    background_seed = _estimate_local_background(x, y_raw)
    y_signal = np.clip(y_raw - background_seed, 0.0, None)
    if not np.any(y_signal > 0):
        return []

    detection_profile = _build_detection_profile(y_signal, smooth=smooth_detection)
    n = len(cluster.peaks)
    x0 = float(np.median(x))
    dx = max(float(np.median(np.diff(x))), 1e-6)
    noise = max(robust_sigma(y_signal - detection_profile), 1.0)
    signal_span = max(float(np.ptp(y_raw)), robust_sigma(y_raw), 1.0)

    initial = []
    lower = []
    upper = []
    for pos, spec in enumerate(cluster.peaks):
        center0 = _anchored_center(spec, anchor_model)
        local_left = max(float(x[0]), center0 - min(spec.center_tol, 0.018))
        local_right = min(float(x[-1]), center0 + min(spec.center_tol, 0.018))
        local_mask = (x >= local_left) & (x <= local_right)
        if np.any(local_mask):
            local_idx = np.flatnonzero(local_mask)
            apex_idx = int(local_idx[np.argmax(detection_profile[local_idx])])
            center0 = float(x[apex_idx])
            height0 = max(float(detection_profile[apex_idx]), float(np.nanmax(y_signal)) * 0.025, noise)
        else:
            height0 = max(float(np.interp(center0, x, detection_profile)), float(np.nanmax(y_signal)) * 0.025, noise)
        sigma0 = float(np.clip(spec.sigma, spec.sigma_min, min(spec.sigma_max, 0.035)))
        initial.extend([height0, center0, sigma0])
        lower.extend([0.0, max(float(x[0]), center0 - spec.center_tol), spec.sigma_min])
        upper.extend([max(signal_span * 3.0, height0 * 20.0), min(float(x[-1]), center0 + spec.center_tol), min(spec.sigma_max, 0.055)])

    slope0, intercept0 = np.polyfit(x - x0, background_seed, deg=1)
    initial.extend([float(intercept0), float(slope0)])
    lower.extend([float(intercept0 - 3.0 * signal_span), float(slope0 - 3.0 * signal_span / max(np.ptp(x), dx))])
    upper.extend([float(intercept0 + 3.0 * signal_span), float(slope0 + 3.0 * signal_span / max(np.ptp(x), dx))])

    p0 = np.asarray(initial, dtype=float)
    lower = np.asarray(lower, dtype=float)
    upper = np.asarray(upper, dtype=float)
    p0 = np.clip(p0, lower + 1e-9, upper - 1e-9)
    weights = 1.0 / np.sqrt(np.clip(y_signal, 0.0, None) + noise)

    def unpack(params: np.ndarray):
        components = []
        for pos in range(n):
            base = 3 * pos
            components.append((float(params[base]), float(params[base + 1]), float(params[base + 2])))
        intercept = float(params[3 * n])
        slope = float(params[3 * n + 1])
        return components, intercept, slope

    def residuals(params: np.ndarray) -> np.ndarray:
        components, intercept, slope = unpack(params)
        y_fit = intercept + slope * (x - x0)
        for height, center, sigma in components:
            y_fit = y_fit + _gaussian_curve(x, height, center, sigma)
        residual = (y_fit - y_raw) * weights
        penalties = []
        centers = [center for _, center, _ in components]
        for left_center, right_center in zip(centers[:-1], centers[1:]):
            gap = right_center - left_center
            if gap < 0.004:
                penalties.append((0.004 - gap) * 1e4)
        for (height, center, sigma), spec in zip(components, cluster.peaks):
            reference_center = _anchored_center(spec, anchor_model)
            penalties.append(max(0.0, abs(center - reference_center) - spec.center_tol * 0.65) * 300.0)
            if height > 0 and sigma > 0.045:
                penalties.append((sigma - 0.045) * 200.0)
        if penalties:
            return np.concatenate([residual, np.asarray(penalties, dtype=float)])
        return residual

    try:
        fit = least_squares(
            residuals,
            p0,
            bounds=(lower, upper),
            loss="soft_l1",
            f_scale=1.5,
            max_nfev=3000,
        )
    except Exception:
        return []
    if not fit.success and fit.cost > 0:
        return []

    components, intercept, slope = unpack(fit.x)
    rows = []
    ordered = list(zip(cluster.peaks, components))
    if any(right[1][1] - left[1][1] <= 0.002 for left, right in zip(ordered[:-1], ordered[1:])):
        return []

    for pos, (spec, (height, center, sigma)) in enumerate(ordered):
        curve = np.clip(_gaussian_curve(x, height, center, sigma), 0.0, None)
        area = float(max(height, 0.0) * max(sigma, 1e-9) * math.sqrt(2.0 * math.pi))
        threshold = max(float(np.nanmax(curve)) * 0.015, noise * 0.20)
        support = np.where(curve > threshold)[0]
        if support.size:
            start_x = float(x[int(support[0])])
            end_x = float(x[int(support[-1])])
        else:
            start_x = float(center - 2.5 * sigma)
            end_x = float(center + 2.5 * sigma)
        rows.append({
            "code": spec.code,
            "found_rt": float(center),
            "area": area,
            "integration_start_x": start_x,
            "integration_end_x": end_x,
            "status": f"anchored_gauss_{cluster.name}",
            "matched_peak_id": pos,
        })
    return rows


def _build_matched_from_rows(reference_targets: pd.DataFrame, rows_by_code: dict[str, dict]) -> pd.DataFrame:
    target_rows = []
    for _, target in reference_targets.reset_index(drop=True).iterrows():
        code = str(target["code"])
        row = rows_by_code.get(code)
        if row is None:
            row = {
                "code": code,
                "found_rt": np.nan,
                "area": np.nan,
                "integration_start_x": np.nan,
                "integration_end_x": np.nan,
                "status": "not_found",
                "matched_peak_id": np.nan,
            }
        target_rows.append({
            "display_name": target.get("display_name", code),
            "code": code,
            "expected_rt": target.get("expected_rt", np.nan),
            "found_rt": row["found_rt"],
            "area": row["area"],
            "integration_start_x": row["integration_start_x"],
            "integration_end_x": row["integration_end_x"],
            "status": row["status"],
            "matched_peak_id": row["matched_peak_id"],
        })
    matched = pd.DataFrame(target_rows)
    valid_total = float(pd.to_numeric(matched["area"], errors="coerce").fillna(0.0).sum())
    matched["percent_area"] = pd.to_numeric(matched["area"], errors="coerce") / valid_total * 100.0 if valid_total > 0 else np.nan
    return matched


def _add_raw_profile_and_derivatives(processed: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    out = processed.copy()
    x = out["x_corrected"].to_numpy(dtype=float)
    y = out["y_corrected"].to_numpy(dtype=float)
    out["y_smooth"] = y.copy()
    out["dy"] = np.gradient(y, x)
    out["d2y"] = np.gradient(out["dy"].to_numpy(dtype=float), x)
    return out, 1


def process_chromatogram_core_nosmooth(dataframe: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    processed = core_signal.add_baseline(dataframe, **core_signal.BASELINE_KWARGS)
    processed, best_window = _add_raw_profile_and_derivatives(processed)
    peaks = core_signal.detect_peak_candidates(processed, best_window=best_window)
    matched_targets, rt_shift = core_matching.match_targets_to_peaks(reference_targets, peaks)
    peaks, matched_targets = core_clusters.refine_cluster_matches(processed, peaks, matched_targets)
    omega = core_metrics.compute_omega(matched_targets)
    return core_metrics.annotate_result(
        {
            "processed_df": processed,
            "best_window": best_window,
            "peaks_df": peaks,
            "matched_targets_df": matched_targets,
            "rt_shift": rt_shift,
            "omega": omega,
            "omega_report": omega["omega3_trio"],
            "signal_mode": "chebyshev_no_smooth",
        },
        baseline_mode="chebyshev_no_smooth",
    )


def process_chromatogram_core_raw_boundaries(dataframe: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    old_peak_weight = core_signal.PEAK_BOUNDARY_SIGNAL_SMOOTH_WEIGHT
    old_cluster_weight = core_clusters.CLUSTER_BOUNDARY_METRIC_SMOOTH_WEIGHT
    try:
        core_signal.PEAK_BOUNDARY_SIGNAL_SMOOTH_WEIGHT = 0.0
        core_clusters.CLUSTER_BOUNDARY_METRIC_SMOOTH_WEIGHT = 0.0
        result = omega_core.process_batch(dataframe, reference_targets)
    finally:
        core_signal.PEAK_BOUNDARY_SIGNAL_SMOOTH_WEIGHT = old_peak_weight
        core_clusters.CLUSTER_BOUNDARY_METRIC_SMOOTH_WEIGHT = old_cluster_weight

    out = dict(result)
    out["signal_mode"] = f"{result.get('baseline_mode', 'production')}_raw_boundaries"
    return out


def process_chromatogram_clustered(
    dataframe: pd.DataFrame,
    reference_targets: pd.DataFrame,
    smooth_detection: bool = True,
) -> dict:
    processed = prepare_raw_signal(dataframe)
    rows_by_code: dict[str, dict] = {}
    for cluster in CLUSTERS:
        for row in _fit_cluster(processed, cluster, smooth_detection=smooth_detection):
            rows_by_code[row["code"]] = row

    matched = _build_matched_from_rows(reference_targets, rows_by_code)
    strict = compute_strict_omega(matched)
    corrected_omega = core_metrics.compute_omega(matched)
    return {
        "processed_df": processed,
        "matched_targets_df": matched,
        "omega_strict_cluster": strict,
        "omega_report": strict,
        "omega_corrected": corrected_omega,
        "signal_mode": "raw_direct" if smooth_detection else "raw_direct_no_smooth_detection",
    }


def process_chromatogram_clustered_nosmooth(dataframe: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    return process_chromatogram_clustered(dataframe, reference_targets, smooth_detection=False)


def process_chromatogram_anchored_gaussian(
    dataframe: pd.DataFrame,
    reference_targets: pd.DataFrame,
    smooth_detection: bool = True,
) -> dict:
    base = omega_core.process_batch(dataframe, reference_targets)
    anchor_model = estimate_anchor_model(base["matched_targets_df"])
    out = dict(base)
    out["base"] = base
    suffix = "anchored_gaussian" if smooth_detection else "anchored_gaussian_no_smooth_detection"
    out["signal_mode"] = f"{base.get('baseline_mode', 'production')}_base_with_{suffix}"
    out["accepted_clusters"] = []
    out["anchor_count"] = anchor_model.count
    out["anchor_rms_shift"] = anchor_model.rms_shift
    out["anchor_max_abs_shift"] = anchor_model.max_abs_shift
    out["anchor_median_shift"] = anchor_model.median_shift
    out["anchor_spread"] = anchor_model.spread
    out["anchor_left_shift"] = anchor_model.left_shift
    out["anchor_right_shift"] = anchor_model.right_shift
    out["anchor_gradient_shift"] = anchor_model.gradient_shift
    out["anchor_mode"] = anchor_model.mode
    out["anchor_rows"] = list(anchor_model.anchor_rows)
    if anchor_model.mode not in {"coherent", "gradient"}:
        return out

    processed = prepare_raw_signal(dataframe)
    merged = base["matched_targets_df"].copy()
    base_omega = float(base["omega_report"])
    fit_rows_by_cluster: dict[str, list[dict]] = {}
    for cluster in CLUSTERS:
        if cluster.name not in {"C18", "C20", "C22"}:
            continue
        if not _anchor_allows_cluster(anchor_model, cluster):
            continue
        rows = _fit_cluster_gaussian_anchored(processed, cluster, anchor_model, smooth_detection=smooth_detection)
        if rows:
            fit_rows_by_cluster[cluster.name] = rows

    out["fit_candidate_clusters"] = ",".join(sorted(fit_rows_by_cluster))
    for cluster in CLUSTERS:
        rows = fit_rows_by_cluster.get(cluster.name)
        if not rows:
            continue
        fit_matched = _build_matched_from_rows(reference_targets, {row["code"]: row for row in rows})
        if not _cluster_fit_is_plausible(merged, fit_matched, cluster):
            continue

        tentative = merged.copy()
        for row in rows:
            target_idx = tentative.index[tentative["code"] == row["code"]]
            if len(target_idx) == 0:
                continue
            idx = int(target_idx[0])
            for column in ["found_rt", "area", "integration_start_x", "integration_end_x", "status", "matched_peak_id"]:
                tentative.at[idx, column] = row[column]
            tentative.at[idx, "status"] = f"{row['status']}_hybrid"
        tentative_total = float(pd.to_numeric(tentative["area"], errors="coerce").fillna(0.0).sum())
        tentative["percent_area"] = pd.to_numeric(tentative["area"], errors="coerce") / tentative_total * 100.0 if tentative_total > 0 else np.nan
        tentative_omega = float(core_metrics.compute_omega(tentative)["omega3_trio"])
        if np.isfinite(base_omega) and np.isfinite(tentative_omega):
            omega_shift = abs(tentative_omega - base_omega)
            if omega_shift < ANCHOR_CLUSTER_MIN_OMEGA_SHIFT:
                continue
            if omega_shift > ANCHOR_CLUSTER_MAX_OMEGA_SHIFT:
                continue
        merged = tentative
        out["accepted_clusters"].append(cluster.name)

    valid_total = float(pd.to_numeric(merged["area"], errors="coerce").fillna(0.0).sum())
    merged["percent_area"] = pd.to_numeric(merged["area"], errors="coerce") / valid_total * 100.0 if valid_total > 0 else np.nan
    omega = core_metrics.compute_omega(merged)
    out["processed_df"] = base["processed_df"]
    out["matched_targets_df"] = merged
    out["omega"] = omega
    out["omega_report"] = omega["omega3_trio"]
    out["accepted_clusters"] = out["accepted_clusters"]
    return out


def process_chromatogram_anchored_gaussian_nosmooth(
    dataframe: pd.DataFrame,
    reference_targets: pd.DataFrame,
) -> dict:
    return process_chromatogram_anchored_gaussian(dataframe, reference_targets, smooth_detection=False)


def _cluster_area(frame: pd.DataFrame, codes: Iterable[str]) -> float:
    values = pd.to_numeric(frame.loc[frame["code"].isin(list(codes)), "area"], errors="coerce")
    return float(values.fillna(0.0).sum())


def _cluster_fit_is_plausible(
    base_matched: pd.DataFrame,
    fit_matched: pd.DataFrame,
    cluster: ClusterSpec,
) -> bool:
    codes = [spec.code for spec in cluster.peaks]
    base_total = _cluster_area(base_matched, codes)
    fit_total = _cluster_area(fit_matched, codes)
    if base_total <= 0 or fit_total <= 0:
        return False
    total_ratio = fit_total / base_total
    if not 0.72 <= total_ratio <= 1.24:
        return False

    fit_rows = fit_matched[fit_matched["code"].isin(codes)].copy()
    if len(fit_rows) != len(codes):
        return False
    fit_rows["found_rt"] = pd.to_numeric(fit_rows["found_rt"], errors="coerce")
    fit_rows["integration_start_x"] = pd.to_numeric(fit_rows["integration_start_x"], errors="coerce")
    fit_rows["integration_end_x"] = pd.to_numeric(fit_rows["integration_end_x"], errors="coerce")
    if fit_rows[["found_rt", "integration_start_x", "integration_end_x"]].isna().any().any():
        return False

    by_code = fit_rows.set_index("code")
    rts = np.asarray([float(by_code.at[code, "found_rt"]) for code in codes], dtype=float)
    if np.any(np.diff(rts) <= 0.002):
        return False
    widths = np.asarray([
        float(by_code.at[code, "integration_end_x"] - by_code.at[code, "integration_start_x"])
        for code in codes
    ], dtype=float)
    if np.any(widths <= 0.004) or np.any(widths > 0.13):
        return False
    return True


def process_chromatogram_hybrid(dataframe: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    base = omega_core.process_batch(dataframe, reference_targets)
    confidence = base.get("confidence", {})
    confidence_score = float(confidence.get("score", 100.0)) if isinstance(confidence, dict) else 100.0
    if confidence_score >= HYBRID_CONFIDENCE_TRIGGER:
        out = dict(base)
        out["accepted_clusters"] = []
        out["fit_candidate"] = None
        out["base"] = base
        out["signal_mode"] = f"{base.get('baseline_mode', 'production')}_base"
        return out

    fit = process_chromatogram_clustered(dataframe, reference_targets)
    merged = base["matched_targets_df"].copy()
    fit_matched = fit["matched_targets_df"]
    accepted_clusters = []
    base_omega = float(base["omega_report"])

    for cluster in CLUSTERS:
        if cluster.name not in {"C18", "C20", "C22"}:
            continue
        if not _cluster_fit_is_plausible(merged, fit_matched, cluster):
            continue
        tentative = merged.copy()
        for spec in cluster.peaks:
            source = fit_matched[fit_matched["code"] == spec.code]
            target_idx = tentative.index[tentative["code"] == spec.code]
            if source.empty or len(target_idx) == 0:
                continue
            source_row = source.iloc[0]
            idx = target_idx[0]
            for column in ["found_rt", "area", "integration_start_x", "integration_end_x", "status", "matched_peak_id"]:
                tentative.at[idx, column] = source_row[column]
            tentative.at[idx, "status"] = f"{source_row['status']}_hybrid"
        tentative_total = float(pd.to_numeric(tentative["area"], errors="coerce").fillna(0.0).sum())
        tentative["percent_area"] = pd.to_numeric(tentative["area"], errors="coerce") / tentative_total * 100.0 if tentative_total > 0 else np.nan
        tentative_omega = float(core_metrics.compute_omega(tentative)["omega3_trio"])
        if np.isfinite(base_omega) and np.isfinite(tentative_omega):
            if abs(tentative_omega - base_omega) > HYBRID_MAX_CLUSTER_OMEGA_SHIFT:
                continue
        merged = tentative
        accepted_clusters.append(cluster.name)

    valid_total = float(pd.to_numeric(merged["area"], errors="coerce").fillna(0.0).sum())
    merged["percent_area"] = pd.to_numeric(merged["area"], errors="coerce") / valid_total * 100.0 if valid_total > 0 else np.nan
    omega = core_metrics.compute_omega(merged)
    return {
        "processed_df": base["processed_df"],
        "matched_targets_df": merged,
        "omega_report": omega["omega3_trio"],
        "omega": omega,
        "accepted_clusters": accepted_clusters,
        "fit_candidate": fit,
        "base": base,
        "signal_mode": f"{base.get('baseline_mode', 'production')}_base_with_raw_candidates",
    }


def compute_strict_omega(matched_targets_df: pd.DataFrame) -> float:
    areas = pd.to_numeric(matched_targets_df["area"], errors="coerce").fillna(0.0)
    total = float(areas.sum())
    if total <= 0:
        return float("nan")

    def area_of(code: str) -> float:
        row = matched_targets_df.loc[matched_targets_df["code"] == code, "area"]
        return float(pd.to_numeric(row, errors="coerce").fillna(0.0).iloc[0]) if not row.empty else 0.0

    return 100.0 * (area_of("C20:5") + area_of("C22:5") + area_of("C22:6")) / total


def load_excel_refs(xlsx_path: Path) -> list[tuple[int, float]]:
    df = pd.read_excel(xlsx_path, header=None)
    refs: list[tuple[int, float]] = []
    for row in df.itertuples(index=False):
        sample_no = None
        ref_value = None
        for value in row:
            if pd.isna(value):
                continue
            if sample_no is None:
                text = str(value).strip()
                if text.isdigit():
                    sample_no = int(text)
                    continue
            if ref_value is None:
                try:
                    ref_value = float(str(value).strip().replace(",", "."))
                    continue
                except ValueError:
                    pass
        if sample_no is not None and ref_value is not None:
            refs.append((sample_no, ref_value))
    return refs


def validate_batches(
    data_dir: Path = DEFAULT_DATA_DIR,
    reference_path: Path = DEFAULT_REFERENCE_PATH,
    mode: str = "hybrid",
) -> pd.DataFrame:
    reference = omega_core.load_reference_targets(reference_path)
    process_map = {
        "core-nosmooth": process_chromatogram_core_nosmooth,
        "core-raw-boundaries": process_chromatogram_core_raw_boundaries,
        "hybrid": process_chromatogram_hybrid,
        "clustered": process_chromatogram_clustered,
        "clustered-nosmooth": process_chromatogram_clustered_nosmooth,
        "anchored": process_chromatogram_anchored_gaussian,
        "anchored-nosmooth": process_chromatogram_anchored_gaussian_nosmooth,
    }
    process = process_map[mode]
    rows = []
    for xlsx_path in sorted(data_dir.glob("test_bigbatch_*.xlsx")):
        date = xlsx_path.stem.split("_")[-1]
        csv_path = data_dir / f"{date}.CSV"
        if not csv_path.exists():
            continue
        refs = load_excel_refs(xlsx_path)
        batches = omega_core.load_batches(csv_path, cutoff_minutes=4.0)
        for sample_no, ref_value in refs:
            if sample_no < 1 or sample_no > len(batches):
                continue
            batch = batches[sample_no - 1]
            result = process(batch["dataframe"], reference)
            calc = float(result["omega_report"])
            delta = calc - ref_value
            rows.append({
                "date": date,
                "sample_no": sample_no,
                "sample_name": batch.get("sample_name", ""),
                "reference": ref_value,
                "calculated": calc,
                "delta": delta,
                "abs_delta": abs(delta),
                "mode": mode,
                "signal_mode": result.get("signal_mode", "legacy_hybrid_base"),
                "accepted_clusters": ",".join(result.get("accepted_clusters", [])),
                "anchor_count": result.get("anchor_count", np.nan),
                "anchor_rms_shift": result.get("anchor_rms_shift", np.nan),
                "anchor_max_abs_shift": result.get("anchor_max_abs_shift", np.nan),
                "anchor_median_shift": result.get("anchor_median_shift", np.nan),
                "anchor_spread": result.get("anchor_spread", np.nan),
                "anchor_left_shift": result.get("anchor_left_shift", np.nan),
                "anchor_right_shift": result.get("anchor_right_shift", np.nan),
                "anchor_gradient_shift": result.get("anchor_gradient_shift", np.nan),
                "anchor_mode": result.get("anchor_mode", ""),
                "fit_candidate_clusters": result.get("fit_candidate_clusters", ""),
            })
    return pd.DataFrame(rows)


def summarize(results: pd.DataFrame) -> dict:
    if results.empty:
        return {}
    abs_delta = results["abs_delta"]
    delta = results["delta"]
    return {
        "n": int(len(results)),
        "MAE": float(abs_delta.mean()),
        "RMSE": float(math.sqrt(float(np.mean(np.square(delta))))),
        "within_0_2": int((abs_delta <= 0.2).sum()),
        "within_0_3": int((abs_delta <= 0.3).sum()),
        "within_0_4": int((abs_delta <= 0.4).sum()),
        "within_0_5": int((abs_delta <= 0.5).sum()),
        "within_0_6": int((abs_delta <= 0.6).sum()),
        "max_abs": float(abs_delta.max()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Experimental library-based Omega cluster engine.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE_PATH)
    parser.add_argument(
        "--mode",
        choices=[
            "hybrid",
            "core-nosmooth",
            "core-raw-boundaries",
            "clustered",
            "clustered-nosmooth",
            "anchored",
            "anchored-nosmooth",
        ],
        default="hybrid",
    )
    parser.add_argument("--fit-backend", choices=["valley", "pv"], default="valley")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    global PseudoVoigtModel
    if args.fit_backend == "valley":
        PseudoVoigtModel = None

    results = validate_batches(args.data_dir, args.reference, mode=args.mode)
    print("ALL", summarize(results))
    old45 = results[results["date"].isin({"13032026", "14012026", "20032026"})]
    print("OLD45", summarize(old45))
    if not results.empty:
        print("Worst")
        print(results.sort_values("abs_delta", ascending=False).head(12).to_string(index=False))
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        results.to_excel(args.out, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
