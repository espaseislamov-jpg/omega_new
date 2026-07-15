"""Guarded component-fit recovery used by the production cluster pipeline.

This module contains the small subset of the historical GUI monolith that
is still numerically active.  It deliberately has no GUI dependency.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from scipy.stats import median_abs_deviation

try:
    from lmfit.models import LinearModel, PseudoVoigtModel
except Exception:
    LinearModel = None
    PseudoVoigtModel = None

PEAK_SUPPORT_THRESHOLD_SIGMA = 0.80
PEAK_SUPPORT_THRESHOLD_FRACTION = 0.020
PEAK_SUPPORT_CONSECUTIVE_POINTS = 4
ENABLE_OVERWIDE_C22_PVFIT_REFINEMENT = True
C22_PVFIT_OVERWIDE_MEAN_WIDTH_MIN = 0.032
C22_PVFIT_OVERWIDE_DHA_WIDTH_MIN = 0.039
C22_PVFIT_OVERWIDE_C22_4_WIDTH_MIN = 0.033
C22_PVFIT_AREA_RATIO_MIN = 0.92
C22_PVFIT_AREA_RATIO_MAX = 0.995
C18_OVERLAP_START_TOLERANCE = 0.001
C20_FIT_EPA_AREA_MAX = 450.0
C20_FIT_EPA_PROMINENCE_MAX = 1500.0
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


def _recompute_matched_percent_area(matched_targets_df: pd.DataFrame) -> pd.DataFrame:
    out = matched_targets_df.copy()
    total_area = float(pd.to_numeric(out["area"], errors="coerce").fillna(0.0).sum())
    if total_area > 0:
        out["percent_area"] = 100.0 * pd.to_numeric(out["area"], errors="coerce") / total_area
    else:
        out["percent_area"] = np.nan
    return out


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

    nominal_centers = {
        "C22:6": 9.252,
        "C22:5": 9.285,
        "C22:4": 9.316,
    }
    observed_shifts = []
    for code in c22_codes:
        found_rt = pd.to_numeric(
            cluster.loc[cluster["code"] == code, "found_rt"], errors="coerce"
        ).iloc[0]
        if np.isfinite(found_rt):
            observed_shifts.append(float(found_rt) - nominal_centers[code])
    cluster_shift = float(np.median(observed_shifts)) if observed_shifts else 0.0
    shifted_centers = [nominal_centers[code] + cluster_shift for code in c22_codes]

    fit_out, _ = _refine_cluster_with_deconvolution(
        df=df,
        peaks_df=peaks_df,
        matched_targets_df=out,
        cluster_codes=c22_codes,
        default_centers=shifted_centers,
        window_left=9.22,
        window_right=9.33,
        center_tolerances=[0.012, 0.012, 0.012],
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
