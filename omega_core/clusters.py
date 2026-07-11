from __future__ import annotations

import itertools
import os

import numpy as np
import pandas as pd

from . import legacy_fit, metrics, rt_profile
from .signal import (
    CHEMSTATION_INITIAL_AREA_REJECT,
    CHEMSTATION_INITIAL_THRESHOLD,
    _extract_peak_geometry,
    _get_x_column_name,
    _merge_peak_records,
    _robust_sigma,
)


PEAK_SUPPORT_THRESHOLD_SIGMA = 0.80
PEAK_SUPPORT_THRESHOLD_FRACTION = 0.020
PEAK_SUPPORT_CONSECUTIVE_POINTS = 4
C18_OVERLAP_START_TOLERANCE = 0.001
CLUSTER_BOUNDARY_METRIC_SMOOTH_WEIGHT = 0.60
C20_LOCAL_AREA_RATIO_TRIGGER = 1.05
C20_LOCAL_BOUNDARY_EXTENSION = 0.004
ENABLE_C22_TAIL_TIGHTENING = True
C22_TAIL_TIGHTENING_MEAN_WIDTH = 0.036
C22_TAIL_TIGHTENING_WIDTH_SCALE = 0.95
C22_TAIL_TIGHTENING_AREA_RATIO_MIN = 0.92
C22_TAIL_TIGHTENING_AREA_RATIO_MAX = 0.995
C22_TAIL_TIGHTENING_DPA_RATIO_TRIGGER = 0.90
C22_TAIL_TIGHTENING_DHA_WIDTH_TRIGGER = 0.040
ENABLE_FINAL_BOUNDARY_EXPANSION = True
EXPERIMENTAL_FORCE_CLUSTER_VALLEYS = os.environ.get("OMEGA_FORCE_CLUSTER_VALLEYS", "0").strip() == "1"
FINAL_BOUNDARY_MAX_EXTENSION = 0.028
FINAL_BOUNDARY_MAX_WIDTH = 0.120
FINAL_BOUNDARY_MAX_AREA_RATIO = 1.22
FINAL_BOUNDARY_MIN_AREA_RATIO = 1.003
FINAL_BOUNDARY_THRESHOLD_FRACTION = 0.003
FINAL_BOUNDARY_THRESHOLD_SIGMA = 0.20
FINAL_BOUNDARY_CONSECUTIVE_POINTS = 3
FINAL_BOUNDARY_MAX_CHANGED_PEAKS = 16
FINAL_BOUNDARY_MAX_OMEGA_SHIFT = 0.400
FINAL_BOUNDARY_MAX_STRICT_SPREAD_INCREASE = 0.250
FINAL_BOUNDARY_FALLBACK_MAX_CHANGED_PEAKS = 6
FINAL_BOUNDARY_FALLBACK_MAX_AREA_RATIO = 1.080
FINAL_BOUNDARY_FALLBACK_MAX_OMEGA_SHIFT = 0.100
FINAL_BOUNDARY_FALLBACK_MAX_STRICT_SPREAD_INCREASE = 0.080
ENABLE_TARGET_RT_CORRIDOR_GUARD = os.environ.get("OMEGA_TARGET_RT_CORRIDOR_GUARD", "0").strip() == "1"
TARGET_RT_CORRIDOR_GUARD_CODES = {"C18:2N6C", "C18:1N9C", "C18:3N3", "C18:0", "C20:4N6", "C20:5", "C20:3N8", "C22:6", "C22:5", "C22:4"}
TARGET_RT_CORRIDOR_MIN_WIDTH = 0.004
JUDGE_DECISIONS_ATTR = "judge_decisions"
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


def _recompute_matched_percent_area(matched_targets: pd.DataFrame) -> pd.DataFrame:
    out = matched_targets.copy()
    total_area = float(pd.to_numeric(out["area"], errors="coerce").fillna(0.0).sum())
    if total_area > 0:
        out["percent_area"] = 100.0 * pd.to_numeric(out["area"], errors="coerce") / total_area
    else:
        out["percent_area"] = np.nan
    return out


def _cluster_has_integration_overlap(matched_targets: pd.DataFrame, cluster_codes) -> bool:
    cluster = matched_targets[matched_targets["code"].isin(cluster_codes)].copy()
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


def _cluster_has_duplicate_peak_ids(matched_targets: pd.DataFrame, cluster_codes) -> bool:
    cluster = matched_targets[matched_targets["code"].isin(cluster_codes)].copy()
    if cluster.empty:
        return False
    peak_ids = pd.to_numeric(cluster["matched_peak_id"], errors="coerce").dropna().astype(int)
    return bool(peak_ids.duplicated().any())


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
    min_prominence = max(float(min_prominence), CHEMSTATION_INITIAL_THRESHOLD)
    min_area = max(float(min_area), CHEMSTATION_INITIAL_AREA_REJECT)
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
    min_apex_gaps=None,
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
        if min_apex_gaps is not None:
            required_gaps = [float(value) for value in min_apex_gaps]
            observed_gaps = np.diff(chosen["apex_x"].to_numpy(dtype=float))
            if any(float(gap) < required_gaps[min(idx, len(required_gaps) - 1)] for idx, gap in enumerate(observed_gaps)):
                continue
        score = (
            float(sum(distances))
            - 1e-6 * float(chosen["prominence"].sum())
            - 1e-7 * float(chosen["area"].sum())
        )
        if best_choice is None or score < best_choice[0]:
            best_choice = (score, chosen)

    return None if best_choice is None else best_choice[1]


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


def _assign_local_geometry_to_row(out: pd.DataFrame, row_idx: int, geom: pd.Series, status: str) -> None:
    out.at[row_idx, "found_rt"] = float(geom["apex_x"])
    out.at[row_idx, "area"] = float(geom["area"])
    out.at[row_idx, "integration_start_x"] = float(geom["start_x"])
    out.at[row_idx, "integration_end_x"] = float(geom["end_x"])
    out.at[row_idx, "status"] = status
    out.at[row_idx, "match_score"] = np.nan
    out.at[row_idx, "matched_peak_id"] = np.nan


def _assign_local_geometry_bounds_to_row(out: pd.DataFrame, row_idx: int, geom: pd.Series, status: str) -> None:
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
) -> np.ndarray:
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
) -> int:
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

    def score(idx: int) -> float:
        value_score = float(metric[idx]) / scale
        distance_score = 0.12 * abs(float(idx) - target_idx) / span
        return value_score + distance_score

    return int(min(local_candidates, key=score))


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


def _reintegrate_cluster_by_local_minima(
    df: pd.DataFrame,
    matched_targets: pd.DataFrame,
    cluster_codes,
    window_left: float,
    window_right: float,
    status_suffix: str,
    force: bool = False,
) -> pd.DataFrame:
    out = matched_targets.copy()
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


def _should_force_c18_valley_split(matched_targets: pd.DataFrame) -> bool:
    cluster = matched_targets[matched_targets["code"].isin(["C18:2N6C", "C18:1N9C", "C18:3N3", "C18:0"])].copy()
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
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets.copy()
    if df is None or df.empty or out is None or out.empty:
        return out

    out = _reintegrate_cluster_by_local_minima(
        df=df,
        matched_targets=out,
        cluster_codes=["C18:2N6C", "C18:1N9C", "C18:3N3", "C18:0"],
        window_left=7.56,
        window_right=7.79,
        status_suffix="valley",
        force=EXPERIMENTAL_FORCE_CLUSTER_VALLEYS or _should_force_c18_valley_split(out),
    )
    out = _reintegrate_cluster_by_local_minima(
        df=df,
        matched_targets=out,
        cluster_codes=["C20:4N6", "C20:5", "C20:3N8"],
        window_left=8.34,
        window_right=8.50,
        status_suffix="valley",
        force=EXPERIMENTAL_FORCE_CLUSTER_VALLEYS,
    )
    out = _reintegrate_cluster_by_local_minima(
        df=df,
        matched_targets=out,
        cluster_codes=["C22:6", "C22:5", "C22:4"],
        window_left=9.22,
        window_right=9.34,
        status_suffix="valley",
        force=EXPERIMENTAL_FORCE_CLUSTER_VALLEYS,
    )
    return _recompute_matched_percent_area(out)


def refine_overlapped_c22_cluster_areas(
    df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets.copy()
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

    return _recompute_matched_percent_area(out)


def refine_c18_c20_cluster_matches(
    df: pd.DataFrame,
    peaks_df: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    out = matched_targets.copy()
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
        max_distances=[0.025, 0.018, 0.025],
        min_apex_gaps=[0.016, 0.020],
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


def refine_small_peak_integrations(
    df: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets.copy()
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
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets.copy()
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
    dha_width = float(
        ordered.loc[ordered["code"] == "C22:6", "integration_end_x"].iloc[0]
        - ordered.loc[ordered["code"] == "C22:6", "integration_start_x"].iloc[0]
    )
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


def _find_signal_floor_boundary(
    signal: np.ndarray,
    start_idx: int,
    limit_idx: int,
    direction: int,
    threshold: float,
    consecutive_points: int = FINAL_BOUNDARY_CONSECUTIVE_POINTS,
) -> int:
    if direction < 0:
        start_idx = int(max(start_idx, limit_idx))
        for idx in range(start_idx, int(limit_idx) - 1, -1):
            seg_start = max(int(limit_idx), idx - int(consecutive_points) + 1)
            if np.all(signal[seg_start:idx + 1] <= threshold):
                return int(idx)
        return int(limit_idx)

    start_idx = int(min(start_idx, limit_idx))
    for idx in range(start_idx, int(limit_idx) + 1):
        seg_end = min(int(limit_idx) + 1, idx + int(consecutive_points))
        if np.all(signal[idx:seg_end] <= threshold):
            return int(idx)
    return int(limit_idx)


def _omega_strict_spread(omega: dict) -> float:
    final = float(omega.get("omega3_trio", np.nan))
    strict = float(omega.get("omega3_trio_strict", np.nan))
    if not (np.isfinite(final) and np.isfinite(strict)):
        return np.nan
    return abs(final - strict)


def _with_judge_decisions(frame: pd.DataFrame, decisions: list[dict]) -> pd.DataFrame:
    out = frame.copy()
    existing = list(getattr(frame, "attrs", {}).get(JUDGE_DECISIONS_ATTR, []))
    out.attrs[JUDGE_DECISIONS_ATTR] = existing + list(decisions)
    return out


def _apply_final_boundary_candidate(frame: pd.DataFrame, item: dict) -> pd.DataFrame:
    out = frame.copy()
    row_idx = int(item["row_idx"])
    out.at[row_idx, "area"] = float(item["new_area"])
    out.at[row_idx, "integration_start_x"] = float(item["new_start_x"])
    out.at[row_idx, "integration_end_x"] = float(item["new_end_x"])
    out.at[row_idx, "status"] = _append_status_suffix(out.at[row_idx, "status"], "baseexpand")
    return _recompute_matched_percent_area(out)


def _judge_final_boundary_candidate(
    current: pd.DataFrame,
    candidate: pd.DataFrame,
    changed_count: int,
) -> tuple[bool, str, dict]:
    details = {
        "changed_count": int(changed_count),
        "current_omega": np.nan,
        "candidate_omega": np.nan,
        "omega_shift": np.nan,
        "current_strict_spread": np.nan,
        "candidate_strict_spread": np.nan,
        "strict_spread_shift": np.nan,
    }
    if changed_count <= 0:
        return False, "no_changed_peaks", details
    if changed_count > FINAL_BOUNDARY_MAX_CHANGED_PEAKS:
        return False, "too_many_changed_peaks", details

    current_omega = metrics.compute_omega(current)
    candidate_omega = metrics.compute_omega(candidate)
    current_value = float(current_omega.get("omega3_trio", np.nan))
    candidate_value = float(candidate_omega.get("omega3_trio", np.nan))
    details["current_omega"] = current_value
    details["candidate_omega"] = candidate_value
    if np.isfinite(current_value) and np.isfinite(candidate_value):
        omega_shift = candidate_value - current_value
        details["omega_shift"] = float(omega_shift)
        if abs(omega_shift) > FINAL_BOUNDARY_MAX_OMEGA_SHIFT:
            return False, "omega_shift_too_large", details

    current_spread = _omega_strict_spread(current_omega)
    candidate_spread = _omega_strict_spread(candidate_omega)
    details["current_strict_spread"] = current_spread
    details["candidate_strict_spread"] = candidate_spread
    if np.isfinite(current_spread) and np.isfinite(candidate_spread):
        spread_shift = candidate_spread - current_spread
        details["strict_spread_shift"] = float(spread_shift)
        if spread_shift > FINAL_BOUNDARY_MAX_STRICT_SPREAD_INCREASE:
            return False, "strict_spread_increase_too_large", details

    return True, "accepted", details


def _judge_fallback_boundary_candidate(
    current: pd.DataFrame,
    candidate: pd.DataFrame,
    changed_count: int,
    area_ratio: float,
) -> tuple[bool, str, dict]:
    accepted, reason, details = _judge_final_boundary_candidate(current, candidate, changed_count)
    if not accepted:
        return accepted, reason, details
    if changed_count > FINAL_BOUNDARY_FALLBACK_MAX_CHANGED_PEAKS:
        return False, "fallback_too_many_changed_peaks", details
    if area_ratio > FINAL_BOUNDARY_FALLBACK_MAX_AREA_RATIO:
        return False, "fallback_area_ratio_too_large", details
    omega_shift = float(details.get("omega_shift", np.nan))
    if np.isfinite(omega_shift) and abs(omega_shift) > FINAL_BOUNDARY_FALLBACK_MAX_OMEGA_SHIFT:
        return False, "fallback_omega_shift_too_large", details
    spread_shift = float(details.get("strict_spread_shift", np.nan))
    if np.isfinite(spread_shift) and spread_shift > FINAL_BOUNDARY_FALLBACK_MAX_STRICT_SPREAD_INCREASE:
        return False, "fallback_strict_spread_increase_too_large", details
    return True, "accepted", details


def expand_final_peak_boundaries(
    df: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    out = matched_targets.copy()
    if (
        not ENABLE_FINAL_BOUNDARY_EXPANSION
        or df is None or df.empty
        or out is None or out.empty
    ):
        return out

    x_col = _get_x_column_name(df)
    x = df[x_col].to_numpy(dtype=float)
    y_raw = df["y_corrected"].to_numpy(dtype=float)
    y = np.clip(y_raw, 0.0, None)
    if len(x) < 8 or not np.any(y > 0):
        return out

    work = out.copy()
    for column in ["found_rt", "area", "integration_start_x", "integration_end_x"]:
        work[column] = pd.to_numeric(work.get(column), errors="coerce")
    work = work.dropna(subset=["found_rt", "area", "integration_start_x", "integration_end_x"]).sort_values("found_rt")
    if work.empty:
        return out

    ordered_indices = list(work.index)
    ordered_rts = work["found_rt"].to_numpy(dtype=float)
    global_noise = max(_robust_sigma(y_raw), 1.0)

    candidate_rows: list[dict] = []
    for pos, row_idx in enumerate(ordered_indices):
        row = work.loc[row_idx]
        code = str(row.get("code", ""))
        rt = float(row["found_rt"])
        current_start = float(row["integration_start_x"])
        current_end = float(row["integration_end_x"])
        current_area = float(row["area"])
        if current_area <= 0 or current_end <= current_start:
            continue

        apex_idx = int(np.argmin(np.abs(x - rt)))
        current_start_idx = int(np.searchsorted(x, current_start, side="left"))
        current_end_idx = int(np.searchsorted(x, current_end, side="right") - 1)
        current_start_idx = max(0, min(current_start_idx, len(x) - 2))
        current_end_idx = max(current_start_idx + 1, min(current_end_idx, len(x) - 1))

        left_limit_x = rt - FINAL_BOUNDARY_MAX_EXTENSION
        right_limit_x = rt + FINAL_BOUNDARY_MAX_EXTENSION
        if pos > 0 and np.isfinite(ordered_rts[pos - 1]):
            left_limit_x = max(left_limit_x, 0.5 * (ordered_rts[pos - 1] + rt))
        if pos < len(ordered_rts) - 1 and np.isfinite(ordered_rts[pos + 1]):
            right_limit_x = min(right_limit_x, 0.5 * (rt + ordered_rts[pos + 1]))

        left_limit_idx = int(np.searchsorted(x, left_limit_x, side="left"))
        right_limit_idx = int(np.searchsorted(x, right_limit_x, side="right") - 1)
        left_limit_idx = max(0, min(left_limit_idx, apex_idx))
        right_limit_idx = max(apex_idx, min(right_limit_idx, len(x) - 1))
        if left_limit_idx >= apex_idx or right_limit_idx <= apex_idx:
            continue

        local_left = max(0, left_limit_idx - 16)
        local_right = min(len(x) - 1, right_limit_idx + 16)
        local_noise = max(_robust_sigma(y_raw[local_left:local_right + 1]), global_noise * 0.20, 1.0)
        apex_height = max(float(y[apex_idx]), float(np.nanmax(y[left_limit_idx:right_limit_idx + 1])), 1.0)
        threshold = max(local_noise * FINAL_BOUNDARY_THRESHOLD_SIGMA, apex_height * FINAL_BOUNDARY_THRESHOLD_FRACTION)

        new_start_idx = _find_signal_floor_boundary(
            signal=y,
            start_idx=current_start_idx,
            limit_idx=left_limit_idx,
            direction=-1,
            threshold=threshold,
        )
        new_end_idx = _find_signal_floor_boundary(
            signal=y,
            start_idx=current_end_idx,
            limit_idx=right_limit_idx,
            direction=1,
            threshold=threshold,
        )
        if new_start_idx >= apex_idx or new_end_idx <= apex_idx or new_end_idx <= new_start_idx:
            continue

        new_width = float(x[new_end_idx] - x[new_start_idx])
        if new_width > FINAL_BOUNDARY_MAX_WIDTH:
            continue

        new_area = float(np.trapezoid(y[new_start_idx:new_end_idx + 1], x[new_start_idx:new_end_idx + 1]))
        if current_area <= 0:
            continue
        area_ratio = new_area / current_area
        if not (FINAL_BOUNDARY_MIN_AREA_RATIO <= area_ratio <= FINAL_BOUNDARY_MAX_AREA_RATIO):
            continue

        if x[new_start_idx] >= current_start and x[new_end_idx] <= current_end:
            continue

        candidate_rows.append({
            "judge": "final_boundary_v0",
            "candidate": "baseexpand",
            "row_idx": int(row_idx),
            "code": code,
            "found_rt": rt,
            "old_start_x": current_start,
            "old_end_x": current_end,
            "new_start_x": float(x[new_start_idx]),
            "new_end_x": float(x[new_end_idx]),
            "old_area": current_area,
            "new_area": new_area,
            "area_ratio": area_ratio,
            "old_width": float(current_end - current_start),
            "new_width": new_width,
            "width_ratio": new_width / max(float(current_end - current_start), 1e-9),
            "support_threshold": threshold,
        })

    if not candidate_rows:
        return out

    full_candidate = _recompute_matched_percent_area(out)
    for item in candidate_rows:
        full_candidate = _apply_final_boundary_candidate(full_candidate, item)
    full_accepted, full_reason, full_details = _judge_final_boundary_candidate(
        matched_targets,
        full_candidate,
        len(candidate_rows),
    )
    decisions = []
    if full_accepted:
        for item in candidate_rows:
            decision = dict(item)
            decision.update(full_details)
            decision["decision"] = "accepted"
            decision["reason"] = full_reason
            decisions.append(decision)
        return _with_judge_decisions(full_candidate, decisions)

    accepted_frame = _recompute_matched_percent_area(out)
    accepted_count = 0
    for item in sorted(candidate_rows, key=lambda value: (float(value["area_ratio"]), str(value["code"]))):
        trial = _apply_final_boundary_candidate(accepted_frame, item)
        accepted, reason, judge_details = _judge_fallback_boundary_candidate(
            matched_targets,
            trial,
            accepted_count + 1,
            float(item["area_ratio"]),
        )
        decision = dict(item)
        decision.update(judge_details)
        decision["decision"] = "accepted" if accepted else "rejected"
        decision["reason"] = reason
        decisions.append(decision)
        if accepted:
            accepted_frame = trial
            accepted_count += 1

    return _with_judge_decisions(accepted_frame, decisions)



def enforce_target_rt_corridors(
    df: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    """Clip integration intervals to RT corridors derived from target order.

    Boundary expansions and local fits can occasionally let one target absorb the
    shoulder or middle of a neighboring target. For the stable C18/C20/C22 regions,
    the corrected target RTs are much more stable than the discovered integration
    width, so use midpoints between neighboring target RTs as hard split guards.
    """
    out = matched_targets.copy()
    if (
        not ENABLE_TARGET_RT_CORRIDOR_GUARD
        or df is None or df.empty
        or out is None or out.empty
    ):
        return out

    x_col = _get_x_column_name(df)
    x = df[x_col].to_numpy(dtype=float)
    y = np.clip(df["y_corrected"].to_numpy(dtype=float), 0.0, None)
    if len(x) < 4:
        return out

    work = out.copy()
    for column in ["found_rt", "corrected_target_rt", "integration_start_x", "integration_end_x", "area"]:
        work[column] = pd.to_numeric(work.get(column), errors="coerce")
    anchor_coefficient = rt_profile.estimate_anchor_coefficient(work)
    manual_centers = work["code"].map(lambda code: rt_profile.MANUAL_TABLE_RTS.get(str(code)))
    manual_centers = pd.to_numeric(manual_centers, errors="coerce")
    if np.isfinite(anchor_coefficient) and anchor_coefficient > 0:
        manual_centers = manual_centers / float(anchor_coefficient)
    work["_corridor_center"] = manual_centers.where(
        manual_centers.notna(),
        work["corrected_target_rt"].where(work["corrected_target_rt"].notna(), work["found_rt"]),
    )
    work = work.dropna(subset=["_corridor_center", "found_rt", "integration_start_x", "integration_end_x"]).sort_values("_corridor_center")
    if len(work) < 2:
        return out

    ordered_indices = list(work.index)
    centers = work["_corridor_center"].to_numpy(dtype=float)
    decisions: list[dict] = []
    for pos, row_idx in enumerate(ordered_indices):
        code = str(work.at[row_idx, "code"])
        if code not in TARGET_RT_CORRIDOR_GUARD_CODES:
            continue
        found_rt = float(work.at[row_idx, "found_rt"])
        start_x = float(work.at[row_idx, "integration_start_x"])
        end_x = float(work.at[row_idx, "integration_end_x"])
        if not (np.isfinite(found_rt) and np.isfinite(start_x) and np.isfinite(end_x) and end_x > start_x):
            continue

        left_guard = -np.inf
        right_guard = np.inf
        if pos > 0 and np.isfinite(centers[pos - 1]):
            left_guard = 0.5 * (centers[pos - 1] + centers[pos])
        if pos < len(centers) - 1 and np.isfinite(centers[pos + 1]):
            right_guard = 0.5 * (centers[pos] + centers[pos + 1])
        if not (left_guard < found_rt < right_guard):
            continue

        new_start_x = max(start_x, left_guard) if np.isfinite(left_guard) else start_x
        new_end_x = min(end_x, right_guard) if np.isfinite(right_guard) else end_x
        if new_end_x - new_start_x < TARGET_RT_CORRIDOR_MIN_WIDTH or not (new_start_x < found_rt < new_end_x):
            continue
        if new_start_x <= start_x + 1e-9 and new_end_x >= end_x - 1e-9:
            continue

        start_idx = int(np.searchsorted(x, new_start_x, side="left"))
        end_idx = int(np.searchsorted(x, new_end_x, side="right") - 1)
        start_idx = max(0, min(start_idx, len(x) - 2))
        end_idx = max(start_idx + 1, min(end_idx, len(x) - 1))
        new_area = float(np.trapezoid(y[start_idx:end_idx + 1], x[start_idx:end_idx + 1]))
        if not np.isfinite(new_area) or new_area <= 0:
            continue

        old_area = float(work.at[row_idx, "area"]) if np.isfinite(work.at[row_idx, "area"]) else np.nan
        out.at[row_idx, "integration_start_x"] = float(x[start_idx])
        out.at[row_idx, "integration_end_x"] = float(x[end_idx])
        out.at[row_idx, "area"] = new_area
        out.at[row_idx, "status"] = _append_status_suffix(out.at[row_idx, "status"], "rtcorridor")
        decisions.append({
            "judge": "target_rt_corridor_guard",
            "candidate": "rtcorridor",
            "decision": "accepted",
            "reason": "clipped_to_target_midpoint_corridor",
            "row_idx": int(row_idx),
            "code": code,
            "found_rt": found_rt,
            "old_start_x": start_x,
            "old_end_x": end_x,
            "new_start_x": float(x[start_idx]),
            "new_end_x": float(x[end_idx]),
            "old_area": old_area,
            "new_area": new_area,
            "area_ratio": new_area / max(old_area, 1e-9) if np.isfinite(old_area) else np.nan,
        })

    if not decisions:
        return out
    out = _recompute_matched_percent_area(out)
    return _with_judge_decisions(out, decisions)


DPA_OVERINTEGRATION_RATIO_MIN = 1.35
DPA_TIGHT_WIDTH_MAX = 0.028
DPA_TIGHT_HALF_WINDOW = 0.014


def tighten_dpa_overintegration_by_local_bounds(
    df: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    """Trim over-wide C22:5/DPA integrations when they dominate C22:4.

    The recurrent field failure is a DPA interval absorbing too much neighbor
    shoulder.  When DPA is much larger than C22:4 and its current interval is
    wider than a normal resolved DPA peak, re-bound it around the apex by local
    minima in a narrow window.  This changes integration boundaries directly;
    metrics-level debit remains only a secondary safety net.
    """
    out = matched_targets.copy()
    if df is None or df.empty or out is None or out.empty:
        return out
    required = {"C22:5", "C22:4"}
    if not required.issubset(set(out.get("code", []))):
        return out

    dpa_idx_list = out.index[out["code"] == "C22:5"].tolist()
    c224_idx_list = out.index[out["code"] == "C22:4"].tolist()
    if not dpa_idx_list or not c224_idx_list:
        return out
    dpa_idx = int(dpa_idx_list[0])
    c224_idx = int(c224_idx_list[0])
    dpa_area = pd.to_numeric(pd.Series([out.at[dpa_idx, "area"]]), errors="coerce").iloc[0]
    c224_area = pd.to_numeric(pd.Series([out.at[c224_idx, "area"]]), errors="coerce").iloc[0]
    found_rt = pd.to_numeric(pd.Series([out.at[dpa_idx, "found_rt"]]), errors="coerce").iloc[0]
    start_x = pd.to_numeric(pd.Series([out.at[dpa_idx, "integration_start_x"]]), errors="coerce").iloc[0]
    end_x = pd.to_numeric(pd.Series([out.at[dpa_idx, "integration_end_x"]]), errors="coerce").iloc[0]
    if not all(np.isfinite(v) for v in [dpa_area, c224_area, found_rt, start_x, end_x]):
        return out
    if c224_area <= 0 or dpa_area / c224_area <= DPA_OVERINTEGRATION_RATIO_MIN:
        return out
    if end_x - start_x <= DPA_TIGHT_WIDTH_MAX:
        return out

    x_col = _get_x_column_name(df)
    x = df[x_col].to_numpy(dtype=float)
    y = np.clip(df["y_corrected"].to_numpy(dtype=float), 0.0, None)
    if len(x) < 4:
        return out
    apex_idx = int(np.argmin(np.abs(x - float(found_rt))))
    left_limit = max(float(start_x), float(found_rt) - DPA_TIGHT_HALF_WINDOW)
    right_limit = min(float(end_x), float(found_rt) + DPA_TIGHT_HALF_WINDOW)
    left_idx = int(np.searchsorted(x, left_limit, side="left"))
    right_idx = int(np.searchsorted(x, right_limit, side="right") - 1)
    left_idx = max(0, min(left_idx, apex_idx))
    right_idx = min(len(x) - 1, max(right_idx, apex_idx))
    if right_idx <= left_idx or apex_idx <= left_idx or apex_idx >= right_idx:
        return out
    new_start_idx = int(left_idx + np.argmin(y[left_idx:apex_idx + 1]))
    new_end_idx = int(apex_idx + np.argmin(y[apex_idx:right_idx + 1]))
    if new_end_idx <= new_start_idx or not (new_start_idx < apex_idx < new_end_idx):
        return out
    new_width = float(x[new_end_idx] - x[new_start_idx])
    if new_width < TARGET_RT_CORRIDOR_MIN_WIDTH or new_width >= float(end_x - start_x):
        return out
    new_area = float(np.trapezoid(y[new_start_idx:new_end_idx + 1], x[new_start_idx:new_end_idx + 1]))
    if not np.isfinite(new_area) or new_area <= 0 or new_area >= dpa_area:
        return out
    out.at[dpa_idx, "area"] = new_area
    out.at[dpa_idx, "integration_start_x"] = float(x[new_start_idx])
    out.at[dpa_idx, "integration_end_x"] = float(x[new_end_idx])
    out.at[dpa_idx, "status"] = _append_status_suffix(out.at[dpa_idx, "status"], "dpatight")
    return _recompute_matched_percent_area(out)

def refine_cluster_matches(
    processed: pd.DataFrame,
    peaks: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    peaks, matched_targets = refine_c18_c20_cluster_matches(processed, peaks, matched_targets)
    matched_targets = refine_overlapped_c22_cluster_areas(processed, peaks, matched_targets)
    matched_targets = refine_cluster_areas_by_local_valleys(processed, matched_targets)
    matched_targets = legacy_fit.recover_missing_c22_components_with_fit(processed, peaks, matched_targets)
    matched_targets = legacy_fit.recover_underintegrated_c20_components_with_fit(processed, peaks, matched_targets)
    matched_targets = legacy_fit.recover_overlapped_c18_components_with_fit(processed, peaks, matched_targets)
    matched_targets = tighten_overwide_c22_cluster_tails(processed, matched_targets)
    matched_targets = legacy_fit.refine_overwide_c22_cluster_with_pvfit(processed, peaks, matched_targets)
    matched_targets = refine_small_peak_integrations(processed, matched_targets)
    matched_targets = expand_final_peak_boundaries(processed, matched_targets)
    matched_targets = tighten_dpa_overintegration_by_local_bounds(processed, matched_targets)
    matched_targets = enforce_target_rt_corridors(processed, matched_targets)
    return peaks, matched_targets
