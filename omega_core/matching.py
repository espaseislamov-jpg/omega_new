from __future__ import annotations

import numpy as np
import pandas as pd


RELIABLE_RT_WINDOW = 0.035
RELIABLE_RT_DOMINANT_DISTANCE_MAX = 0.025
RELIABLE_RT_DOMINANT_AREA_MULTIPLIER = 8.0
RELIABLE_RT_DOMINANT_AREA_MIN_DELTA = 400.0


def estimate_rt_shift(
    expected_rts: np.ndarray,
    observed_rts: np.ndarray,
    max_shift: float = 0.15,
    step: float = 0.001,
    sigma: float = 0.025,
) -> float:
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


def _clear_match(out: pd.DataFrame, row_idx: int) -> None:
    out.at[row_idx, "found_rt"] = np.nan
    out.at[row_idx, "area"] = np.nan
    out.at[row_idx, "percent_area"] = np.nan
    out.at[row_idx, "matched_peak_id"] = np.nan
    out.at[row_idx, "match_score"] = np.nan
    out.at[row_idx, "status"] = "not_found"


def _assign_peak(out: pd.DataFrame, row_idx: int, peak_row: pd.Series, status: str, match_score: float) -> None:
    out.at[row_idx, "found_rt"] = float(peak_row["apex_x"])
    out.at[row_idx, "area"] = float(peak_row["area"])
    out.at[row_idx, "percent_area"] = float(peak_row["percent_area"])
    out.at[row_idx, "matched_peak_id"] = int(peak_row["peak_id"])
    out.at[row_idx, "match_score"] = float(match_score)
    out.at[row_idx, "status"] = status
    out.at[row_idx, "integration_start_x"] = float(peak_row["start_x"])
    out.at[row_idx, "integration_end_x"] = float(peak_row["end_x"])


def _get_order_bounds(out: pd.DataFrame, row_idx: int) -> tuple[float, float]:
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
    matched_targets: pd.DataFrame,
    peaks: pd.DataFrame,
    cluster_codes,
    target_apexes,
    max_distance: float,
    status: str,
    rt_shift: float,
) -> pd.DataFrame:
    out = matched_targets.copy()
    if out.empty or peaks is None or peaks.empty:
        return out

    peaks_apex = peaks["apex_x"].to_numpy(dtype=float)
    peaks_prominence = peaks["prominence"].to_numpy(dtype=float)
    peaks_area = peaks["area"].to_numpy(dtype=float)
    used_mask = np.zeros(len(peaks), dtype=bool)
    chosen = {}
    for code, target_apex in zip(cluster_codes, target_apexes):
        adjusted_target = float(target_apex + rt_shift)
        distances = np.abs(peaks_apex - adjusted_target)
        candidate_positions = np.flatnonzero((~used_mask) & (distances <= max_distance))
        if candidate_positions.size == 0:
            continue
        best_pos = _select_best_peak_position(candidate_positions, distances, peaks_prominence, peaks_area)
        peak_row = peaks.iloc[best_pos]
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
    matched_targets: pd.DataFrame,
    peaks: pd.DataFrame,
    rt_shift: float = 0.0,
) -> pd.DataFrame:
    out = matched_targets.copy()
    if out.empty or peaks is None or peaks.empty:
        return out

    peaks_apex = peaks["apex_x"].to_numpy(dtype=float)
    peaks_prominence = peaks["prominence"].to_numpy(dtype=float)
    peaks_area = peaks["area"].to_numpy(dtype=float)
    cluster_codes = ["C22:6", "C22:5", "C22:4"]
    target_apexes = [9.252, 9.285, 9.316]
    max_distances = [0.020, 0.020, 0.020]

    chosen = {}
    used_mask = np.zeros(len(peaks), dtype=bool)
    for code, target_apex, max_distance in zip(cluster_codes, target_apexes, max_distances):
        adjusted_target = float(target_apex + rt_shift)
        distances = np.abs(peaks_apex - adjusted_target)
        candidate_positions = np.flatnonzero((~used_mask) & (distances <= max_distance))
        if candidate_positions.size == 0:
            continue
        best_pos = _select_best_peak_position(candidate_positions, distances, peaks_prominence, peaks_area)
        peak_row = peaks.iloc[best_pos]
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


def match_targets_to_peaks(targets_df: pd.DataFrame, peaks_df: pd.DataFrame) -> tuple[pd.DataFrame, float]:
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
