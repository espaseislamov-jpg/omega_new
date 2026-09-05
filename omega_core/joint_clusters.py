from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy.signal import savgol_filter

from . import metrics, rt_profile
from .signal import _get_x_column_name


C22_CODES = ("C22:6", "C22:5", "C22:4")
SMOOTH_WINDOWS = (5, 7, 9, 11, 15)
MIN_CONSENSUS_FRACTION = 0.60
MAX_VALLEY_SPREAD = 0.0045
MIN_PEAK_WIDTH = 0.008
MAX_PEAK_WIDTH = 0.100
MIN_LOCAL_TOTAL_RATIO = 0.85
MAX_LOCAL_TOTAL_RATIO = 1.15
MIN_COMPONENT_AREA_RATIO = 0.60
MAX_COMPONENT_AREA_RATIO = 1.50
MAX_PARTITION_L1_SHIFT = 0.20
MAX_OMEGA_SHIFT = 0.22


def _finite(value, default=np.nan) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float(default)
    return number if np.isfinite(number) else float(default)


def _valid_windows(length: int) -> list[int]:
    return [window for window in SMOOTH_WINDOWS if window <= length and window % 2 == 1]


def _consensus_valley(
    x: np.ndarray,
    signal: np.ndarray,
    left_apex: int,
    right_apex: int,
) -> tuple[int | None, dict]:
    if right_apex <= left_apex + 4:
        return None, {"support": 0, "windows": 0, "spread": np.nan}
    segment = signal[left_apex:right_apex + 1]
    x_segment = x[left_apex:right_apex + 1]
    windows = _valid_windows(len(segment))
    candidates: list[int] = []
    for window in windows:
        smooth = savgol_filter(segment, window_length=window, polyorder=min(3, window - 2), mode="interp")
        derivative = np.gradient(smooth, x_segment)
        curvature = np.gradient(derivative, x_segment)
        positive = curvature[curvature > 0]
        curvature_floor = max(
            float(np.quantile(positive, 0.10)) * 0.05 if positive.size else 0.0,
            1e-12,
        )
        local_candidates = []
        for idx in range(2, len(segment) - 2):
            local_curvature = float(np.median(curvature[idx - 1:idx + 2]))
            if derivative[idx - 1] < 0.0 <= derivative[idx] and local_curvature > curvature_floor:
                local_candidates.append(int(idx))
        if local_candidates:
            best = min(local_candidates, key=lambda idx: float(smooth[idx]))
            candidates.append(int(left_apex + best))

    required = max(1, int(math.ceil(len(windows) * MIN_CONSENSUS_FRACTION)))
    if len(candidates) < required:
        return None, {"support": len(candidates), "windows": len(windows), "spread": np.nan}
    candidate_x = x[np.asarray(candidates, dtype=int)]
    center = float(np.median(candidate_x))
    spread = float(np.max(np.abs(candidate_x - center)))
    dx = float(np.median(np.diff(x_segment)))
    if spread > max(MAX_VALLEY_SPREAD, 2.0 * dx):
        return None, {"support": len(candidates), "windows": len(windows), "spread": spread}
    return int(np.argmin(np.abs(x - center))), {
        "support": len(candidates),
        "windows": len(windows),
        "spread": spread,
    }


def _local_apex(x: np.ndarray, signal: np.ndarray, rt: float, radius: float = 0.008) -> int:
    left = int(np.searchsorted(x, float(rt) - radius, side="left"))
    right = int(np.searchsorted(x, float(rt) + radius, side="right"))
    left = max(0, min(left, len(x) - 1))
    right = max(left + 1, min(right, len(x)))
    return int(left + np.argmax(signal[left:right]))


def refine_c22_jointly(
    processed: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> tuple[pd.DataFrame, dict]:
    """Repartition C22 as one transaction while preserving its total area."""
    if processed is None or processed.empty or matched_targets is None or matched_targets.empty:
        return matched_targets.copy(), {"accepted": False, "reason": "missing_input"}
    rows = matched_targets[matched_targets["code"].isin(C22_CODES)].copy()
    rows = rows.set_index("code").reindex(C22_CODES)
    required_columns = ("found_rt", "area", "integration_start_x", "integration_end_x")
    if rows[list(required_columns)].apply(pd.to_numeric, errors="coerce").isna().any().any():
        return matched_targets.copy(), {"accepted": False, "reason": "incomplete_c22_cluster"}

    x = processed[_get_x_column_name(processed)].to_numpy(dtype=float)
    y = processed["y_corrected"].to_numpy(dtype=float)
    y_smooth = processed.get("y_smooth", processed["y_corrected"]).to_numpy(dtype=float)
    found_rts = [_finite(rows.at[code, "found_rt"]) for code in C22_CODES]
    apexes = [_local_apex(x, y_smooth, rt) for rt in found_rts]
    if not (apexes[0] < apexes[1] < apexes[2]):
        return matched_targets.copy(), {"accepted": False, "reason": "invalid_apex_order"}

    cluster_start = int(np.argmin(np.abs(x - _finite(rows.at[C22_CODES[0], "integration_start_x"]))))
    cluster_end = int(np.argmin(np.abs(x - _finite(rows.at[C22_CODES[-1], "integration_end_x"]))))
    if not (cluster_start < apexes[0] < apexes[1] < apexes[2] < cluster_end):
        return matched_targets.copy(), {"accepted": False, "reason": "invalid_outer_bounds"}

    valley_1, evidence_1 = _consensus_valley(x, y_smooth, apexes[0], apexes[1])
    valley_2, evidence_2 = _consensus_valley(x, y_smooth, apexes[1], apexes[2])
    if valley_1 is None or valley_2 is None or not (apexes[0] < valley_1 < apexes[1] < valley_2 < apexes[2]):
        return matched_targets.copy(), {
            "accepted": False,
            "reason": "no_consensus_valleys",
            "valley_1": evidence_1,
            "valley_2": evidence_2,
        }

    boundaries = (cluster_start, valley_1, valley_2, cluster_end)
    widths = [float(x[right] - x[left]) for left, right in zip(boundaries[:-1], boundaries[1:])]
    if any(width < MIN_PEAK_WIDTH or width > MAX_PEAK_WIDTH for width in widths):
        return matched_targets.copy(), {"accepted": False, "reason": "implausible_component_width", "widths": widths}

    edge_points = max(3, min(9, (cluster_end - cluster_start) // 12))
    left_level = float(np.median(y[cluster_start:cluster_start + edge_points]))
    right_level = float(np.median(y[cluster_end - edge_points + 1:cluster_end + 1]))
    baseline = np.linspace(left_level, right_level, cluster_end - cluster_start + 1)
    local_signal = np.clip(y[cluster_start:cluster_end + 1] - baseline, 0.0, None)
    raw_areas = []
    for left, right in zip(boundaries[:-1], boundaries[1:]):
        raw_areas.append(float(np.trapezoid(
            local_signal[left - cluster_start:right - cluster_start + 1],
            x[left:right + 1],
        )))
    current_areas = np.asarray([_finite(rows.at[code, "area"]) for code in C22_CODES], dtype=float)
    current_total = float(np.sum(current_areas))
    raw_total = float(np.sum(raw_areas))
    local_total_ratio = raw_total / max(current_total, 1e-12)
    if not MIN_LOCAL_TOTAL_RATIO <= local_total_ratio <= MAX_LOCAL_TOTAL_RATIO:
        return matched_targets.copy(), {
            "accepted": False,
            "reason": "local_total_disagrees",
            "local_total_ratio": local_total_ratio,
        }
    scaled_areas = np.asarray(raw_areas, dtype=float) * (current_total / max(raw_total, 1e-12))
    component_ratios = scaled_areas / np.maximum(current_areas, 1e-12)
    if np.any(component_ratios < MIN_COMPONENT_AREA_RATIO) or np.any(component_ratios > MAX_COMPONENT_AREA_RATIO):
        return matched_targets.copy(), {
            "accepted": False,
            "reason": "component_change_too_large",
            "component_ratios": component_ratios.tolist(),
        }
    current_fractions = current_areas / current_total
    new_fractions = scaled_areas / current_total
    partition_l1_shift = float(np.sum(np.abs(new_fractions - current_fractions)))
    if partition_l1_shift > MAX_PARTITION_L1_SHIFT:
        return matched_targets.copy(), {
            "accepted": False,
            "reason": "partition_shift_too_large",
            "partition_l1_shift": partition_l1_shift,
        }

    out = matched_targets.copy()
    for code, left, right, area in zip(C22_CODES, boundaries[:-1], boundaries[1:], scaled_areas):
        row_index = out.index[out["code"] == code]
        if len(row_index) != 1:
            return matched_targets.copy(), {"accepted": False, "reason": "duplicate_c22_target"}
        idx = row_index[0]
        out.at[idx, "integration_start_x"] = float(x[left])
        out.at[idx, "integration_end_x"] = float(x[right])
        out.at[idx, "area"] = float(area)
        out.at[idx, "status"] = "joint_c22_partition"
    total_area = float(pd.to_numeric(out["area"], errors="coerce").fillna(0.0).sum())
    out["percent_area"] = 100.0 * pd.to_numeric(out["area"], errors="coerce") / total_area if total_area > 0 else np.nan
    return out, {
        "accepted": True,
        "reason": "joint_partition_ready",
        "boundaries": [float(x[index]) for index in boundaries],
        "widths": widths,
        "current_areas": current_areas.tolist(),
        "new_areas": scaled_areas.tolist(),
        "component_ratios": component_ratios.tolist(),
        "local_total_ratio": local_total_ratio,
        "partition_l1_shift": partition_l1_shift,
        "valley_1": evidence_1,
        "valley_2": evidence_2,
    }


def build_joint_c22_candidate(current_result: dict) -> tuple[dict | None, dict]:
    current_matched = current_result.get("matched_targets_df")
    if rt_profile.uses_custom_instrument_profile(current_matched):
        return None, {"accepted": False, "reason": "custom_profile_not_validated"}
    matched, decision = refine_c22_jointly(
        current_result.get("processed_df"),
        current_matched,
    )
    if not decision.get("accepted"):
        return None, decision
    matched = metrics.annotate_peak_heights(current_result["processed_df"], matched)
    omega = metrics.compute_omega(matched)
    cluster_quality = metrics.compute_cluster_quality(matched)
    confidence = metrics.assess_confidence(
        matched_targets=matched,
        peaks=current_result["peaks_df"],
        omega=omega,
        baseline_mode=str(current_result.get("baseline_mode", "chebyshev")),
        cluster_quality_score=cluster_quality,
    )
    multiplier = _finite(current_result.get("result_multiplier"), 1.0)
    omega_report = _finite(omega.get("omega3_trio")) * multiplier
    return {
        "matched_targets_df": matched,
        "omega": omega,
        "omega_report": omega_report,
        "cluster_quality_score": cluster_quality,
        "confidence": confidence,
    }, decision


def judge_joint_c22_candidate(current_result: dict, candidate: dict | None, decision: dict) -> tuple[bool, str]:
    if candidate is None or not decision.get("accepted"):
        return False, str(decision.get("reason", "no_candidate"))
    current_omega = _finite(current_result.get("omega_report"))
    candidate_omega = _finite(candidate.get("omega_report"))
    current_risk = current_result.get("confidence", {}).get("high_error_risk", {})
    if int(current_risk.get("score", 0)) < 85:
        return False, "current_not_high_risk"
    if not np.isfinite(candidate_omega) or abs(candidate_omega - current_omega) > MAX_OMEGA_SHIFT:
        return False, "omega_shift_too_large"
    candidate_risk = candidate.get("confidence", {}).get("high_error_risk", {})
    if int(candidate_risk.get("score", 0)) > int(current_risk.get("score", 0)):
        return False, "high_error_risk_increased"
    return True, "accepted"


def maybe_apply_joint_c22_retry(current_result: dict) -> dict:
    """Attempt one bounded C22 re-partition only for a legacy high-risk result."""
    current_risk = current_result.get("confidence", {}).get("high_error_risk", {})
    current_matched = current_result.get("matched_targets_df")
    if int(current_risk.get("score", 0)) < 85:
        out = dict(current_result)
        out["joint_c22_retry"] = {
            "attempted": False,
            "accepted": False,
            "reason": "current_not_high_risk",
            "current_omega": _finite(current_result.get("omega_report")),
            "candidate_omega": np.nan,
        }
        return out
    if rt_profile.uses_custom_instrument_profile(current_matched):
        out = dict(current_result)
        out["joint_c22_retry"] = {
            "attempted": False,
            "accepted": False,
            "reason": "custom_profile_not_validated",
            "current_omega": _finite(current_result.get("omega_report")),
            "candidate_omega": np.nan,
        }
        return out
    candidate, decision = build_joint_c22_candidate(current_result)
    accepted, reason = judge_joint_c22_candidate(current_result, candidate, decision)
    audit = {
        **decision,
        "attempted": True,
        "accepted": bool(accepted),
        "reason": reason,
        "current_omega": _finite(current_result.get("omega_report")),
        "candidate_omega": _finite(candidate.get("omega_report")) if candidate is not None else np.nan,
    }
    if not accepted or candidate is None:
        out = dict(current_result)
        out["joint_c22_retry"] = audit
        return out
    out = dict(current_result)
    out.update(candidate)
    out["joint_c22_retry"] = audit
    return out
