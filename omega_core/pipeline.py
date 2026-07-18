from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import chromatopy_adapter, clusters, io, matching, metrics, rt_profile, signal

def process_from_baseline(
    processed: pd.DataFrame,
    reference_targets: pd.DataFrame,
    strict_matching: bool = False,
) -> dict:
    processed, best_window = signal.add_smoothing_and_derivatives(processed)
    peaks = signal.detect_peak_candidates(processed, best_window=best_window)
    matched_targets, rt_shift = matching.match_targets_to_peaks(
        reference_targets,
        peaks,
        strict=strict_matching,
    )
    matched_targets = chromatopy_adapter.apply_chromatopy_target_integration(processed, matched_targets)
    peaks, matched_targets = clusters.refine_cluster_matches(processed, peaks, matched_targets)
    matched_targets = rt_profile.annotate_rt_profile(matched_targets)
    omega = metrics.compute_omega(matched_targets)
    judge_decisions = matched_targets.attrs.get(clusters.JUDGE_DECISIONS_ATTR, [])
    return {
        "processed_df": processed,
        "best_window": best_window,
        "peaks_df": peaks,
        "matched_targets_df": matched_targets,
        "judge_decisions_df": pd.DataFrame(judge_decisions),
        "rt_shift": rt_shift,
        "omega": omega,
        "omega_report": omega["omega3_trio"],
    }
def _target_width(matched_targets: pd.DataFrame, code: str) -> float:
    row = matched_targets[matched_targets["code"] == code]
    if row.empty:
        return np.nan
    start_x = pd.to_numeric(row["integration_start_x"], errors="coerce").iloc[0]
    end_x = pd.to_numeric(row["integration_end_x"], errors="coerce").iloc[0]
    if not (np.isfinite(start_x) and np.isfinite(end_x)):
        return np.nan
    return float(end_x - start_x)


def _target_area(matched_targets: pd.DataFrame, code: str) -> float:
    row = matched_targets[matched_targets["code"] == code]
    if row.empty:
        return np.nan
    return float(pd.to_numeric(row["area"], errors="coerce").iloc[0])


def _target_status(matched_targets: pd.DataFrame, code: str) -> str:
    row = matched_targets[matched_targets["code"] == code]
    if row.empty:
        return ""
    value = row.iloc[0].get("status", "")
    return "" if pd.isna(value) else str(value)


def _is_low_ratio_narrow_dha_shape(result: dict) -> bool:
    matched_targets = result.get("matched_targets_df")
    if matched_targets is None or matched_targets.empty:
        return False
    dha_width = _target_width(matched_targets, "C22:6")
    c22_4_width = _target_width(matched_targets, "C22:4")
    dpa_area = _target_area(matched_targets, "C22:5")
    c22_4_area = _target_area(matched_targets, "C22:4")
    strict_value = float(result.get("omega", {}).get("omega3_trio_strict", np.nan))
    high_omega_low_ratio = bool(
        np.all(np.isfinite([dha_width, c22_4_width, dpa_area, c22_4_area, strict_value]))
        and c22_4_area > 0
        and dha_width <= 0.030
        and c22_4_width >= 0.040
        and 0.35 <= dpa_area / c22_4_area <= 0.60
        and strict_value >= 6.0
    )
    mid_omega_balanced_ratio = bool(
        np.all(np.isfinite([dha_width, c22_4_width, dpa_area, c22_4_area, strict_value]))
        and c22_4_area > 0
        and dha_width <= 0.026
        and c22_4_width >= 0.045
        and 0.49 <= dpa_area / c22_4_area <= 0.55
        and 4.5 <= strict_value <= 5.2
    )
    return high_omega_low_ratio or mid_omega_balanced_ratio


def _should_try_asls_shape_fallback(result: dict) -> bool:
    matched_targets = result.get("matched_targets_df")
    if matched_targets is None or matched_targets.empty:
        return False

    dha_width = _target_width(matched_targets, "C22:6")
    c22_4_width = _target_width(matched_targets, "C22:4")
    c22_status = " ".join(_target_status(matched_targets, code) for code in ["C22:6", "C22:5", "C22:4"])
    strict_value = float(result.get("omega", {}).get("omega3_trio_strict", np.nan))

    broad_cluster_shape = bool(
        np.isfinite(dha_width)
        and np.isfinite(c22_4_width)
        and dha_width >= 0.040
        and c22_4_width >= 0.035
        and "matched_c22_pvfit_tail" not in c22_status
        and (not np.isfinite(strict_value) or strict_value >= 4.0)
    )
    low_ratio_narrow_dha_shape = _is_low_ratio_narrow_dha_shape(result)
    return broad_cluster_shape or low_ratio_narrow_dha_shape


def _c22_mean_width(result: dict) -> float:
    matched_targets = result.get("matched_targets_df")
    if matched_targets is None or matched_targets.empty:
        return np.nan
    widths = [_target_width(matched_targets, code) for code in ["C22:6", "C22:5", "C22:4"]]
    finite = [value for value in widths if np.isfinite(value)]
    return float(np.mean(finite)) if finite else np.nan


def _omega_final_strict_spread(result: dict) -> float:
    omega = result.get("omega", {})
    final = float(omega.get("omega3_trio", np.nan))
    strict = float(omega.get("omega3_trio_strict", np.nan))
    if not (np.isfinite(final) and np.isfinite(strict)):
        return np.nan
    return abs(final - strict)


def _confidence_score(result: dict) -> float:
    confidence = result.get("confidence", {})
    if not isinstance(confidence, dict):
        return np.nan
    # The high-recall production judge may lower the displayed score to force a
    # review, but it must never change baseline/model selection and therefore
    # must never change the reported omega value.
    return float(confidence.get("geometry_score", confidence.get("score", np.nan)))


def _accept_asls_shape_fallback(current: dict, candidate: dict) -> bool:
    current_width = _c22_mean_width(current)
    candidate_width = _c22_mean_width(candidate)
    current_quality = float(current.get("cluster_quality_score", np.nan))
    candidate_quality = float(candidate.get("cluster_quality_score", np.nan))
    current_spread = _omega_final_strict_spread(current)
    candidate_spread = _omega_final_strict_spread(candidate)
    current_confidence = _confidence_score(current)
    candidate_confidence = _confidence_score(candidate)
    current_omega = float(current.get("omega_report", np.nan))
    candidate_omega = float(candidate.get("omega_report", np.nan))
    candidate_omega_data = candidate.get("omega", {})
    candidate_c22_ratio = float(candidate_omega_data.get("c22_reference_ratio", np.nan))
    candidate_c22_debit = float(candidate_omega_data.get("c22_overintegration_debit_points", 0.0))

    if not (np.isfinite(current_width) and np.isfinite(candidate_width)):
        return False
    if (
        np.isfinite(current_omega)
        and np.isfinite(candidate_omega)
        and abs(candidate_omega - current_omega) > 0.45
    ):
        allow_high_omega_c22_downshift = bool(
            current_omega >= 9.0
            and candidate_omega < current_omega
            and current_omega - candidate_omega <= 1.10
            and np.isfinite(candidate_c22_ratio)
            and 1.35 <= candidate_c22_ratio <= 1.90
            and candidate_c22_debit > 0
        )
        if not allow_high_omega_c22_downshift:
            return False
    if (
        np.isfinite(current_omega)
        and np.isfinite(candidate_omega)
        and candidate_omega < current_omega
        and current_omega < 8.0
        and np.isfinite(candidate_c22_ratio)
        and candidate_c22_ratio >= 1.50
        and candidate_c22_debit > 0
    ):
        return False
    if candidate_width > current_width + 0.002 and not _is_low_ratio_narrow_dha_shape(current):
        return False
    if np.isfinite(current_quality) and np.isfinite(candidate_quality) and candidate_quality < current_quality:
        return False
    if np.isfinite(candidate_spread) and np.isfinite(current_spread) and candidate_spread > max(0.70, current_spread + 0.20):
        return False
    if np.isfinite(candidate_confidence) and np.isfinite(current_confidence) and candidate_confidence < current_confidence - 20.0:
        return False
    return True


def _maybe_apply_asls_shape_fallback(
    dataframe: pd.DataFrame,
    reference_targets: pd.DataFrame,
    current_result: dict,
) -> dict:
    if not (
        signal.ENABLE_ASLS_SHAPE_FALLBACK
        and signal.Baseline is not None
        and _should_try_asls_shape_fallback(current_result)
    ):
        return current_result

    alt_processed = signal.add_asls_baseline(dataframe)
    alt_result = metrics.annotate_result(
        process_from_baseline(alt_processed, reference_targets),
        baseline_mode="asls_shape_fallback",
    )
    if _accept_asls_shape_fallback(current_result, alt_result):
        return alt_result
    return current_result


def _is_structural_stop(result: dict) -> bool:
    confidence = result.get("confidence", {})
    risk = confidence.get("high_error_risk", {}) if isinstance(confidence, dict) else {}
    reason_codes = set(risk.get("reason_codes", [])) if isinstance(risk, dict) else set()
    return bool(reason_codes & {"missing_key_peak", "duplicate_peak_assignment"})


def _strict_assignments_are_sane(result: dict) -> tuple[bool, str]:
    matched = result.get("matched_targets_df")
    if matched is None or matched.empty:
        return False, "empty_retry_result"
    work = matched.copy()
    for column in ["found_rt", "area", "matched_peak_id"]:
        work[column] = pd.to_numeric(work.get(column), errors="coerce")

    key_codes = ["C20:5", "C22:6", "C22:5"]
    for code in key_codes:
        row = work[work["code"] == code]
        if row.empty or not np.isfinite(row["area"].iloc[0]) or float(row["area"].iloc[0]) <= 0:
            return False, f"key_peak_still_missing:{code}"

    peak_ids = work["matched_peak_id"].dropna().astype(int)
    if peak_ids.duplicated().any():
        return False, "duplicate_peak_assignment_remains"

    anchor = rt_profile.estimate_anchor_coefficient(work)
    guarded_groups = [
        ["C18:2N6C", "C18:1N9C", "C18:3N3", "C18:0"],
        ["C20:4N6", "C20:5", "C20:3N8"],
        ["C22:6", "C22:5", "C22:4"],
    ]
    for codes in guarded_groups:
        rows = work[work["code"].isin(codes)].set_index("code")
        found_values = []
        for code in codes:
            if code not in rows.index:
                continue
            found_rt = float(rows.at[code, "found_rt"])
            if not np.isfinite(found_rt):
                continue
            expected = rt_profile.expected_rt(code, anchor)
            if abs(found_rt - expected) > 0.055:
                return False, f"implausible_rt:{code}"
            found_values.append(found_rt)
        if len(found_values) >= 2 and np.any(np.diff(found_values) <= 0.006):
            return False, "cluster_order_or_spacing_invalid"
    return True, "strict_assignments_valid"


def _with_structural_retry_note(result: dict, retry: dict) -> dict:
    out = dict(result)
    out["structural_retry"] = dict(retry)
    confidence = dict(out.get("confidence", {}))
    confidence["structural_retry"] = dict(retry)
    reasons = list(confidence.get("reasons", []))
    if retry.get("accepted"):
        note = "Автоматическая перепроверка заново сопоставила пики и прошла строгий контроль."
        confidence["button_text"] = "Готово после автоматической перепроверки"
    else:
        note = "Автоматическая перепроверка не исправила назначение пиков — решение оставлено оператору."
    if note not in reasons:
        reasons.insert(0, note)
    confidence["reasons"] = reasons
    out["confidence"] = confidence
    return out


def _retry_structural_stop(
    dataframe: pd.DataFrame,
    reference_targets: pd.DataFrame,
    current_result: dict,
) -> dict:
    """Run exactly one stricter full rematch after a structural STOP."""
    if not _is_structural_stop(current_result):
        return current_result

    current_mode = str(current_result.get("baseline_mode", ""))
    try:
        if not current_mode.startswith("asls") and signal.Baseline is not None:
            retry_processed = signal.add_asls_baseline(dataframe)
            retry_mode = "asls_structural_retry"
        else:
            retry_processed = signal.add_baseline(dataframe, **signal.BASELINE_KWARGS)
            retry_mode = "chebyshev_structural_retry"
        candidate = metrics.annotate_result(
            process_from_baseline(retry_processed, reference_targets, strict_matching=True),
            baseline_mode=retry_mode,
        )
        assignments_sane, reason = _strict_assignments_are_sane(candidate)
        accepted = bool(assignments_sane and not _is_structural_stop(candidate))
        retry = {
            "attempted": True,
            "accepted": accepted,
            "baseline_mode": retry_mode,
            "reason": reason,
        }
        if accepted:
            return _with_structural_retry_note(candidate, retry)
        return _with_structural_retry_note(current_result, retry)
    except Exception as exc:
        return _with_structural_retry_note(current_result, {
            "attempted": True,
            "accepted": False,
            "baseline_mode": "retry_failed",
            "reason": f"retry_exception:{type(exc).__name__}",
        })


def process_batch(dataframe: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    processed = signal.add_baseline(dataframe, **signal.BASELINE_KWARGS)
    result = metrics.annotate_result(
        process_from_baseline(processed, reference_targets),
        baseline_mode="chebyshev",
    )
    result = _maybe_apply_asls_shape_fallback(dataframe, reference_targets, result)

    if (
        signal.ENABLE_ARPLS_BASELINE_FALLBACK
        and signal.Baseline is not None
        and result["cluster_quality_score"] < signal.CLUSTER_QUALITY_COMPLETE_SCORE
    ):
        alt_processed = signal.add_arpls_baseline(dataframe)
        alt_result = metrics.annotate_result(
            process_from_baseline(alt_processed, reference_targets),
            baseline_mode="arpls_fallback",
        )
        if (
            alt_result["cluster_quality_score"] >= signal.CLUSTER_QUALITY_COMPLETE_SCORE
            and alt_result["cluster_quality_score"] > result["cluster_quality_score"]
        ):
            result = alt_result

    return _retry_structural_stop(dataframe, reference_targets, result)


def process_file(
    file_path: Path,
    reference_path: Path = io.DEFAULT_REFERENCE_PATH,
    cutoff_minutes: float = 4.0,
) -> list[dict]:
    reference_targets = io.load_reference_targets(reference_path)
    batches = io.load_batches(file_path, cutoff_minutes=cutoff_minutes)
    results = []
    for batch in batches:
        result = process_batch(batch["dataframe"], reference_targets)
        results.append({**batch, **result})
    return results
