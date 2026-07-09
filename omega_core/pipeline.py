from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import chromatopy_adapter, clusters, io, matching, metrics, rt_profile, signal


def process_from_baseline(processed: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    processed, best_window = signal.add_smoothing_and_derivatives(processed)
    peaks = signal.detect_peak_candidates(processed, best_window=best_window)
    matched_targets, rt_shift = matching.match_targets_to_peaks(reference_targets, peaks)
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


def _target_status(matched_targets: pd.DataFrame, code: str) -> str:
    row = matched_targets[matched_targets["code"] == code]
    if row.empty:
        return ""
    value = row.iloc[0].get("status", "")
    return "" if pd.isna(value) else str(value)


def _should_try_asls_shape_fallback(result: dict) -> bool:
    matched_targets = result.get("matched_targets_df")
    if matched_targets is None or matched_targets.empty:
        return False

    dha_width = _target_width(matched_targets, "C22:6")
    c22_4_width = _target_width(matched_targets, "C22:4")
    c22_status = " ".join(_target_status(matched_targets, code) for code in ["C22:6", "C22:5", "C22:4"])
    strict_value = float(result.get("omega", {}).get("omega3_trio_strict", np.nan))

    return bool(
        np.isfinite(dha_width)
        and np.isfinite(c22_4_width)
        and dha_width >= 0.040
        and c22_4_width >= 0.035
        and "matched_c22_pvfit_tail" not in c22_status
        and (not np.isfinite(strict_value) or strict_value >= 4.0)
    )


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
    return float(confidence.get("score", np.nan))


def _accept_asls_shape_fallback(current: dict, candidate: dict) -> bool:
    current_width = _c22_mean_width(current)
    candidate_width = _c22_mean_width(candidate)
    current_quality = float(current.get("cluster_quality_score", np.nan))
    candidate_quality = float(candidate.get("cluster_quality_score", np.nan))
    current_spread = _omega_final_strict_spread(current)
    candidate_spread = _omega_final_strict_spread(candidate)
    current_confidence = _confidence_score(current)
    candidate_confidence = _confidence_score(candidate)

    if not (np.isfinite(current_width) and np.isfinite(candidate_width)):
        return False
    if candidate_width > current_width + 0.002:
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
            return alt_result

    return result


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
