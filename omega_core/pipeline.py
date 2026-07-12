from __future__ import annotations

from pathlib import Path
import os
import warnings

import numpy as np
import pandas as pd

from . import chromatopy_adapter, clusters, io, matching, metrics, rt_profile, signal

FULL_CHROMATOPY_ENGINE = os.environ.get("OMEGA_ENGINE", "current").strip().lower() in {"chromatopy", "chromatopy_clean", "full_chromatopy"}


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



def _prepare_chromatopy_matched_targets(matched: pd.DataFrame) -> pd.DataFrame:
    out = matched.copy()
    if "target_rt" in out and "corrected_target_rt" not in out:
        out["corrected_target_rt"] = out["target_rt"]
    if "matched_peak_id" not in out:
        out["matched_peak_id"] = np.nan
    if "match_score" not in out:
        out["match_score"] = np.nan
    for peak_id, row_idx in enumerate(out.index[out["found_rt"].notna()], start=1):
        out.at[row_idx, "matched_peak_id"] = peak_id
        out.at[row_idx, "match_score"] = 0.0
    return rt_profile.annotate_rt_profile(out)


def process_chromatopy_batch(dataframe: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    import omega_chromatopy_clean

    config = omega_chromatopy_clean.IntegrationConfig(use_chromatopy_fit=False)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        clean_result = omega_chromatopy_clean.integrate_batch(dataframe, reference_targets, config)

    processed = clean_result["processed_df"]
    processed = signal.add_smoothing_and_derivatives(processed)[0] if "dy" not in processed else processed
    matched_targets = _prepare_chromatopy_matched_targets(clean_result["matched_targets_df"])
    omega = metrics.compute_omega(matched_targets)
    peaks = pd.DataFrame({
        "peak_id": pd.to_numeric(matched_targets["matched_peak_id"], errors="coerce"),
        "apex_x": pd.to_numeric(matched_targets["found_rt"], errors="coerce"),
        "area": pd.to_numeric(matched_targets["area"], errors="coerce"),
    }).dropna(subset=["peak_id", "apex_x"]).reset_index(drop=True)
    cluster_quality = metrics.compute_cluster_quality(matched_targets)
    confidence = metrics.assess_confidence(matched_targets, peaks, omega, "chromatopy_clean", cluster_quality)
    return {
        "engine": "chromatopy_clean",
        "processed_df": processed,
        "best_window": config.smoothing_window,
        "peaks_df": peaks,
        "matched_targets_df": matched_targets,
        "judge_decisions_df": pd.DataFrame(),
        "rt_shift": clean_result.get("rt_shift", 0.0),
        "omega": omega,
        "omega_report": omega["omega3_trio"],
        "boundary_mode": clean_result.get("boundary_mode", "chromatopy"),
        "cluster_quality_score": cluster_quality,
        "confidence": confidence,
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
    current_omega = float(current.get("omega_report", np.nan))
    candidate_omega = float(candidate.get("omega_report", np.nan))

    if not (np.isfinite(current_width) and np.isfinite(candidate_width)):
        return False
    if (
        np.isfinite(current_omega)
        and np.isfinite(candidate_omega)
        and abs(candidate_omega - current_omega) > 0.45
    ):
        return False
    candidate_omega_data = candidate.get("omega", {})
    candidate_c22_ratio = float(candidate_omega_data.get("c22_reference_ratio", np.nan))
    candidate_c22_debit = float(candidate_omega_data.get("c22_overintegration_debit_points", 0.0))
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


def process_batch(dataframe: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    if FULL_CHROMATOPY_ENGINE:
        try:
            return process_chromatopy_batch(dataframe, reference_targets)
        except Exception:
            if os.environ.get("OMEGA_REQUIRE_CHROMATOPY", "0").strip() == "1":
                raise

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
