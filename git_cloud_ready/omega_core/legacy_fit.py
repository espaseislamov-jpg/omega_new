from __future__ import annotations

import numpy as np
import pandas as pd


C20_FIT_EPA_AREA_MAX = 450.0
C20_FIT_EPA_PROMINENCE_MAX = 1500.0
C18_OVERLAP_START_TOLERANCE = 0.001
C22_PVFIT_OVERWIDE_MEAN_WIDTH_MIN = 0.032
C22_PVFIT_OVERWIDE_DHA_WIDTH_MIN = 0.039
C22_PVFIT_OVERWIDE_C22_4_WIDTH_MIN = 0.033


def _safe_copy(frame: pd.DataFrame | None) -> pd.DataFrame | None:
    return frame.copy() if frame is not None else None


def _cluster_has_integration_overlap(matched_targets: pd.DataFrame, cluster_codes: list[str]) -> bool:
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


def _cluster_has_duplicate_peak_ids(matched_targets: pd.DataFrame, cluster_codes: list[str]) -> bool:
    cluster = matched_targets[matched_targets["code"].isin(cluster_codes)].copy()
    if cluster.empty:
        return False
    peak_ids = pd.to_numeric(cluster["matched_peak_id"], errors="coerce").dropna().astype(int)
    return bool(peak_ids.duplicated().any())


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


def _has_complete_c22_for_missing_recovery(matched_targets: pd.DataFrame | None) -> bool:
    if matched_targets is None or matched_targets.empty:
        return False

    required_codes = ["C22:6", "C22:5", "C22:4"]
    cluster = matched_targets[matched_targets["code"].isin(required_codes)].copy()
    if len(cluster) != len(required_codes):
        return False

    cluster["area"] = pd.to_numeric(cluster.get("area"), errors="coerce")
    area_by_code = cluster.set_index("code")["area"]
    dpa_area = float(area_by_code.get("C22:5", np.nan))

    # The legacy routine is only needed when DPA was genuinely missed.
    return bool(np.isfinite(dpa_area) and dpa_area > 0)


def _should_run_c20_recovery(peaks: pd.DataFrame | None, matched_targets: pd.DataFrame | None) -> bool:
    if peaks is None or peaks.empty or matched_targets is None or matched_targets.empty:
        return False

    epa_row = matched_targets[matched_targets["code"] == "C20:5"]
    if epa_row.empty:
        return False

    epa_area = pd.to_numeric(epa_row["area"], errors="coerce").iloc[0]
    matched_peak_id = pd.to_numeric(epa_row["matched_peak_id"], errors="coerce").iloc[0]
    if not np.isfinite(epa_area) or not np.isfinite(matched_peak_id):
        return False
    if float(epa_area) > C20_FIT_EPA_AREA_MAX:
        return False

    peak_row = peaks[peaks["peak_id"] == int(matched_peak_id)]
    if peak_row.empty:
        return False
    epa_prominence = float(peak_row.iloc[0].get("raw_prominence", peak_row.iloc[0]["prominence"]))
    return bool(epa_prominence <= C20_FIT_EPA_PROMINENCE_MAX)


def _should_run_c18_recovery(matched_targets: pd.DataFrame | None) -> bool:
    if matched_targets is None or matched_targets.empty:
        return False

    c18_codes = ["C18:1N9C", "C18:3N3", "C18:0"]
    return bool(
        _should_force_c18_valley_split(matched_targets)
        or _cluster_has_duplicate_peak_ids(matched_targets, c18_codes)
        or _cluster_has_integration_overlap(matched_targets, c18_codes)
    )


def _should_run_overwide_c22_pvfit(matched_targets: pd.DataFrame | None) -> bool:
    if matched_targets is None or matched_targets.empty:
        return False

    c22_codes = ["C22:6", "C22:5", "C22:4"]
    cluster = matched_targets[matched_targets["code"].isin(c22_codes)].copy()
    if len(cluster) != len(c22_codes):
        return False

    cluster["integration_start_x"] = pd.to_numeric(cluster["integration_start_x"], errors="coerce")
    cluster["integration_end_x"] = pd.to_numeric(cluster["integration_end_x"], errors="coerce")
    cluster["area"] = pd.to_numeric(cluster["area"], errors="coerce")
    if cluster[["integration_start_x", "integration_end_x", "area"]].isna().any().any():
        return False

    status_text = " ".join(cluster["status"].fillna("").astype(str).tolist())
    if "split" in status_text or "tailtight" not in status_text:
        return False

    ordered = cluster.set_index("code").loc[c22_codes].reset_index()
    widths = (ordered["integration_end_x"] - ordered["integration_start_x"]).to_numpy(dtype=float)
    if not np.all(np.isfinite(widths)):
        return False

    mean_width = float(np.mean(widths))
    dha_width = float(widths[0])
    c22_4_width = float(widths[2])
    return bool(
        mean_width > C22_PVFIT_OVERWIDE_MEAN_WIDTH_MIN
        and dha_width > C22_PVFIT_OVERWIDE_DHA_WIDTH_MIN
        and c22_4_width > C22_PVFIT_OVERWIDE_C22_4_WIDTH_MIN
    )


def recover_missing_c22_components_with_fit(
    processed: pd.DataFrame,
    peaks: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    if _has_complete_c22_for_missing_recovery(matched_targets):
        return _safe_copy(matched_targets)

    from New_idea import recover_missing_c22_components_with_fit as legacy_func

    return legacy_func(processed, peaks, matched_targets)


def recover_underintegrated_c20_components_with_fit(
    processed: pd.DataFrame,
    peaks: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    if not _should_run_c20_recovery(peaks, matched_targets):
        return _safe_copy(matched_targets)

    from New_idea import recover_underintegrated_c20_components_with_fit as legacy_func

    return legacy_func(processed, peaks, matched_targets)


def recover_overlapped_c18_components_with_fit(
    processed: pd.DataFrame,
    peaks: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    if not _should_run_c18_recovery(matched_targets):
        return _safe_copy(matched_targets)

    from New_idea import recover_overlapped_c18_components_with_fit as legacy_func

    return legacy_func(processed, peaks, matched_targets)


def refine_overwide_c22_cluster_with_pvfit(
    processed: pd.DataFrame,
    peaks: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    if not _should_run_overwide_c22_pvfit(matched_targets):
        return _safe_copy(matched_targets)

    from New_idea import refine_overwide_c22_cluster_with_pvfit as legacy_func

    return legacy_func(processed, peaks, matched_targets)
