from __future__ import annotations

import math
from collections.abc import Mapping

import numpy as np
import pandas as pd

from . import rt_profile


C22_OVERLAP_TRIGGER_OMEGA_MIN = 4.0
C22_OVERLAP_RATIO_OFFSET = 1.25
C22_OVERLAP_RATIO_SLOPE = 1.0
C22_OVERLAP_FRACTION_CAP = 0.95
ENABLE_DATA_DRIVEN_C22_OVERLAP_MODEL = True
C22_OVERLAP_MODEL_BLEND = 0.60
C22_OVERLAP_MODEL_APPLY_FRACTION_MIN = 0.90
C22_OVERLAP_WIDE_CLUSTER_MEAN_WIDTH = 0.030
C22_OVERLAP_WIDE_CLUSTER_SCALE = 0.65
# Conservative C22/DPA over-integration guard: when DPA dwarfs C22:4 in the
# same local cluster, a bounded part of DPA is treated as likely shared tail area.
C22_DPA_OVERINTEGRATION_RATIO_MIN = 1.35
C22_DPA_OVERINTEGRATION_DPA_FRACTION = 0.30
C22_DPA_OVERINTEGRATION_MAX_OMEGA_POINTS = 0.45
# The narrow-DPA/broad-C22:4 shape is a reliable resolved-cluster guard only
# inside the RT profile on which it was validated.  A later high-anchor profile
# has the same apparent widths but systematically over-assigns the shared tail
# to DPA, so it must retain the bounded over-integration debit.
C22_RESOLVED_PROFILE_ANCHOR_MIN = 0.9999
C22_RESOLVED_PROFILE_ANCHOR_MAX = 1.00065
# Bounded C22 width-balance calibration learned from regression diagnostics.
# It nudges narrow DPA/C22:4 cluster cases by at most a few tenths of an omega point.
C22_WIDTH_BALANCE_DPA_WIDTH_MAX = 0.030
C22_WIDTH_BALANCE_C22_4_NARROW_MAX = 0.020
C22_WIDTH_BALANCE_DHA_WIDTH_MAX = 0.040
C22_WIDTH_BALANCE_SMALL_DPA_AREA_MAX = 163.69
C22_WIDTH_BALANCE_DHA_AREA_MAX = 1618.65
C22_WIDTH_BALANCE_LOW_OVERLAP_FRACTION_MAX = 0.61
C22_WIDTH_BALANCE_POSITIVE_POINTS = 0.20
C22_WIDTH_BALANCE_NEGATIVE_POINTS = -0.10
# Guard against C22 credit overshoot in borderline DPA/C22:4-ratio cases.
# These are exactly the cases where an integration shoulder can look like shared
# C22:4 tail area, but the strict trio value is already near the clinical target.
C22_LOW_RATIO_CREDIT_GUARD_MIN = 0.45
C22_LOW_RATIO_CREDIT_GUARD_MAX = 0.55
C22_LOW_RATIO_CREDIT_GUARD_STRICT_MAX = 5.50
C22_MISSING_DHA_COELUTION_RATIO_MIN = 1.50
C22_MISSING_DHA_COELUTION_CREDIT_FRACTION = 0.75
C22_LOW_RATIO_CREDIT_GUARD_MAX_POINTS = 0.25

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

# Re-enabled behind a tighter ratio gate after July-regression validation.
# The broad C20/EPA model caused over-estimation, but a narrow severe-underfit
# gate fixes the remaining large EPA under-integration failures.
ENABLE_DATA_DRIVEN_C20_EPA_MODEL = True
C20_EPA_MODEL_GATE_RATIO_MAX = 0.25

# EPA credit is valid for the older RT profile where the C20:5 shoulder is
# routinely under-integrated into C20:3N8.  In the newer profile the same model
# over-corrects badly, so gate it by the stable-anchor coefficient derived from
# the manual workbook RTs.
C20_EPA_MODEL_ANCHOR_COEFFICIENT_MAX = 0.99975
C20_EPA_MODEL_BLEND = 0.25
C20_EPA_OVERLAP_WIDE_NEIGHBOR_RATIO = 1.30
C20_EPA_OVERLAP_EXTRA_SCALE = 1.60
C20_EPA_UNDERFIT_RATIO_MAX = 0.13
C20_EPA_UNDERFIT_WIDTH_RATIO = 1.80
C20_EPA_UNDERFIT_STRICT_MAX = 4.60
C20_EPA_UNDERFIT_CREDIT_MIN = 20.0
C20_EPA_UNDERFIT_EXTRA_SCALE = 1.70
# A very broad C20:3N8 component is not evidence that most of its area belongs
# to EPA.  The historical 65% floor turned exactly this geometry into the
# largest positive C20 outliers in the legacy corpus.  Keep the floor only for
# a compact neighbour; broad neighbours retain the bounded model estimate.
C20_EPA_FORCED_CREDIT_NEIGHBOR_WIDTH_MAX = 0.072
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

CLUSTER_QUALITY_COMPLETE_SCORE = 50.0
GEOMETRY_READY_SCORE = 85.0
GEOMETRY_STOP_SCORE = 55.0

# Conservative, reference-free warning envelope calibrated on every available
# non-sealed manual batch.  It is intentionally tuned for recall: all known
# abs(error) > 0.5 cases must be sent to manual review.  The sealed 14072026
# batch was not used to select these limits.
HIGH_ERROR_LEGACY_FRACTION_SPLIT = 0.0041
HIGH_ERROR_DPA_RT_MAX = 9.2777
HIGH_ERROR_C20_TARGET_RT_MIN = 7.7845
HIGH_ERROR_EPA_RT_MAX = 8.4140
HIGH_ERROR_EARLY_C20_TARGET_RT_MAX = 7.7795
HIGH_ERROR_C22_4_LEFT_WIDTH_MIN = 0.01746
# This threshold is intentionally a review-order hint, not a warning rule.
# On the non-sealed legacy corpus it ranked large errors slightly better than
# the corresponding area ratio, but replacing the area rule missed one held-out
# batch error.  A custom instrument profile must provide its own validated
# threshold before the height ratio is interpreted.
LEGACY_C22_HEIGHT_PRIORITY_RATIO = 1.3635

# Emergency, deliberately broad geometry gates for a custom instrument profile.
# They are not a substitute for profile calibration.  Their only purpose is to
# keep obviously misplaced/truncated numerator peaks out of an unattended
# result while the second-instrument reference set is still small.
EMERGENCY_KEY_RT_RESIDUAL_MAX = 0.020
EMERGENCY_KEY_WIDTH_MIN = 0.008
EMERGENCY_KEY_WIDTH_MAX = 0.085
EMERGENCY_APEX_EDGE_MIN = 0.0025
EMERGENCY_C20_HEIGHT_RATIO_MIN = 0.20
EMERGENCY_C20_HEIGHT_RATIO_MAX = 1.20
EMERGENCY_C22_HEIGHT_RATIO_MIN = 0.55
EMERGENCY_C22_HEIGHT_RATIO_MAX = 2.05


def annotate_peak_heights(
    processed: pd.DataFrame,
    matched_targets: pd.DataFrame,
) -> pd.DataFrame:
    """Attach smoothed baseline-corrected peak heights inside accepted bounds."""
    if matched_targets is None:
        return matched_targets
    out = matched_targets.copy()
    out["peak_height_smooth"] = np.nan
    if processed is None or processed.empty or out.empty:
        return out
    if not {"integration_start_x", "integration_end_x"}.issubset(out.columns):
        return out

    x_column = "x_corrected" if "x_corrected" in processed else "x"
    if x_column not in processed or "y_smooth" not in processed:
        return out
    x = pd.to_numeric(processed[x_column], errors="coerce").to_numpy(dtype=float)
    y_smooth = pd.to_numeric(processed["y_smooth"], errors="coerce").to_numpy(dtype=float)
    finite_signal = np.isfinite(x) & np.isfinite(y_smooth)
    if not finite_signal.any():
        return out

    for row_index, row in out.iterrows():
        start_x = pd.to_numeric(pd.Series([row.get("integration_start_x")]), errors="coerce").iloc[0]
        end_x = pd.to_numeric(pd.Series([row.get("integration_end_x")]), errors="coerce").iloc[0]
        if not (np.isfinite(start_x) and np.isfinite(end_x) and end_x > start_x):
            continue
        mask = finite_signal & (x >= float(start_x)) & (x <= float(end_x))
        if int(mask.sum()) < 3:
            continue
        out.at[row_index, "peak_height_smooth"] = max(float(np.nanmax(y_smooth[mask])), 0.0)
    return out


def compute_omega(matched_targets: pd.DataFrame) -> dict:
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
        "c22_overintegration_debit_area": 0.0,
        "c22_overintegration_debit_points": 0.0,
        "c22_overintegration_model_applied": False,
        "c22_width_balance_points": 0.0,
        "c22_width_balance_model_applied": False,
        "c22_missing_dha_coelution_applied": False,
        "c18_denominator_scale": 1.0,
    }
    if matched_targets is None or matched_targets.empty:
        return result

    valid = matched_targets[pd.notna(matched_targets["area"])].copy()
    if valid.empty:
        return result
    total_area = float(valid["area"].sum())
    if total_area <= 0 or not np.isfinite(total_area):
        return result

    def area_of(code: str) -> float:
        row = valid[valid["code"] == code]
        return float(row["area"].iloc[0]) if not row.empty else 0.0

    def width_of(code: str) -> float:
        row = valid[valid["code"] == code]
        if row.empty:
            return np.nan
        start_x = pd.to_numeric(row["integration_start_x"], errors="coerce").iloc[0] if "integration_start_x" in row else np.nan
        end_x = pd.to_numeric(row["integration_end_x"], errors="coerce").iloc[0] if "integration_end_x" in row else np.nan
        if not (np.isfinite(start_x) and np.isfinite(end_x)):
            return np.nan
        return float(end_x - start_x)

    def status_of(code: str) -> str:
        row = matched_targets[matched_targets["code"] == code]
        if row.empty or "status" not in row:
            return ""
        value = row["status"].iloc[0]
        return "" if pd.isna(value) else str(value)

    calculation_mode = "legacy"
    if "instrument_profile_calculation_mode" in valid:
        calculation_mode = str(valid["instrument_profile_calculation_mode"].iloc[0])
    if calculation_mode == "direct_components":
        if "instrument_profile_omega_component" not in valid:
            return result
        component_rows = valid[
            valid["instrument_profile_omega_component"].fillna(False).astype(bool)
        ]
        component_area = float(component_rows["area"].sum())
        direct_value = 100.0 * component_area / total_area
        epa = area_of("C20:5")
        dha = area_of("C22:6")
        dpa = area_of("C22:5")
        c20_3 = area_of("C20:3N8")
        c22_4 = area_of("C22:4")
        result.update({
            "omega3_trio": direct_value,
            "omega3_trio_strict": direct_value,
            "omega3_trio_corrected": direct_value,
            "total_area": total_area,
            "effective_total_area": total_area,
            "epa_area": epa,
            "dha_area": dha,
            "dpa_area": dpa,
            "epa_neighbor_area": c20_3,
            "epa_effective_area": epa,
            "c22_overlap_source_area": c22_4,
            "c22_reference_ratio": dpa / c22_4 if c22_4 > 0 else np.nan,
            "epa_anchor_coefficient": rt_profile.estimate_anchor_coefficient(valid),
            "epa_anchor_gate_max": C20_EPA_MODEL_ANCHOR_COEFFICIENT_MAX,
            "profile_component_area": component_area,
            "profile_calculation_applied": True,
        })
        return result

    epa, dha, dpa = area_of("C20:5"), area_of("C22:6"), area_of("C22:5")
    anchor_coefficient = rt_profile.estimate_anchor_coefficient(valid)
    c20_3 = area_of("C20:3N8")
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
    epa_status = status_of("C20:5")
    epa_was_deconvolved_from_missing = "recovered_c20_missing_epa_component" in epa_status
    epa_to_c20_3_ratio = epa / c20_3 if c20_3 > 0 else np.nan
    if (
        ENABLE_DATA_DRIVEN_C20_EPA_MODEL
        and not epa_was_deconvolved_from_missing
        and np.isfinite(anchor_coefficient)
        and anchor_coefficient <= C20_EPA_MODEL_ANCHOR_COEFFICIENT_MAX
        and c20_3 > 0
        and np.isfinite(epa_to_c20_3_ratio)
        and epa_to_c20_3_ratio < C20_EPA_MODEL_GATE_RATIO_MAX
    ):
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

    # Keep the fitted-EPA correction monotonic at the two well-separated ends.
    # An extremely small fitted EPA is a genuine unresolved shoulder, while a
    # moderate EPA/C20:3 ratio does not need an additional overlap credit.
    if "matched_c20_fit" in status_of("C20:5") and np.isfinite(epa_to_c20_3_ratio):
        if (
            epa_to_c20_3_ratio <= 0.10
            and np.isfinite(c20_width_epa)
            and np.isfinite(c20_width_neighbor)
            and c20_width_epa <= 0.025
            and c20_width_neighbor >= 2.0 * c20_width_epa
            and c20_width_neighbor <= C20_EPA_FORCED_CREDIT_NEIGHBOR_WIDTH_MAX
        ):
            epa_credit_area = max(epa_credit_area, 0.65 * c20_3)
            epa_overlap_fraction = epa_credit_area / c20_3 if c20_3 > 0 else 0.0
            epa_model_applied = True
        elif epa_to_c20_3_ratio >= 0.18:
            epa_credit_area = 0.0
            epa_overlap_fraction = 0.0
            epa_model_applied = False
            epa_extra_scale = 1.0

    c22_ratio = dpa / c22_4 if c22_4 > 0 else np.nan
    c22_width_values = np.asarray([width_of("C22:6"), width_of("C22:5"), width_of("C22:4")], dtype=float)
    c22_missing_dha_coelution = bool(
        dha <= 0
        and dpa > 0
        and c22_4 > 0
        and np.isfinite(c22_ratio)
        and c22_ratio >= C22_MISSING_DHA_COELUTION_RATIO_MIN
        and "not_found" in status_of("C22:6")
    )
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
    if c22_missing_dha_coelution:
        c22_fraction = max(c22_fraction, C22_MISSING_DHA_COELUTION_CREDIT_FRACTION)
    c22_credit_area = c22_4 * c22_fraction
    if (
        np.isfinite(c22_ratio)
        and C22_LOW_RATIO_CREDIT_GUARD_MIN <= c22_ratio <= C22_LOW_RATIO_CREDIT_GUARD_MAX
        and strict_value < C22_LOW_RATIO_CREDIT_GUARD_STRICT_MAX
        and c22_credit_area > 0
    ):
        max_guarded_credit = effective_total_area * C22_LOW_RATIO_CREDIT_GUARD_MAX_POINTS / 100.0
        c22_credit_area = float(min(c22_credit_area, max_guarded_credit))
        c22_fraction = c22_credit_area / c22_4 if c22_4 > 0 else 0.0

    c22_debit_area = 0.0
    c22_debit_points = 0.0
    c22_debit_applied = False
    w_dha, w_dpa, w_c22_4 = c22_width_values
    resolved_broad_c22_reference = bool(
        np.all(np.isfinite([w_dpa, w_c22_4]))
        and np.isfinite(anchor_coefficient)
        and C22_RESOLVED_PROFILE_ANCHOR_MIN <= anchor_coefficient <= C22_RESOLVED_PROFILE_ANCHOR_MAX
        and w_dpa <= 0.030
        and w_c22_4 >= 0.035
        and w_c22_4 >= w_dpa + 0.009
    )
    if (
        c22_4 > 0
        and dpa > 0
        and np.isfinite(c22_ratio)
        and c22_ratio > C22_DPA_OVERINTEGRATION_RATIO_MIN
        and not c22_missing_dha_coelution
        and not resolved_broad_c22_reference
    ):
        max_debit_area = effective_total_area * C22_DPA_OVERINTEGRATION_MAX_OMEGA_POINTS / 100.0
        c22_debit_area = float(np.clip(
            dpa * C22_DPA_OVERINTEGRATION_DPA_FRACTION,
            0.0,
            max_debit_area,
        ))
        c22_debit_points = 100.0 * c22_debit_area / effective_total_area
        c22_debit_applied = c22_debit_area > 0

    c22_width_balance_points = 0.0
    c22_width_balance_applied = False
    if np.all(np.isfinite([w_dha, w_dpa, w_c22_4])) and w_dpa <= C22_WIDTH_BALANCE_DPA_WIDTH_MAX:
        if w_c22_4 <= C22_WIDTH_BALANCE_C22_4_NARROW_MAX:
            if w_dha <= C22_WIDTH_BALANCE_DHA_WIDTH_MAX:
                c22_width_balance_points = C22_WIDTH_BALANCE_POSITIVE_POINTS
        elif dpa <= C22_WIDTH_BALANCE_SMALL_DPA_AREA_MAX:
            if c22_fraction <= C22_WIDTH_BALANCE_LOW_OVERLAP_FRACTION_MAX:
                c22_width_balance_points = C22_WIDTH_BALANCE_POSITIVE_POINTS
            else:
                c22_width_balance_points = C22_WIDTH_BALANCE_NEGATIVE_POINTS
        elif dha <= C22_WIDTH_BALANCE_DHA_AREA_MAX:
            c22_width_balance_points = C22_WIDTH_BALANCE_NEGATIVE_POINTS
    c22_width_balance_applied = abs(c22_width_balance_points) > 1e-12

    corrected_value = (
        100.0 * (epa + dha + dpa + epa_credit_area + c22_credit_area - c22_debit_area) / effective_total_area
        + c22_width_balance_points
    )

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
        "epa_anchor_coefficient": anchor_coefficient,
        "epa_anchor_gate_max": C20_EPA_MODEL_ANCHOR_COEFFICIENT_MAX,
        "c22_overlap_source_area": c22_4,
        "c22_overlap_credit_area": c22_credit_area,
        "c22_overlap_fraction": c22_fraction,
        "c22_overlap_legacy_fraction": legacy_fraction,
        "c22_overlap_model_fraction": model_fraction,
        "c22_overlap_model_applied": model_applied,
        "c22_reference_ratio": c22_ratio,
        "c22_width_scale": c22_width_scale,
        "c22_overintegration_debit_area": c22_debit_area,
        "c22_overintegration_debit_points": c22_debit_points,
        "c22_overintegration_model_applied": c22_debit_applied,
        "c22_width_balance_points": c22_width_balance_points,
        "c22_width_balance_model_applied": c22_width_balance_applied,
        "c22_missing_dha_coelution_applied": c22_missing_dha_coelution,
        "c18_denominator_scale": c18_denominator_scale,
    })
    return result


def compute_cluster_quality(matched_targets: pd.DataFrame) -> float:
    if matched_targets is None or matched_targets.empty:
        return -np.inf

    cluster_groups = [
        ["C18:2N6C", "C18:1N9C", "C18:3N3", "C18:0"],
        ["C20:4N6", "C20:5", "C20:3N8"],
        ["C22:6", "C22:5", "C22:4"],
    ]
    score = 0.0
    for cluster_codes in cluster_groups:
        cluster = matched_targets.loc[matched_targets["code"].isin(cluster_codes), ["found_rt", "area", "matched_peak_id"]]
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


def classify_high_error_risk(features: Mapping[str, object]) -> dict:
    """Classify reference-free risk of an omega error above 0.5 percentage point.

    This is deliberately a small, inspectable rule set rather than another
    opaque score.  The numeric envelope was selected on the non-sealed manual
    batches with a hard 100% recall constraint.
    """

    def finite(name: str) -> float:
        try:
            value = float(features.get(name, np.nan))
        except (TypeError, ValueError):
            return np.nan
        return value if np.isfinite(value) else np.nan

    def truthy(name: str, default: bool = False) -> bool:
        value = features.get(name, default)
        if isinstance(value, str):
            return value.strip().casefold() in {"1", "true", "yes", "да"}
        return bool(value)

    score = 0
    reason_codes: list[str] = []
    reasons: list[str] = []
    peak_codes: set[str] = set()

    custom_profile = truthy("instrument_profile_custom_rt")
    judge_calibrated = truthy(
        "instrument_profile_judge_calibrated",
        default=not custom_profile,
    )
    emergency_enabled = truthy("instrument_profile_judge_emergency_enabled")
    emergency_target = finite("instrument_profile_judge_target_abs_error")
    c22_height_ratio = finite("C22_height_ratio")
    c22_height_threshold = finite("C22_height_priority_threshold")

    if custom_profile and not judge_calibrated and not emergency_enabled:
        score = max(score, 88)
        reason_codes.append("instrument_profile_uncalibrated")
        peak_codes.update(["C20:5", "C22:6", "C22:5", "C22:4"])
        reasons.append(
            "Судья ещё не проверен на ручных результатах этого прибора. "
            "До профильной калибровки результат нельзя выпускать только по зелёному статусу."
        )

    emergency_codes = [
        str(code) for code in (features.get("emergency_geometry_codes", []) or []) if str(code)
    ]
    emergency_reasons = [
        str(reason) for reason in (features.get("emergency_geometry_reasons", []) or []) if str(reason)
    ]
    if custom_profile and emergency_enabled and emergency_codes:
        score = max(score, 88)
        reason_codes.append("emergency_profile_geometry")
        peak_codes.update(emergency_codes)
        reasons.append(
            "Экстренная проверка нашла подозрительное назначение или неполные границы "
            f"пиков: {', '.join(sorted(set(emergency_codes)))}."
        )
        reasons.extend(emergency_reasons)

    missing_key_peaks = [
        str(code) for code in (features.get("missing_key_peaks", []) or []) if str(code)
    ]
    if missing_key_peaks:
        score = max(score, 98)
        reason_codes.append("missing_key_peak")
        peak_codes.update(missing_key_peaks)
        reasons.append(
            f"Не найден важный пик: {', '.join(missing_key_peaks)}. "
            "Результат нельзя выпускать без ручной проверки."
        )

    duplicate_codes = [
        str(code) for code in (features.get("duplicate_peak_codes", []) or []) if str(code)
    ]
    if duplicate_codes:
        score = max(score, 98)
        reason_codes.append("duplicate_peak_assignment")
        peak_codes.update(duplicate_codes)
        reasons.append(
            "Один участок сигнала назначен нескольким пикам. "
            f"Проверьте: {', '.join(duplicate_codes)}."
        )

    legacy_fraction = finite("omega_c22_overlap_legacy_fraction")
    dpa_rt = finite("C22_5_found_rt")
    c20_target_rt = finite("C20_3N8_corrected_target_rt")
    epa_rt = finite("C20_5_found_rt")
    c22_4_left_width = finite("C22_4_left_width")

    unstable_profile = bool(
        not custom_profile
        and np.all(np.isfinite([
                legacy_fraction,
                dpa_rt,
                c20_target_rt,
                epa_rt,
        ]))
        and legacy_fraction <= HIGH_ERROR_LEGACY_FRACTION_SPLIT
        and dpa_rt <= HIGH_ERROR_DPA_RT_MAX
        and c20_target_rt > HIGH_ERROR_C20_TARGET_RT_MIN
        and epa_rt <= HIGH_ERROR_EPA_RT_MAX
    )
    if unstable_profile:
        score = max(score, 88)
        reason_codes.append("c20_c22_sensitive_profile")
        peak_codes.update(["C20:5", "C22:5", "C22:4"])
        reasons.append(
            "Такое сочетание пиков C20 и C22 раньше давало ошибки больше 0,5. "
            "Сверьте площади и обе границы C20:5, C22:5 и C22:4."
        )

    early_retention_profile = bool(
        not custom_profile
        and np.isfinite(c20_target_rt)
        and c20_target_rt <= HIGH_ERROR_EARLY_C20_TARGET_RT_MAX
    )
    if early_retention_profile:
        score = max(score, 88)
        reason_codes.append("early_retention_profile")
        peak_codes.update(["C20:5", "C20:3N8", "C22:6", "C22:5", "C22:4"])
        reasons.append(
            "Пики вышли заметно раньше обычного, поэтому автоматические границы менее надёжны. "
            "Сверьте площади и границы пиков C20 и C22."
        )

    broad_shared_c22_tail = bool(
        not custom_profile
        and np.all(np.isfinite([legacy_fraction, c22_4_left_width]))
        and legacy_fraction > HIGH_ERROR_LEGACY_FRACTION_SPLIT
        and c22_4_left_width > HIGH_ERROR_C22_4_LEFT_WIDTH_MIN
    )
    if broad_shared_c22_tail:
        score = max(score, 88)
        reason_codes.append("broad_shared_c22_tail")
        peak_codes.update(["C22:5", "C22:4"])
        reasons.append(
            "Левая граница C22:4 захватывает широкий общий хвост. "
            "Сверьте площади и разделение C22:5/C22:4."
        )

    review_priority = "normal"
    if (
        score >= 85
        and judge_calibrated
        and np.isfinite(c22_height_ratio)
        and np.isfinite(c22_height_threshold)
        and c22_height_ratio >= c22_height_threshold
    ):
        review_priority = "c22_first"
        reason_codes.append("c22_height_priority")
        peak_codes.update(["C22:5", "C22:4"])
        reasons.append(
            "Сначала проверьте разделение C22:5 и C22:4: соотношение высот указывает на повышенный риск."
        )

    if score >= 95:
        band = "HIGH_RISK_STRUCTURAL"
    elif score >= 85:
        band = "HIGH_RISK_GT_0_5"
    else:
        band = "LOW_RISK"
    return {
        "score": int(score),
        "band": band,
        "reason_codes": reason_codes,
        "reasons": reasons,
        "peak_codes": sorted(peak_codes),
        "review_priority": review_priority,
        "c22_height_ratio": c22_height_ratio,
        "c22_height_priority_threshold": c22_height_threshold,
        "instrument_profile_judge_calibrated": judge_calibrated,
        "instrument_profile_judge_emergency_enabled": emergency_enabled,
        "instrument_profile_judge_target_abs_error": emergency_target,
    }


def assess_high_error_risk(matched_targets: pd.DataFrame, omega: Mapping[str, object]) -> dict:
    """Build production features and run the high-recall integration judge."""
    if matched_targets is None or matched_targets.empty:
        return classify_high_error_risk({"missing_key_peaks": ["C20:5", "C22:6", "C22:5"]})

    valid = matched_targets.copy()
    for column in [
        "area",
        "found_rt",
        "corrected_target_rt",
        "integration_start_x",
        "integration_end_x",
        "matched_peak_id",
        "peak_height_smooth",
    ]:
        valid[column] = (
            pd.to_numeric(valid[column], errors="coerce")
            if column in valid
            else np.nan
        )

    def target_value(code: str, column: str) -> float:
        if column not in valid:
            return np.nan
        row = valid.loc[valid["code"] == code, column]
        if row.empty or not np.isfinite(row.iloc[0]):
            return np.nan
        return float(row.iloc[0])

    key_codes = ["C20:5", "C22:6", "C22:5"]
    if "instrument_profile_omega_component" in valid:
        component_mask = valid["instrument_profile_omega_component"].fillna(False).astype(bool)
        selected_codes = valid.loc[component_mask, "code"].astype(str).tolist()
        if selected_codes:
            key_codes = selected_codes

    missing_key_peaks = []
    for code in key_codes:
        area = target_value(code, "area")
        found_rt = target_value(code, "found_rt")
        if not np.isfinite(area) or area <= 0 or not np.isfinite(found_rt):
            missing_key_peaks.append(code)

    matched_ids = valid.dropna(subset=["matched_peak_id"])[["code", "matched_peak_id"]]
    duplicate_mask = matched_ids["matched_peak_id"].duplicated(keep=False)
    duplicate_codes = matched_ids.loc[duplicate_mask, "code"].astype(str).tolist()

    c22_4_rt = target_value("C22:4", "found_rt")
    c22_4_start = target_value("C22:4", "integration_start_x")
    c22_5_height = target_value("C22:5", "peak_height_smooth")
    c22_4_height = target_value("C22:4", "peak_height_smooth")
    c22_height_ratio = (
        c22_5_height / c22_4_height
        if np.isfinite(c22_5_height) and np.isfinite(c22_4_height) and c22_4_height > 0
        else np.nan
    )
    custom_profile = rt_profile.uses_custom_instrument_profile(valid)
    judge_calibrated = not custom_profile
    if "instrument_profile_judge_calibrated" in valid:
        calibration_value = valid["instrument_profile_judge_calibrated"].iloc[0]
        if pd.notna(calibration_value):
            judge_calibrated = bool(calibration_value)
    emergency_enabled = False
    if "instrument_profile_judge_emergency_enabled" in valid:
        emergency_value = valid["instrument_profile_judge_emergency_enabled"].iloc[0]
        if pd.notna(emergency_value):
            emergency_enabled = bool(emergency_value)
    emergency_target = np.nan
    if "instrument_profile_judge_target_abs_error" in valid:
        emergency_target = pd.to_numeric(
            valid["instrument_profile_judge_target_abs_error"], errors="coerce"
        ).iloc[0]

    emergency_codes: set[str] = set()
    emergency_reasons: list[str] = []

    def emergency_issue(codes: list[str], reason: str) -> None:
        emergency_codes.update(codes)
        if reason not in emergency_reasons:
            emergency_reasons.append(reason)

    if custom_profile and emergency_enabled:
        for code in key_codes:
            found_rt = target_value(code, "found_rt")
            target_rt = target_value(code, "corrected_target_rt")
            start_x = target_value(code, "integration_start_x")
            end_x = target_value(code, "integration_end_x")
            if not np.all(np.isfinite([found_rt, target_rt, start_x, end_x])):
                emergency_issue([code], f"У {code} не удалось проверить положение и обе границы.")
                continue
            width = end_x - start_x
            if abs(found_rt - target_rt) > EMERGENCY_KEY_RT_RESIDUAL_MAX:
                emergency_issue([code], f"Вершина {code} слишком далеко от ожидаемого положения.")
            if not EMERGENCY_KEY_WIDTH_MIN <= width <= EMERGENCY_KEY_WIDTH_MAX:
                emergency_issue([code], f"Ширина {code} нетипична для расчётного пика.")
            if (
                found_rt - start_x < EMERGENCY_APEX_EDGE_MIN
                or end_x - found_rt < EMERGENCY_APEX_EDGE_MIN
            ):
                emergency_issue([code], f"Граница {code} проходит почти по вершине пика.")

            status_row = (
                valid.loc[valid["code"] == code, "status"]
                if "status" in valid
                else pd.Series(dtype=object)
            )
            status = (
                ""
                if status_row.empty or pd.isna(status_row.iloc[0])
                else str(status_row.iloc[0]).casefold()
            )
            if any(
                token in status
                for token in ("not_found", "unresolved", "identity_center_locked", "recovered_")
            ):
                emergency_issue([code], f"{code} был восстановлен неуверенно и требует просмотра.")

        c20_5_height = target_value("C20:5", "peak_height_smooth")
        c20_3_height = target_value("C20:3N8", "peak_height_smooth")
        c20_height_ratio = (
            c20_5_height / c20_3_height
            if np.isfinite(c20_5_height) and np.isfinite(c20_3_height) and c20_3_height > 0
            else np.nan
        )
        if np.isfinite(c20_height_ratio) and not (
            EMERGENCY_C20_HEIGHT_RATIO_MIN
            <= c20_height_ratio
            <= EMERGENCY_C20_HEIGHT_RATIO_MAX
        ):
            emergency_issue(
                ["C20:5", "C20:3N8"],
                "Соотношение высот C20 нетипично: сначала проверьте, что выбраны правильные вершины.",
            )

        if np.isfinite(c22_height_ratio) and not (
            EMERGENCY_C22_HEIGHT_RATIO_MIN
            <= c22_height_ratio
            <= EMERGENCY_C22_HEIGHT_RATIO_MAX
        ):
            emergency_issue(
                ["C22:5", "C22:4"],
                "Соотношение высот C22 нетипично: проверьте назначение и разделение пиков.",
            )
    height_threshold = LEGACY_C22_HEIGHT_PRIORITY_RATIO if not custom_profile else np.nan
    if "instrument_profile_c22_height_priority_threshold" in valid:
        profile_threshold = pd.to_numeric(
            valid["instrument_profile_c22_height_priority_threshold"], errors="coerce"
        ).iloc[0]
        if np.isfinite(profile_threshold):
            height_threshold = float(profile_threshold)
    c20_judge_target_rt = target_value("C20:3N8", "corrected_target_rt")
    if not rt_profile.uses_custom_instrument_profile(valid):
        # Legacy high-risk thresholds were calibrated on expected_rt + the
        # global matcher shift.  The GUI now displays the real resolved C20
        # cluster centres, but changing that display must not silently retune
        # the production judge or create extra warnings.
        c20_row = valid[valid["code"] == "C20:3N8"]
        if not c20_row.empty:
            expected_rt = pd.to_numeric(c20_row.get("expected_rt"), errors="coerce").iloc[0]
            local_shift = pd.to_numeric(c20_row.get("rt_local_shift"), errors="coerce").iloc[0]
            if np.isfinite(expected_rt) and np.isfinite(local_shift):
                c20_judge_target_rt = float(expected_rt + local_shift)
    features = {
        "missing_key_peaks": missing_key_peaks,
        "duplicate_peak_codes": duplicate_codes,
        "instrument_profile_custom_rt": custom_profile,
        "instrument_profile_judge_calibrated": judge_calibrated,
        "instrument_profile_judge_emergency_enabled": emergency_enabled,
        "instrument_profile_judge_target_abs_error": emergency_target,
        "emergency_geometry_codes": sorted(emergency_codes),
        "emergency_geometry_reasons": emergency_reasons,
        "C22_height_ratio": c22_height_ratio,
        "C22_height_priority_threshold": height_threshold,
        "omega_c22_overlap_legacy_fraction": omega.get("c22_overlap_legacy_fraction", np.nan),
        "C22_5_found_rt": target_value("C22:5", "found_rt"),
        "C20_3N8_corrected_target_rt": c20_judge_target_rt,
        "C20_5_found_rt": target_value("C20:5", "found_rt"),
        "C22_4_left_width": c22_4_rt - c22_4_start
        if np.isfinite(c22_4_rt) and np.isfinite(c22_4_start)
        else np.nan,
    }
    return classify_high_error_risk(features)


def assess_confidence(
    matched_targets: pd.DataFrame,
    peaks: pd.DataFrame,
    omega: dict,
    baseline_mode: str,
    cluster_quality_score: float,
) -> dict:
    result = {
        "score": np.nan,
        "level": "unknown",
        "label": "—",
        "button_text": "Уверенность: —",
        "reasons": [],
        "metrics": [],
    }
    if matched_targets is None or matched_targets.empty or not isinstance(omega, dict):
        return result

    valid = matched_targets.copy()
    valid["area"] = pd.to_numeric(valid.get("area"), errors="coerce")
    valid["found_rt"] = pd.to_numeric(valid.get("found_rt"), errors="coerce")
    valid["matched_peak_id"] = pd.to_numeric(valid.get("matched_peak_id"), errors="coerce")
    valid["integration_start_x"] = pd.to_numeric(valid.get("integration_start_x"), errors="coerce")
    valid["integration_end_x"] = pd.to_numeric(valid.get("integration_end_x"), errors="coerce")

    score = 100.0
    reason_items: list[tuple[float, str]] = []

    def penalize(points: float, reason: str) -> None:
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

    high_error_risk = assess_high_error_risk(valid, omega)
    risk_reasons = list(high_error_risk.get("reasons", []))
    if high_error_risk.get("score", 0) >= 95 and risk_reasons:
        reason_items.append((65.0, risk_reasons[0]))
        for reason in risk_reasons[1:]:
            reason_items.append((0.0, reason))
    elif high_error_risk.get("score", 0) >= 85 and risk_reasons:
        reason_items.append((46.0, risk_reasons[0]))
        for reason in risk_reasons[1:]:
            reason_items.append((0.0, reason))

    omega_value = float(omega.get("omega3_trio", np.nan))
    strict_value = float(omega.get("omega3_trio_strict", np.nan))
    spread = abs(omega_value - strict_value) if np.isfinite(omega_value) and np.isfinite(strict_value) else np.nan

    if baseline_mode != "chebyshev":
        penalize(4.0, "Фон сигнала оказался сложным — проверьте границы пиков")

    if cluster_quality_score is not None and np.isfinite(cluster_quality_score) and cluster_quality_score < CLUSTER_QUALITY_COMPLETE_SCORE:
        gap = float(CLUSTER_QUALITY_COMPLETE_SCORE - cluster_quality_score)
        penalize(min(18.0, 6.0 + gap * 0.45), "Не все группы пиков распознаны уверенно")

    matched_count = int(valid["matched_peak_id"].notna().sum())
    matched_ids = valid.dropna(subset=["matched_peak_id"])[["code", "matched_peak_id"]].copy()
    duplicate_ids = matched_ids[matched_ids["matched_peak_id"].duplicated(keep=False)]
    duplicate_codes = sorted(duplicate_ids["code"].astype(str).unique().tolist())
    if duplicate_codes:
        penalize(
            min(30.0, 18.0 + 3.0 * len(duplicate_codes)),
            "Один участок сигнала похож сразу на несколько пиков",
        )
    missing_count = int(max(0, len(valid) - matched_count))
    if missing_count > 0:
        penalize(min(18.0, 4.0 * missing_count), "Некоторые пики не удалось определить")

    trio_codes = ["C20:5", "C22:6", "C22:5"]
    trio_missing = [code for code in trio_codes if area_of(code) <= 0]
    if trio_missing:
        if bool(omega.get("c22_missing_dha_coelution_applied", False)) and trio_missing == ["C22:6"]:
            penalize(22.0, "DHA и DPA не разделились — этот участок нужно проверить вручную")
        else:
            penalize(35.0, f"Не найден важный пик: {', '.join(trio_missing)}")

    if np.isfinite(spread):
        if spread > 0.65:
            penalize(14.0, "Автоматическая поправка заметно изменила результат")
        elif spread > 0.45:
            penalize(7.0, "Результат потребовал дополнительной автоматической поправки")

    c18_scale = float(omega.get("c18_denominator_scale", 1.0))
    c18_1 = area_of("C18:1N9C")
    c18_2 = area_of("C18:2N6C")
    c18_3 = area_of("C18:3N3")
    c18_width = width_of("C18:1N9C")
    c18_ratio = c18_1 / c18_2 if c18_2 > 0 else np.nan
    c18_n3_fraction = c18_3 / c18_1 if c18_1 > 0 else np.nan
    c18_status = status_of("C18:1N9C")
    if c18_scale < 0.999:
        penalize(8.0 if c18_scale >= 0.80 else 12.0, "Группа C18 потребовала дополнительной поправки")
    if (
        np.isfinite(c18_ratio)
        and np.isfinite(c18_n3_fraction)
        and np.isfinite(c18_width)
        and c18_ratio > C18_DENOMINATOR_EXTREME_RATIO
        and c18_n3_fraction < C18_DENOMINATOR_EXTREME_SMALL_N3_FRACTION
        and c18_width > C18_DENOMINATOR_EXTREME_WIDTH_MIN
    ):
        penalize(10.0, "Один из пиков C18 выглядит слишком широким")
    if "matched_c18_local_bounds" in c18_status or "matched_c18_pvfit" in c18_status:
        penalize(6.0, "Границы пиков C18 были уточнены автоматически")

    epa = area_of("C20:5")
    c20_3 = area_of("C20:3N8")
    epa_status = status_of("C20:5")
    epa_credit = float(omega.get("epa_overlap_credit_area", 0.0))
    epa_ratio = epa / c20_3 if c20_3 > 0 else np.nan
    w_epa = width_of("C20:5")
    w_c20_3 = width_of("C20:3N8")
    if "identity_center_locked" in epa_status:
        penalize(10.0, "Математическая модель C20:5 ушла в сторону — проверьте выбранную вершину")
    elif "recovered_c20_missing_epa_component" in epa_status:
        penalize(10.0, "Пик C20:5 восстановлен из слившегося участка — проверьте его границы")
    elif "matched_c20_fit" in epa_status:
        penalize(10.0, "Пик C20:5 выделен неуверенно — проверьте его границы")
    elif "matched_c20_local" in epa_status:
        penalize(5.0, "Границы пика C20:5 были уточнены автоматически")
    if epa_credit > C20_EPA_UNDERFIT_CREDIT_MIN:
        penalize(5.0, "Пик C20:5 частично сливается с соседним пиком")
    if (
        np.isfinite(epa_ratio)
        and np.isfinite(w_epa)
        and np.isfinite(w_c20_3)
        and epa_ratio < C20_EPA_UNDERFIT_RATIO_MAX
        and w_c20_3 > w_epa * C20_EPA_UNDERFIT_WIDTH_RATIO
    ):
        penalize(10.0, "Пик C20:5 может быть выделен не полностью")
    elif np.isfinite(epa_ratio) and epa_ratio < 0.30 and "matched_c20_local" in epa_status:
        penalize(12.0, "Пик C20:5 слишком мал по сравнению с соседним")

    c22_statuses = [status_of(code) for code in ["C22:6", "C22:5", "C22:4"]]
    c22_status_text = " ".join(status for status in c22_statuses if status)
    c22_credit = float(omega.get("c22_overlap_credit_area", 0.0))
    c22_widths = np.asarray([width_of("C22:6"), width_of("C22:5"), width_of("C22:4")], dtype=float)
    c22_mean_width = float(np.nanmean(c22_widths)) if np.isfinite(c22_widths).any() else np.nan
    if "recovered_c22_local_unresolved" in c22_status_text:
        penalize(14.0, "Пик C22:5 восстановлен автоматически — проверьте его границы")
    elif "matched_c22_pvfit" in c22_status_text:
        penalize(12.0, "Пики C22 сливаются и требуют визуальной проверки")
    elif "matched_c22_fit" in c22_status_text:
        penalize(10.0, "Один из пиков C22 восстановлен автоматически")
    elif "tailtight" in c22_status_text:
        penalize(6.0, "Хвосты пиков C22 пришлось уточнить автоматически")
    if c22_credit > 0:
        penalize(2.0 if c22_credit < 80 else 4.0, "Пики C22 частично перекрываются")
    if np.isfinite(c22_mean_width):
        if c22_mean_width > 0.036:
            penalize(12.0, "Пики C22 сильно сливаются — проверьте границы")
        elif c22_mean_width > 0.032:
            penalize(6.0, "Пики C22 расположены близко друг к другу")

    peak_count = len(peaks) if peaks is not None else 0
    if peak_count >= 65:
        penalize(5.0, "На хроматограмме много лишних пиков")

    geometry_score = max(0.0, min(100.0, score))
    score = geometry_score
    if high_error_risk.get("score", 0) >= 95:
        score = min(score, 35.0)
    elif high_error_risk.get("score", 0) >= 85:
        score = min(score, 54.0)
    if score >= GEOMETRY_READY_SCORE:
        level = "Геометрия OK"
    elif score >= 70.0:
        level = "Быстрая проверка"
    elif score >= GEOMETRY_STOP_SCORE:
        level = "Проверить границы"
    else:
        level = "Ручная проверка"

    reason_items.sort(key=lambda item: item[0], reverse=True)
    reasons = [f"• {text}" for _points, text in reason_items]
    metrics = []
    if np.isfinite(strict_value) and np.isfinite(omega_value):
        metrics.append(f"Omega final / strict: {omega_value:.2f}% / {strict_value:.2f}%")
    if np.isfinite(spread):
        metrics.append(f"Разброс strict/final: {spread:.2f}")
    metrics.append(f"Baseline: {baseline_mode}")
    if cluster_quality_score is not None and np.isfinite(cluster_quality_score):
        metrics.append(f"Cluster quality: {cluster_quality_score:.1f}")
    metrics.append(f"Matched peaks: {matched_count}/{len(valid)}")
    metrics.append(f"Уникальные peak ID: {matched_ids['matched_peak_id'].nunique()}/{matched_count}")
    if np.isfinite(epa_ratio):
        metrics.append(f"C20 area EPA/C20:3 = {epa_ratio:.2f}")
    c22_ratio = float(omega.get("c22_reference_ratio", np.nan))
    if np.isfinite(c22_ratio):
        metrics.append(f"C22 area DPA/C22:4 = {c22_ratio:.2f}")
    c22_height_ratio = float(high_error_risk.get("c22_height_ratio", np.nan))
    if np.isfinite(c22_height_ratio):
        metrics.append(f"C22 height DPA/C22:4 = {c22_height_ratio:.2f}")
    c22_debit = float(omega.get("c22_overintegration_debit_points", 0.0))
    metrics.append(f"C22 коррекции: credit area {c22_credit:.1f}; debit {c22_debit:.2f} п.п.")

    result["score"] = score
    result["geometry_score"] = geometry_score
    result["level"] = level
    result["label"] = level
    if high_error_risk.get("score", 0) >= 95:
        result["button_text"] = "СТОП — переинтегрировать пики"
    elif high_error_risk.get("score", 0) >= 85:
        risky_codes = ", ".join(high_error_risk.get("peak_codes", []))
        result["button_text"] = f"ПРОВЕРИТЬ — {risky_codes or 'границы'}"
    elif geometry_score < GEOMETRY_STOP_SCORE:
        result["button_text"] = "СТОП — ручная проверка пиков"
    elif geometry_score < GEOMETRY_READY_SCORE:
        result["button_text"] = "ПРОВЕРИТЬ — геометрию пиков"
    else:
        result["button_text"] = "ГОТОВО — ручная правка не нужна"
    result["reasons"] = reasons
    result["metrics"] = metrics
    result["high_error_risk"] = high_error_risk
    return result


def annotate_result(result: dict, baseline_mode: str) -> dict:
    annotated = dict(result)
    annotated["baseline_mode"] = baseline_mode
    annotated["cluster_quality_score"] = compute_cluster_quality(annotated["matched_targets_df"])
    annotated["confidence"] = assess_confidence(
        matched_targets=annotated["matched_targets_df"],
        peaks=annotated["peaks_df"],
        omega=annotated["omega"],
        baseline_mode=annotated["baseline_mode"],
        cluster_quality_score=annotated["cluster_quality_score"],
    )
    return annotated
