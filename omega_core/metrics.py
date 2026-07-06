from __future__ import annotations

import math

import numpy as np
import pandas as pd


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

# Disabled by default after July-regression validation: this overlap-credit
# model was the dominant source of large C20/EPA over-estimation outliers.
# Keep the parameters below for future bounded/gated replacement experiments.
ENABLE_DATA_DRIVEN_C20_EPA_MODEL = False
C20_EPA_MODEL_GATE_RATIO_MAX = 0.70
C20_EPA_MODEL_BLEND = 0.25
C20_EPA_OVERLAP_WIDE_NEIGHBOR_RATIO = 1.30
C20_EPA_OVERLAP_EXTRA_SCALE = 1.60
C20_EPA_UNDERFIT_RATIO_MAX = 0.13
C20_EPA_UNDERFIT_WIDTH_RATIO = 1.80
C20_EPA_UNDERFIT_STRICT_MAX = 4.20
C20_EPA_UNDERFIT_CREDIT_MIN = 20.0
C20_EPA_UNDERFIT_EXTRA_SCALE = 1.50
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
        row = valid[valid["code"] == code]
        if row.empty or "status" not in row:
            return ""
        value = row["status"].iloc[0]
        return "" if pd.isna(value) else str(value)

    epa, dha, dpa = area_of("C20:5"), area_of("C22:6"), area_of("C22:5")
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
    epa_to_c20_3_ratio = epa / c20_3 if c20_3 > 0 else np.nan
    if (
        ENABLE_DATA_DRIVEN_C20_EPA_MODEL
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

    c22_ratio = dpa / c22_4 if c22_4 > 0 else np.nan
    c22_width_values = np.asarray([width_of("C22:6"), width_of("C22:5"), width_of("C22:4")], dtype=float)
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
    c22_credit_area = c22_4 * c22_fraction

    c22_debit_area = 0.0
    c22_debit_points = 0.0
    c22_debit_applied = False
    if (
        c22_4 > 0
        and dpa > 0
        and np.isfinite(c22_ratio)
        and c22_ratio > C22_DPA_OVERINTEGRATION_RATIO_MIN
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
    w_dha, w_dpa, w_c22_4 = c22_width_values
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

    omega_value = float(omega.get("omega3_trio", np.nan))
    strict_value = float(omega.get("omega3_trio_strict", np.nan))
    spread = abs(omega_value - strict_value) if np.isfinite(omega_value) and np.isfinite(strict_value) else np.nan

    if baseline_mode != "chebyshev":
        penalize(8.0, "Потребовался fallback baseline")

    if cluster_quality_score is not None and np.isfinite(cluster_quality_score) and cluster_quality_score < CLUSTER_QUALITY_COMPLETE_SCORE:
        gap = float(CLUSTER_QUALITY_COMPLETE_SCORE - cluster_quality_score)
        penalize(min(18.0, 6.0 + gap * 0.45), f"Качество кластеров ниже целевого ({cluster_quality_score:.1f})")

    matched_count = int(valid["matched_peak_id"].notna().sum())
    missing_count = int(max(0, len(valid) - matched_count))
    if missing_count > 0:
        penalize(min(18.0, 4.0 * missing_count), f"Есть неполные матчинг-пики ({missing_count})")

    trio_codes = ["C20:5", "C22:6", "C22:5"]
    trio_missing = [code for code in trio_codes if area_of(code) <= 0]
    if trio_missing:
        penalize(30.0, f"Отсутствуют ключевые omega-3 пики: {', '.join(trio_missing)}")

    if np.isfinite(spread):
        if spread > 0.45:
            penalize(18.0, f"Сильный разброс strict/final omega ({spread:.2f})")
        elif spread > 0.25:
            penalize(10.0, f"Заметный разброс strict/final omega ({spread:.2f})")
        elif spread > 0.12:
            penalize(5.0, f"Небольшой разброс strict/final omega ({spread:.2f})")

    c18_scale = float(omega.get("c18_denominator_scale", 1.0))
    c18_1 = area_of("C18:1N9C")
    c18_2 = area_of("C18:2N6C")
    c18_3 = area_of("C18:3N3")
    c18_width = width_of("C18:1N9C")
    c18_ratio = c18_1 / c18_2 if c18_2 > 0 else np.nan
    c18_n3_fraction = c18_3 / c18_1 if c18_1 > 0 else np.nan
    c18_status = status_of("C18:1N9C")
    if c18_scale < 0.999:
        penalize(8.0 if c18_scale >= 0.80 else 12.0, f"Сработала коррекция denominator для C18 ({c18_scale:.2f})")
    if (
        np.isfinite(c18_ratio)
        and np.isfinite(c18_n3_fraction)
        and np.isfinite(c18_width)
        and c18_ratio > C18_DENOMINATOR_EXTREME_RATIO
        and c18_n3_fraction < C18_DENOMINATOR_EXTREME_SMALL_N3_FRACTION
        and c18_width > C18_DENOMINATOR_EXTREME_WIDTH_MIN
    ):
        penalize(10.0, "C18:1N9C выглядит перерастянутым относительно C18-кластера")
    if "matched_c18_local_bounds" in c18_status or "matched_c18_pvfit" in c18_status:
        penalize(6.0, "C18-кластер потребовал локальную переинтеграцию")

    epa = area_of("C20:5")
    c20_3 = area_of("C20:3N8")
    epa_status = status_of("C20:5")
    epa_credit = float(omega.get("epa_overlap_credit_area", 0.0))
    epa_ratio = epa / c20_3 if c20_3 > 0 else np.nan
    w_epa = width_of("C20:5")
    w_c20_3 = width_of("C20:3N8")
    if "matched_c20_fit" in epa_status:
        penalize(10.0, "EPA восстанавливался через C20 fit")
    elif "matched_c20_local" in epa_status:
        penalize(5.0, "EPA потребовал локальную C20-коррекцию")
    if epa_credit > C20_EPA_UNDERFIT_CREDIT_MIN:
        penalize(5.0, f"Существенный overlap-credit у EPA ({epa_credit:.1f})")
    if (
        np.isfinite(epa_ratio)
        and np.isfinite(w_epa)
        and np.isfinite(w_c20_3)
        and epa_ratio < C20_EPA_UNDERFIT_RATIO_MAX
        and w_c20_3 > w_epa * C20_EPA_UNDERFIT_WIDTH_RATIO
    ):
        penalize(10.0, "EPA выглядит недобранным относительно C20:3N8")

    c22_statuses = [status_of(code) for code in ["C22:6", "C22:5", "C22:4"]]
    c22_status_text = " ".join(status for status in c22_statuses if status)
    c22_credit = float(omega.get("c22_overlap_credit_area", 0.0))
    c22_widths = np.asarray([width_of("C22:6"), width_of("C22:5"), width_of("C22:4")], dtype=float)
    c22_mean_width = float(np.nanmean(c22_widths)) if np.isfinite(c22_widths).any() else np.nan
    if "matched_c22_pvfit" in c22_status_text:
        penalize(12.0, "C22-кластер потребовал pvfit refinement")
    elif "matched_c22_fit" in c22_status_text:
        penalize(10.0, "C22-кластер потребовал fit-восстановление")
    elif "tailtight" in c22_status_text:
        penalize(6.0, "C22-кластер потребовал tail tightening")
    if c22_credit > 0:
        penalize(4.0 if c22_credit < 80 else 8.0, f"C22 overlap-credit участвует в расчёте ({c22_credit:.1f})")
    if np.isfinite(c22_mean_width):
        if c22_mean_width > 0.036:
            penalize(12.0, f"C22-кластер всё ещё широкий ({c22_mean_width:.3f} min)")
        elif c22_mean_width > 0.032:
            penalize(6.0, f"C22-кластер умеренно широкий ({c22_mean_width:.3f} min)")

    peak_count = len(peaks) if peaks is not None else 0
    if peak_count >= 65:
        penalize(5.0, f"Необычно много детектированных пиков ({peak_count})")

    score = max(0.0, min(100.0, score))
    if score >= 85.0:
        level = "Авто OK"
    elif score >= 70.0:
        level = "Быстро проверить"
    elif score >= 55.0:
        level = "Проверить руками"
    else:
        level = "Ручная проверка"

    reason_items.sort(key=lambda item: item[0], reverse=True)
    reasons = [f"-{int(round(points))}: {text}" for points, text in reason_items]
    metrics = []
    if np.isfinite(strict_value) and np.isfinite(omega_value):
        metrics.append(f"Omega final / strict: {omega_value:.2f}% / {strict_value:.2f}%")
    if np.isfinite(spread):
        metrics.append(f"Разброс strict/final: {spread:.2f}")
    metrics.append(f"Baseline: {baseline_mode}")
    if cluster_quality_score is not None and np.isfinite(cluster_quality_score):
        metrics.append(f"Cluster quality: {cluster_quality_score:.1f}")
    metrics.append(f"Matched peaks: {matched_count}/{len(valid)}")

    result["score"] = score
    result["level"] = level
    result["label"] = level
    result["button_text"] = f"Уверенность: {int(round(score))}"
    result["reasons"] = reasons
    result["metrics"] = metrics
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
