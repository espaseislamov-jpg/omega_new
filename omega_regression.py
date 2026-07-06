from __future__ import annotations

import argparse
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

import omega_core
from omega_core import metrics, pipeline, signal


PROJECT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = Path(r"C:\Users\marat\Desktop\CSV_Omega")
DEFAULT_REFERENCE_PATH = PROJECT_DIR / "reference_targets_reverted_c22fixed.json"
OLD45_DATES = {"13032026", "14012026", "20032026"}
SAMPLE_NAME_RE = re.compile(r"^O(?P<instrument_no>\d+)_(?P<sample_id>\d+)\.D$", re.IGNORECASE)

# Known special case from the 03072026 instrument export: O1 is a blank/empty
# sample, so a position-only workbook would need +1 indexing. ID-based matching
# is still preferred whenever sample IDs are available.
POSITION_INDEX_OFFSETS = {"03072026": 1}
POSITION_MATCH_DATES = {"03072026"}
OMEGA_CODES = ("C20:5", "C22:5", "C22:6", "C20:3N8", "C22:4", "C18:1N9C", "C18:2N6C", "C18:3N3")
DIAGNOSTIC_TARGET_CODES = ("C20:5", "C20:3N8", "C20:4N6", "C22:6", "C22:5", "C22:4", "C18:1N9C", "C18:2N6C", "C18:3N3")
OMEGA_TARGET_CODES = ("C20:5", "C22:6", "C22:5")


@dataclass(frozen=True)
class ExcelReference:
    row_number: int
    ordinal: int
    token: int
    reference: float


@dataclass(frozen=True)
class BatchRecord:
    index: int
    instrument_no: int | None
    sample_id: int | None
    sample_name: str
    dataframe: pd.DataFrame


def normalize_date_token(date: str) -> str:
    text = str(date).strip()
    if len(text) == 6 and text.isdigit():
        return f"{text[:4]}20{text[4:]}"
    return text


def parse_sample_name(sample_name: str) -> tuple[int | None, int | None]:
    match = SAMPLE_NAME_RE.match(str(sample_name).strip())
    if not match:
        return None, None
    return int(match.group("instrument_no")), int(match.group("sample_id"))


def load_excel_refs(xlsx_path: Path) -> list[ExcelReference]:
    raw = pd.read_excel(xlsx_path, header=None)
    refs: list[ExcelReference] = []
    for excel_idx, row in enumerate(raw.itertuples(index=False), start=1):
        sample_token = None
        ref_value = None
        for value in row:
            if pd.isna(value):
                continue
            if sample_token is None:
                text = str(value).strip()
                if text.isdigit():
                    sample_token = int(text)
                    continue
            if ref_value is None:
                try:
                    ref_value = float(str(value).strip().replace(",", "."))
                except ValueError:
                    pass
        if sample_token is not None and ref_value is not None:
            refs.append(ExcelReference(excel_idx, len(refs) + 1, sample_token, ref_value))
    return refs


def discover_reference_pairs(data_dir: Path) -> list[tuple[Path, Path, str, str]]:
    pairs = []
    for xlsx_path in sorted(data_dir.rglob("test_bigbatch_*.xlsx")):
        raw_date = xlsx_path.stem.split("_")[-1]
        date = normalize_date_token(raw_date)
        csv_candidates = [
            xlsx_path.parent / f"{raw_date}.CSV",
            xlsx_path.parent / f"{date}.CSV",
            xlsx_path.parent / f"{raw_date}.csv",
            xlsx_path.parent / f"{date}.csv",
        ]
        csv_path = next((candidate for candidate in csv_candidates if candidate.exists()), None)
        if csv_path is not None:
            pairs.append((xlsx_path, csv_path, raw_date, date))
    return pairs


def build_batch_records(csv_path: Path) -> list[BatchRecord]:
    batches = omega_core.load_batches(csv_path, cutoff_minutes=4.0)
    records: list[BatchRecord] = []
    for index, batch in enumerate(batches):
        sample_name = str(batch.get("sample_name", ""))
        instrument_no, sample_id = parse_sample_name(sample_name)
        records.append(BatchRecord(index, instrument_no, sample_id, sample_name, batch["dataframe"]))
    return records


def resolve_batch(ref: ExcelReference, batches: list[BatchRecord], date: str) -> tuple[BatchRecord | None, str]:
    if date in POSITION_MATCH_DATES:
        batch_index = ref.ordinal - 1
        if 0 <= batch_index < len(batches):
            return batches[batch_index], "position_date_override"
        return None, "missing_position_date_override"

    by_sample_id = {batch.sample_id: batch for batch in batches if batch.sample_id is not None}
    by_instrument_no = {batch.instrument_no: batch for batch in batches if batch.instrument_no is not None}

    if ref.token > 1000:
        matched = by_sample_id.get(ref.token)
        if matched is not None:
            return matched, "sample_id"
        return None, "missing_sample_id"

    matched = by_instrument_no.get(ref.token)
    if matched is not None:
        return matched, "instrument_no"

    batch_index = ref.token - 1 + POSITION_INDEX_OFFSETS.get(date, 0)
    if 0 <= batch_index < len(batches):
        return batches[batch_index], "position"
    return None, "missing_position"


def _safe_float(value, default=np.nan) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float(default)
    return out if np.isfinite(out) else float(default)


def _target_row(matched_targets: pd.DataFrame, code: str) -> pd.Series | None:
    if matched_targets is None or matched_targets.empty or "code" not in matched_targets:
        return None
    row = matched_targets[matched_targets["code"] == code]
    if row.empty:
        return None
    return row.iloc[0]


def _target_width(matched_targets: pd.DataFrame, code: str) -> float:
    row = _target_row(matched_targets, code)
    if row is None:
        return np.nan
    start_x = _safe_float(row.get("integration_start_x"))
    end_x = _safe_float(row.get("integration_end_x"))
    if not (np.isfinite(start_x) and np.isfinite(end_x)):
        return np.nan
    return float(end_x - start_x)


def _target_status(matched_targets: pd.DataFrame, code: str) -> str:
    row = _target_row(matched_targets, code)
    if row is None:
        return ""
    value = row.get("status", "")
    return "" if pd.isna(value) else str(value)


def _extract_result_features(result: dict) -> dict:
    omega = result.get("omega", {}) if isinstance(result.get("omega"), dict) else {}
    confidence = result.get("confidence", {}) if isinstance(result.get("confidence"), dict) else {}
    matched = result.get("matched_targets_df")
    features = {
        "calculated": _safe_float(result.get("omega_report")),
        "confidence": _safe_float(confidence.get("score")),
        "confidence_level": confidence.get("level", "") if isinstance(confidence, dict) else "",
        "baseline_mode": result.get("baseline_mode", ""),
        "cluster_quality_score": _safe_float(result.get("cluster_quality_score")),
        "rt_shift": _safe_float(result.get("rt_shift")),
        "best_window": result.get("best_window", np.nan),
    }
    for key, value in omega.items():
        if isinstance(value, (bool, np.bool_)):
            features[f"omega_{key}"] = bool(value)
        elif isinstance(value, (int, float, np.integer, np.floating)):
            features[f"omega_{key}"] = _safe_float(value)
    if isinstance(matched, pd.DataFrame):
        for code in DIAGNOSTIC_TARGET_CODES:
            slug = code.replace(":", "_").replace("/", "_")
            row = _target_row(matched, code)
            found_rt = _safe_float(row.get("found_rt")) if row is not None else np.nan
            expected_rt = _safe_float(row.get("expected_rt")) if row is not None else np.nan
            corrected_target_rt = _safe_float(row.get("corrected_target_rt")) if row is not None else np.nan
            start_x = _safe_float(row.get("integration_start_x")) if row is not None else np.nan
            end_x = _safe_float(row.get("integration_end_x")) if row is not None else np.nan
            width = _target_width(matched, code)
            features[f"{slug}_area"] = _safe_float(row.get("area")) if row is not None else np.nan
            features[f"{slug}_found_rt"] = found_rt
            features[f"{slug}_expected_rt"] = expected_rt
            features[f"{slug}_corrected_target_rt"] = corrected_target_rt
            features[f"{slug}_rt_error"] = found_rt - corrected_target_rt if np.isfinite(found_rt) and np.isfinite(corrected_target_rt) else np.nan
            features[f"{slug}_width"] = width
            features[f"{slug}_left_width"] = found_rt - start_x if np.isfinite(found_rt) and np.isfinite(start_x) else np.nan
            features[f"{slug}_right_width"] = end_x - found_rt if np.isfinite(end_x) and np.isfinite(found_rt) else np.nan
            if np.isfinite(width) and width > 0 and np.isfinite(found_rt) and np.isfinite(start_x) and np.isfinite(end_x):
                features[f"{slug}_asymmetry"] = max(found_rt - start_x, end_x - found_rt) / max(min(found_rt - start_x, end_x - found_rt), 1e-9)
            else:
                features[f"{slug}_asymmetry"] = np.nan
            features[f"{slug}_status"] = _target_status(matched, code)
    return features


def _judge_features(result: dict) -> dict:
    judge_decisions = result.get("judge_decisions_df")
    if isinstance(judge_decisions, pd.DataFrame) and not judge_decisions.empty:
        return {
            "judge_accepted": int((judge_decisions.get("decision") == "accepted").sum()),
            "judge_rejected": int((judge_decisions.get("decision") == "rejected").sum()),
            "judge_codes": ",".join(sorted({str(code) for code in judge_decisions.get("code", pd.Series(dtype=str)).dropna()})),
            "judge_reasons": ",".join(sorted({str(reason) for reason in judge_decisions.get("reason", pd.Series(dtype=str)).dropna()})),
        }
    return {"judge_accepted": 0, "judge_rejected": 0, "judge_codes": "", "judge_reasons": ""}


def classify_result(row: dict) -> tuple[str, str, str]:
    """Return a production-style review flag and regression outlier class.

    The flag intentionally avoids using the manual reference. The outlier class
    is only a diagnostic label for regression reports and may include the delta
    direction when a reference is available.
    """
    reasons: list[str] = []
    confidence = _safe_float(row.get("confidence"))
    cluster_quality = _safe_float(row.get("cluster_quality_score"))
    rejected = int(row.get("judge_rejected", 0) or 0)
    strict = _safe_float(row.get("omega_omega3_trio_strict"))
    corrected = _safe_float(row.get("omega_omega3_trio"))
    spread = abs(corrected - strict) if np.isfinite(corrected) and np.isfinite(strict) else np.nan

    if np.isfinite(confidence) and confidence < 40:
        reasons.append("low_confidence")
    elif np.isfinite(confidence) and confidence < 60:
        reasons.append("medium_confidence")
    if np.isfinite(cluster_quality) and cluster_quality < 50:
        reasons.append("low_cluster_quality")
    if rejected > 0:
        reasons.append("judge_rejected")
    if np.isfinite(spread) and spread > 0.75:
        reasons.append("strict_corrected_spread")

    status_text = " ".join(str(row.get(f"{code.replace(':', '_')}_status", "")) for code in ["C20:5", "C22:6", "C22:5", "C22:4"])
    if "not_found" in status_text:
        reasons.append("missing_target")
    if "tail" in status_text or "overlap" in status_text or "fit" in status_text:
        reasons.append("cluster_overlap_or_fit")

    flag = "OK"
    if any(reason in reasons for reason in ["low_confidence", "missing_target"]) or len(reasons) >= 3:
        flag = "REJECT"
    elif reasons:
        flag = "REVIEW"

    delta = _safe_float(row.get("delta"))
    direction = "unknown_delta"
    if np.isfinite(delta):
        direction = "overestimated" if delta > 0 else "underestimated"
    if "cluster_overlap_or_fit" in reasons:
        outlier_class = f"{direction}_cluster"
    elif "low_confidence" in reasons or "medium_confidence" in reasons:
        outlier_class = f"{direction}_low_confidence"
    elif "strict_corrected_spread" in reasons:
        outlier_class = f"{direction}_correction_spread"
    else:
        outlier_class = f"{direction}_unclassified"

    return flag, outlier_class, ",".join(reasons)


def _annotated_from_processed(processed: pd.DataFrame, reference_targets: pd.DataFrame, baseline_mode: str) -> dict:
    return metrics.annotate_result(pipeline.process_from_baseline(processed, reference_targets), baseline_mode=baseline_mode)


def evaluate_variants(dataframe: pd.DataFrame, reference_targets: pd.DataFrame, mode: str = "current") -> list[dict]:
    variants: list[dict] = []
    current = omega_core.process_batch(dataframe, reference_targets)
    current["variant_name"] = "current_pipeline"
    variants.append(current)
    if mode == "current":
        return variants

    candidate_builders = [
        ("chebyshev_q05", lambda df: signal.add_baseline(df, **{**signal.BASELINE_KWARGS, "lower_quantile": 0.05})),
        ("chebyshev_q12", lambda df: signal.add_baseline(df, **{**signal.BASELINE_KWARGS, "lower_quantile": 0.12})),
        ("asls_direct", signal.add_asls_baseline),
        ("arpls_direct", signal.add_arpls_baseline),
    ]
    for name, builder in candidate_builders:
        try:
            result = _annotated_from_processed(builder(dataframe), reference_targets, baseline_mode=name)
            result["variant_name"] = name
            variants.append(result)
        except Exception as exc:
            variants.append({"variant_name": name, "variant_error": f"{type(exc).__name__}: {exc}"})
    return variants


def select_variant(variants: list[dict], reference: float, selector_mode: str) -> dict:
    valid = [variant for variant in variants if np.isfinite(_safe_float(variant.get("omega_report")))]
    if not valid:
        error_text = "; ".join(filter(None, [variant.get("variant_error", "") for variant in variants]))
        raise ValueError(error_text or "No valid processing variants")
    if selector_mode == "oracle":
        # Research-only upper bound: shows how much error can be reduced by
        # choosing among already available processing variants. Do not use for
        # production patient results because it uses the manual reference.
        return min(valid, key=lambda variant: abs(_safe_float(variant.get("omega_report")) - float(reference)))
    return valid[0]


def build_audit_row(date: str, xlsx_path: Path, csv_path: Path, refs: list[ExcelReference], batches: list[BatchRecord]) -> dict:
    ref_tokens = {ref.token for ref in refs}
    batch_ids = {batch.sample_id for batch in batches if batch.sample_id is not None}
    matched_by_id = len({token for token in ref_tokens if token > 1000 and token in batch_ids})
    return {
        "date": date,
        "xlsx_path": str(xlsx_path),
        "csv_path": str(csv_path),
        "reference_rows": len(refs),
        "instrument_batches": len(batches),
        "matched_by_sample_id": matched_by_id,
        "reference_id_rows": sum(1 for token in ref_tokens if token > 1000),
        "position_rows": sum(1 for token in ref_tokens if token <= 1000),
        "missing_reference_ids": max(0, sum(1 for token in ref_tokens if token > 1000 and token not in batch_ids)),
    }


def dump_debug(debug_dir: Path, row: dict, result: dict) -> None:
    sample_dir = debug_dir / str(row["date"]) / str(row["sample_name"] or row["sample_no"])
    sample_dir.mkdir(parents=True, exist_ok=True)
    for key, filename in [
        ("processed_df", "processed_chromatogram.csv"),
        ("peaks_df", "detected_peaks.csv"),
        ("matched_targets_df", "matched_targets.csv"),
        ("judge_decisions_df", "judge_decisions.csv"),
    ]:
        value = result.get(key)
        if isinstance(value, pd.DataFrame):
            value.to_csv(sample_dir / filename, index=False)
    omega = result.get("omega", {}) if isinstance(result.get("omega"), dict) else {}
    metadata = {k: v for k, v in row.items() if isinstance(v, (str, int, float, bool)) or pd.isna(v)}
    (sample_dir / "omega_debug.json").write_text(
        json.dumps({"metadata": metadata, "omega": omega}, ensure_ascii=False, indent=2, default=str),
        encoding="utf-8",
    )
    _write_debug_plot(sample_dir / "plot.png", row, result)


def _write_debug_plot(path: Path, row: dict, result: dict) -> None:
    processed = result.get("processed_df")
    matched = result.get("matched_targets_df")
    if not isinstance(processed, pd.DataFrame) or processed.empty:
        return
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return

    x_col = "x_corrected" if "x_corrected" in processed.columns else "x"
    x = pd.to_numeric(processed[x_col], errors="coerce")
    y = pd.to_numeric(processed.get("y_corrected", processed.get("y")), errors="coerce")
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.plot(x, y, color="#1f77b4", linewidth=0.8, label="corrected signal")
    if "baseline" in processed and "y" in processed:
        ax.plot(x, pd.to_numeric(processed["y"], errors="coerce"), color="#cccccc", linewidth=0.5, alpha=0.8, label="raw signal")
        ax.plot(x, pd.to_numeric(processed["baseline"], errors="coerce"), color="#ff7f0e", linewidth=0.8, alpha=0.8, label="baseline")

    if isinstance(matched, pd.DataFrame) and not matched.empty:
        for _, target in matched.iterrows():
            code = str(target.get("code", ""))
            found_rt = _safe_float(target.get("found_rt"))
            start_x = _safe_float(target.get("integration_start_x"))
            end_x = _safe_float(target.get("integration_end_x"))
            if np.isfinite(start_x) and np.isfinite(end_x):
                ax.axvspan(start_x, end_x, color="#2ca02c", alpha=0.08)
            if np.isfinite(found_rt):
                ax.axvline(found_rt, color="#d62728", linewidth=0.5, alpha=0.65)
                ax.text(found_rt, ax.get_ylim()[1], code, rotation=90, fontsize=7, va="top", ha="center")

    title = (
        f"{row.get('date')} {row.get('sample_name')} "
        f"ref={_safe_float(row.get('reference')):.3f} calc={_safe_float(row.get('calculated')):.3f} "
        f"delta={_safe_float(row.get('delta')):.3f} flag={row.get('review_flag', '')}"
    )
    ax.set_title(title)
    ax.set_xlabel("Retention time, min")
    ax.set_ylabel("Signal")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run_current_engine(
    data_dir: Path = DEFAULT_DATA_DIR,
    reference_path: Path = DEFAULT_REFERENCE_PATH,
    selector_mode: str = "current",
    debug_dir: Path | None = None,
    debug_threshold: float = 0.5,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    reference_targets = omega_core.load_reference_targets(reference_path)
    rows: list[dict] = []
    audit_rows: list[dict] = []
    error_rows: list[dict] = []
    variant_rows: list[dict] = []

    for xlsx_path, csv_path, raw_date, date in discover_reference_pairs(data_dir):
        refs = load_excel_refs(xlsx_path)
        batches = build_batch_records(csv_path)
        audit_rows.append(build_audit_row(date, xlsx_path, csv_path, refs, batches))
        for ref in refs:
            batch, match_method = resolve_batch(ref, batches, date)
            if batch is None:
                error_rows.append({
                    "date": date,
                    "xlsx_path": str(xlsx_path),
                    "csv_path": str(csv_path),
                    "excel_row": ref.row_number,
                    "sample_no": ref.token,
                    "reference": ref.reference,
                    "match_method": match_method,
                    "error": "No matching instrument batch",
                })
                continue

            base_row = {
                "date": date,
                "raw_date": raw_date,
                "excel_row": ref.row_number,
                "sample_no": ref.token,
                "instrument_batch_index": batch.index + 1,
                "instrument_no": batch.instrument_no,
                "sample_id": batch.sample_id,
                "sample_name": batch.sample_name,
                "match_method": match_method,
                "reference": float(ref.reference),
            }
            try:
                variants = evaluate_variants(batch.dataframe, reference_targets, mode="variants" if selector_mode == "oracle" else "current")
                selected = select_variant(variants, ref.reference, selector_mode=selector_mode)
            except Exception as exc:
                error_row = {**base_row, "error": f"{type(exc).__name__}: {exc}"}
                error_rows.append(error_row)
                rows.append({**base_row, "calculated": np.nan, "delta": np.nan, "abs_delta": np.nan, "error": error_row["error"]})
                continue

            for variant in variants:
                features = _extract_result_features(variant) if "variant_error" not in variant else {}
                variant_rows.append({
                    **base_row,
                    "variant_name": variant.get("variant_name", ""),
                    "variant_error": variant.get("variant_error", ""),
                    **features,
                    "delta": _safe_float(features.get("calculated")) - float(ref.reference) if features else np.nan,
                })

            features = _extract_result_features(selected)
            calculated = features["calculated"]
            delta = calculated - float(ref.reference)
            out_row = {
                **base_row,
                **features,
                **_judge_features(selected),
                "selected_variant": selected.get("variant_name", "current_pipeline"),
                "delta": delta,
                "abs_delta": abs(delta),
                "error": "",
            }
            review_flag, outlier_class, review_reasons = classify_result(out_row)
            out_row["review_flag"] = review_flag
            out_row["outlier_class"] = outlier_class
            out_row["review_reasons"] = review_reasons
            rows.append(out_row)
            if debug_dir is not None and np.isfinite(out_row["abs_delta"]) and out_row["abs_delta"] > debug_threshold:
                dump_debug(debug_dir, out_row, selected)

    return pd.DataFrame(rows), pd.DataFrame(audit_rows), pd.DataFrame(error_rows), pd.DataFrame(variant_rows)


def summarize(results: pd.DataFrame) -> dict[str, float | int]:
    if results.empty or "abs_delta" not in results:
        return {"n": 0}
    valid = results.dropna(subset=["abs_delta"])
    if valid.empty:
        return {"n": 0}
    abs_delta = valid["abs_delta"]
    delta = valid["delta"]
    return {
        "n": int(len(valid)),
        "MAE": float(abs_delta.mean()),
        "RMSE": float(math.sqrt(float(np.mean(np.square(delta))))),
        "mean_delta": float(delta.mean()),
        "median_abs": float(abs_delta.median()),
        "std_delta": float(delta.std(ddof=0)),
        "within_0_2": int((abs_delta <= 0.2).sum()),
        "within_0_3": int((abs_delta <= 0.3).sum()),
        "within_0_4": int((abs_delta <= 0.4).sum()),
        "within_0_5": int((abs_delta <= 0.5).sum()),
        "within_0_6": int((abs_delta <= 0.6).sum()),
        "max_abs": float(abs_delta.max()),
    }


def build_summary_table(results: pd.DataFrame) -> pd.DataFrame:
    rows = [{"scope": "ALL", **summarize(results)}]
    old45 = results[results["date"].isin(OLD45_DATES)] if "date" in results else pd.DataFrame()
    rows.append({"scope": "OLD45", **summarize(old45)})
    if "date" in results:
        for date, group in results.groupby("date", sort=True):
            rows.append({"scope": date, **summarize(group)})
    return pd.DataFrame(rows)



def _issue_parts(row: pd.Series) -> list[str]:
    parts: list[str] = []
    confidence = _safe_float(row.get("confidence"))
    delta = _safe_float(row.get("delta"))
    c22_ratio = _safe_float(row.get("omega_c22_reference_ratio"))
    c22_debit_points = _safe_float(row.get("omega_c22_overintegration_debit_points"), 0.0)
    c22_credit_points = 100.0 * _safe_float(row.get("omega_c22_overlap_credit_area"), 0.0) / max(_safe_float(row.get("omega_effective_total_area"), 1.0), 1e-9)
    c20_status = " ".join(str(row.get(f"{code.replace(':', '_')}_status", "")) for code in ["C20:5", "C20:3N8", "C20:4N6"])
    c22_status = " ".join(str(row.get(f"{code.replace(':', '_')}_status", "")) for code in ["C22:6", "C22:5", "C22:4"])
    c18_status = " ".join(str(row.get(f"{code.replace(':', '_')}_status", "")) for code in ["C18:1N9C", "C18:2N6C", "C18:3N3"])

    if np.isfinite(confidence) and confidence < 60:
        parts.append("low_or_medium_confidence")
    if str(row.get("baseline_mode", "")) not in {"", "chebyshev"}:
        parts.append("baseline_fallback")
    if "not_found" in f"{c20_status} {c22_status} {c18_status}":
        parts.append("missing_target")
    if any(token in c22_status for token in ["tail", "baseexpand", "fit"]):
        parts.append("c22_complex_boundaries")
    if any(token in c20_status for token in ["baseexpand", "fit", "local"]):
        parts.append("c20_complex_boundaries")
    if any(token in c18_status for token in ["baseexpand", "valley"]):
        parts.append("c18_complex_boundaries")
    if np.isfinite(c22_ratio) and c22_ratio > 1.35:
        parts.append("high_dpa_to_c22_4_ratio")
    if c22_debit_points > 0.05:
        parts.append("c22_debit_applied")
    if c22_credit_points > 0.05:
        parts.append("c22_credit_applied")

    max_rt_error = 0.0
    max_width = 0.0
    max_asymmetry = 0.0
    for code in DIAGNOSTIC_TARGET_CODES:
        slug = code.replace(":", "_")
        rt_error = abs(_safe_float(row.get(f"{slug}_rt_error")))
        width = _safe_float(row.get(f"{slug}_width"))
        asymmetry = _safe_float(row.get(f"{slug}_asymmetry"))
        if np.isfinite(rt_error):
            max_rt_error = max(max_rt_error, rt_error)
        if np.isfinite(width):
            max_width = max(max_width, width)
        if np.isfinite(asymmetry):
            max_asymmetry = max(max_asymmetry, asymmetry)
    if max_rt_error > 0.035:
        parts.append("large_rt_error")
    if max_width > 0.075:
        parts.append("wide_peak_window")
    if max_asymmetry > 4.0:
        parts.append("asymmetric_peak_window")

    if np.isfinite(delta) and abs(delta) <= 0.3:
        parts.append("within_inter_operator_band")
    elif np.isfinite(delta) and abs(delta) <= 0.5:
        parts.append("within_clinical_band")
    elif np.isfinite(delta) and delta > 0:
        parts.append("overestimated_gt_0_5")
    elif np.isfinite(delta):
        parts.append("underestimated_gt_0_5")
    return parts


def classify_diagnostic_bucket(row: pd.Series) -> str:
    parts = set(_issue_parts(row))
    delta = _safe_float(row.get("delta"))
    if np.isfinite(delta) and abs(delta) <= 0.3:
        return "ok_within_0_3"
    if np.isfinite(delta) and abs(delta) <= 0.5:
        return "watch_within_0_5"
    direction = "over" if np.isfinite(delta) and delta > 0 else "under"
    if "c22_complex_boundaries" in parts or "high_dpa_to_c22_4_ratio" in parts or "c22_debit_applied" in parts:
        return f"{direction}_c22_cluster"
    if "c20_complex_boundaries" in parts:
        return f"{direction}_c20_cluster"
    if "c18_complex_boundaries" in parts or "baseline_fallback" in parts:
        return f"{direction}_baseline_or_c18_denominator"
    if "low_or_medium_confidence" in parts or "large_rt_error" in parts or "missing_target" in parts:
        return f"{direction}_low_confidence_or_rt"
    return f"{direction}_unclassified"


def build_diagnostic_tables(results: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if results.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    sample_rows = []
    peak_rows = []
    for _, row in results.iterrows():
        issue_parts = _issue_parts(row)
        bucket = classify_diagnostic_bucket(row)
        sample_key = {
            "date": row.get("date"),
            "sample_no": row.get("sample_no"),
            "sample_name": row.get("sample_name"),
            "reference": row.get("reference"),
            "calculated": row.get("calculated"),
            "delta": row.get("delta"),
            "abs_delta": row.get("abs_delta"),
            "confidence": row.get("confidence"),
            "review_flag": row.get("review_flag"),
            "outlier_class": row.get("outlier_class"),
            "diagnostic_bucket": bucket,
            "diagnostic_reasons": ",".join(issue_parts),
            "baseline_mode": row.get("baseline_mode"),
            "rt_shift": row.get("rt_shift"),
            "omega_c22_reference_ratio": row.get("omega_c22_reference_ratio"),
            "omega_c22_overintegration_debit_points": row.get("omega_c22_overintegration_debit_points"),
        }
        sample_rows.append(sample_key)
        for code in DIAGNOSTIC_TARGET_CODES:
            slug = code.replace(":", "_")
            peak_rows.append({
                **{key: sample_key[key] for key in ["date", "sample_no", "sample_name", "diagnostic_bucket"]},
                "target_code": code,
                "expected_rt": row.get(f"{slug}_expected_rt"),
                "corrected_target_rt": row.get(f"{slug}_corrected_target_rt"),
                "found_rt": row.get(f"{slug}_found_rt"),
                "rt_error": row.get(f"{slug}_rt_error"),
                "area": row.get(f"{slug}_area"),
                "width": row.get(f"{slug}_width"),
                "left_width": row.get(f"{slug}_left_width"),
                "right_width": row.get(f"{slug}_right_width"),
                "asymmetry": row.get(f"{slug}_asymmetry"),
                "status": row.get(f"{slug}_status"),
            })

    sample_diagnostics = pd.DataFrame(sample_rows)
    peak_diagnostics = pd.DataFrame(peak_rows)
    issue_summary = (
        sample_diagnostics.groupby("diagnostic_bucket", dropna=False)
        .agg(
            n=("diagnostic_bucket", "size"),
            MAE=("abs_delta", "mean"),
            RMSE=("delta", lambda values: math.sqrt(float(np.mean(np.square(values))))),
            max_abs=("abs_delta", "max"),
            within_0_3=("abs_delta", lambda values: int((values <= 0.3).sum())),
            within_0_5=("abs_delta", lambda values: int((values <= 0.5).sum())),
        )
        .reset_index()
        .sort_values(["max_abs", "MAE"], ascending=[False, False])
    )
    return sample_diagnostics, peak_diagnostics, issue_summary

def _markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_empty_"
    cols = [str(col) for col in df.columns]

    def fmt(value) -> str:
        if pd.isna(value):
            return ""
        if isinstance(value, float):
            return f"{value:.6f}"
        return str(value).replace("|", "\\|")

    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join(["---"] * len(cols)) + " |"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(fmt(row[col]) for col in df.columns) + " |")
    return "\n".join(lines)


def write_reports(
    out: Path,
    results: pd.DataFrame,
    summary: pd.DataFrame,
    outliers: pd.DataFrame,
    audit: pd.DataFrame,
    errors: pd.DataFrame,
    variants: pd.DataFrame,
    command_text: str,
) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    sample_diagnostics, peak_diagnostics, issue_summary = build_diagnostic_tables(results)
    with pd.ExcelWriter(out) as writer:
        results.to_excel(writer, sheet_name="Results", index=False)
        summary.to_excel(writer, sheet_name="Summary", index=False)
        outliers.to_excel(writer, sheet_name="Outliers_gt_0_5", index=False)
        audit.to_excel(writer, sheet_name="Input_audit", index=False)
        errors.to_excel(writer, sheet_name="Errors", index=False)
        sample_diagnostics.to_excel(writer, sheet_name="Sample_diagnostics", index=False)
        peak_diagnostics.to_excel(writer, sheet_name="Peak_diagnostics", index=False)
        issue_summary.to_excel(writer, sheet_name="Issue_summary", index=False)
        if not variants.empty:
            variants.to_excel(writer, sheet_name="Variants", index=False)

    stem = out.with_suffix("")
    prefix = "omega_regression" if out.name == "omega_regression_current.xlsx" else out.stem

    summary.to_csv(stem.parent / f"{prefix}_summary.csv", index=False)
    outliers.to_csv(stem.parent / f"{prefix}_outliers_gt_0_5.csv", index=False)
    audit.to_csv(stem.parent / f"{prefix}_input_audit.csv", index=False)
    errors.to_csv(stem.parent / f"{prefix}_errors.csv", index=False)
    if not variants.empty:
        variants.to_csv(stem.parent / f"{prefix}_variants.csv", index=False)
    sample_diagnostics.to_csv(stem.parent / f"{prefix}_sample_diagnostics.csv", index=False)
    peak_diagnostics.to_csv(stem.parent / f"{prefix}_peak_diagnostics.csv", index=False)
    issue_summary.to_csv(stem.parent / f"{prefix}_issue_summary.csv", index=False)

    all_row = summary[summary["scope"] == "ALL"].iloc[0].to_dict() if not summary.empty else {"n": 0}
    review_summary = (
        results.groupby(["review_flag", "outlier_class"], dropna=False)
        .size()
        .reset_index(name="n")
        .sort_values(["review_flag", "n"], ascending=[True, False])
        if {"review_flag", "outlier_class"}.issubset(results.columns)
        else pd.DataFrame()
    )
    report_lines = [
        "# Omega regression report",
        "",
        f"Generated with: `{command_text}`",
        "",
        f"Total evaluated samples: {int(all_row.get('n', 0))}",
        f"Overall MAE: {float(all_row.get('MAE', np.nan)):.6f}" if "MAE" in all_row else "Overall MAE: n/a",
        f"Overall RMSE: {float(all_row.get('RMSE', np.nan)):.6f}" if "RMSE" in all_row else "Overall RMSE: n/a",
        f"Overall max abs delta: {float(all_row.get('max_abs', np.nan)):.6f}" if "max_abs" in all_row else "Overall max abs delta: n/a",
        "",
        "## Summary",
        "",
        _markdown_table(summary),
        "",
        "## Input audit",
        "",
        _markdown_table(audit),
        "",
        "## Errors",
        "",
        _markdown_table(errors),
        "",
        "## Review / outlier classification",
        "",
        _markdown_table(review_summary),
        "",
        "## Diagnostic issue summary",
        "",
        _markdown_table(issue_summary),
        "",
        "## Top diagnostic samples",
        "",
        _markdown_table(sample_diagnostics.sort_values("abs_delta", ascending=False).head(25)[[col for col in ["date", "sample_no", "sample_name", "reference", "calculated", "delta", "confidence", "diagnostic_bucket", "diagnostic_reasons"] if col in sample_diagnostics.columns]]),
        "",
        "## Outliers > 0.5",
        "",
        f"Count: {len(outliers)}",
        "",
        _markdown_table(outliers[[col for col in ["date", "sample_no", "instrument_no", "sample_id", "sample_name", "match_method", "reference", "calculated", "delta", "confidence", "selected_variant"] if col in outliers.columns]]),
    ]
    (stem.parent / f"{prefix}_report.md").write_text("\n".join(report_lines), encoding="utf-8")


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Regression harness for the current Omega engine.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE_PATH)
    parser.add_argument("--out", type=Path, default=PROJECT_DIR / "omega_regression_current.xlsx")
    parser.add_argument("--selector-mode", choices=["current", "oracle"], default="current")
    parser.add_argument("--debug-dir", type=Path, default=None)
    parser.add_argument("--debug-threshold", type=float, default=0.5)
    args = parser.parse_args(list(argv) if argv is not None else None)

    results, audit, errors, variants = run_current_engine(
        args.data_dir,
        args.reference,
        selector_mode=args.selector_mode,
        debug_dir=args.debug_dir,
        debug_threshold=args.debug_threshold,
    )
    summary = build_summary_table(results)
    outliers = results[results["abs_delta"] > 0.5].sort_values("abs_delta", ascending=False) if "abs_delta" in results else pd.DataFrame()

    print(summary.to_string(index=False))
    if not audit.empty:
        print("\nInput audit")
        print(audit.to_string(index=False))
    if not errors.empty:
        print("\nErrors")
        print(errors.to_string(index=False))
    if not outliers.empty:
        print("\nOutliers > 0.5")
        display_cols = [col for col in ["date", "sample_no", "sample_name", "reference", "calculated", "delta", "confidence", "selected_variant"] if col in outliers]
        print(outliers[display_cols].to_string(index=False))

    if args.out:
        write_reports(
            args.out,
            results,
            summary,
            outliers,
            audit,
            errors,
            variants,
            command_text=" ".join(
                [
                    "python omega_regression.py",
                    f"--data-dir {args.data_dir}",
                    *(["--selector-mode oracle"] if args.selector_mode == "oracle" else []),
                    f"--out {args.out}",
                    *(["--debug-dir", str(args.debug_dir), "--debug-threshold", str(args.debug_threshold)] if args.debug_dir is not None else []),
                ]
            ),
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
