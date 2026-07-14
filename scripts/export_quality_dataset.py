from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

import omega_core
import omega_regression as regression


DEFAULT_DATA_DIR = PROJECT_DIR / "data" / "regression"
DEFAULT_REFERENCE_PATH = PROJECT_DIR / "reference_targets_reverted_c22fixed.json"
DEFAULT_OUTPUT = PROJECT_DIR / "artifacts" / "quality_dataset.csv"
SEALED_DATES = frozenset({"14072026"})
STATUS_TOKENS = (
    "not_found",
    "recovered",
    "estimated",
    "local",
    "fit",
    "overlap",
    "tail",
    "baseexpand",
    "split",
    "valley",
    "reassigned",
)


def _git_revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=PROJECT_DIR,
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except Exception:
        return "unknown"


def _finite_float(value, default=np.nan) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return float(default)
    return number if np.isfinite(number) else float(default)


def extract_quality_features(result: dict) -> dict[str, float | int | bool]:
    """Create reference-free numeric features available in the production GUI."""
    raw = {
        **regression._extract_result_features(result),
        **regression._judge_features(result),
    }
    features: dict[str, float | int | bool] = {}
    for key, value in raw.items():
        if isinstance(value, (bool, np.bool_)):
            features[key] = bool(value)
        elif isinstance(value, (int, float, np.integer, np.floating)):
            features[key] = _finite_float(value)

    baseline_mode = str(result.get("baseline_mode", "") or "")
    features["baseline_is_chebyshev"] = baseline_mode == "chebyshev"
    features["baseline_is_fallback"] = baseline_mode not in {"", "chebyshev"}

    matched = result.get("matched_targets_df")
    if not isinstance(matched, pd.DataFrame) or matched.empty:
        features.update({
            "targets_total": 0,
            "targets_found": 0,
            "targets_missing": 0,
            "matched_peak_ids_unique": 0,
            "matched_peak_id_duplicates": 0,
        })
        return features

    found = pd.to_numeric(matched.get("found_rt"), errors="coerce")
    peak_ids = pd.to_numeric(matched.get("matched_peak_id"), errors="coerce").dropna()
    unique_peak_ids = int(peak_ids.nunique())
    features["targets_total"] = int(len(matched))
    features["targets_found"] = int(found.notna().sum())
    features["targets_missing"] = int(found.isna().sum())
    features["matched_peak_ids_unique"] = unique_peak_ids
    features["matched_peak_id_duplicates"] = int(max(0, len(peak_ids) - unique_peak_ids))

    token_counts = {token: 0 for token in STATUS_TOKENS}
    for _, target in matched.iterrows():
        code = str(target.get("code", "")).replace(":", "_").replace("/", "_")
        status = str(target.get("status", "") or "").lower()
        if not code:
            continue
        for token in STATUS_TOKENS:
            present = token in status
            features[f"{code}_status_{token}"] = present
            token_counts[token] += int(present)
    for token, count in token_counts.items():
        features[f"status_count_{token}"] = count
    return features


def build_dataset(
    data_dir: Path,
    reference_path: Path,
    dates: set[str] | None = None,
    limit: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    requested_sealed = SEALED_DATES.intersection(dates or set())
    if requested_sealed:
        raise ValueError(f"Sealed dates cannot be exported: {sorted(requested_sealed)}")

    reference_targets = omega_core.load_reference_targets(reference_path)
    rows: list[dict] = []
    errors: list[dict] = []
    audit: list[dict] = []

    for xlsx_path, csv_path, raw_date, date in regression.discover_reference_pairs(data_dir):
        if date in SEALED_DATES:
            print(f"[sealed] skipping {date}", flush=True)
            continue
        if dates is not None and date not in dates:
            continue

        refs = regression.load_excel_refs(xlsx_path)
        batches = regression.build_batch_records(csv_path)
        audit_row = regression.build_audit_row(date, xlsx_path, csv_path, refs, batches)
        audit.append(audit_row)
        print(f"[{date}] references={len(refs)} chromatograms={len(batches)}", flush=True)

        for ref in refs:
            if limit is not None and len(rows) >= limit:
                break
            batch, match_method = regression.resolve_batch(ref, batches, date)
            if batch is None:
                errors.append({
                    "date": date,
                    "excel_row": ref.row_number,
                    "sample_token": ref.token,
                    "reference": ref.reference,
                    "match_method": match_method,
                    "error": "No matching chromatogram",
                })
                continue

            try:
                result = omega_core.process_batch(batch.dataframe, reference_targets)
                features = extract_quality_features(result)
                calculated = _finite_float(result.get("omega_report"))
                if not np.isfinite(calculated):
                    raise ValueError("Integrator returned a non-finite omega value")
                delta = calculated - float(ref.reference)
                rows.append({
                    "batch_date": date,
                    "raw_date": raw_date,
                    "sample_id": batch.sample_id,
                    "instrument_no": batch.instrument_no,
                    "sample_name": batch.sample_name,
                    "excel_row": ref.row_number,
                    "match_method": match_method,
                    "reference": float(ref.reference),
                    **features,
                    "delta": delta,
                    "abs_error": abs(delta),
                    "error_gt_0_3": int(abs(delta) > 0.3),
                    "error_gt_0_5": int(abs(delta) > 0.5),
                })
                print(
                    f"  {len(rows):04d} {batch.sample_name}: "
                    f"manual={ref.reference:.3f} calculated={calculated:.3f} abs={abs(delta):.3f}",
                    flush=True,
                )
            except Exception as exc:
                errors.append({
                    "date": date,
                    "excel_row": ref.row_number,
                    "sample_token": ref.token,
                    "sample_name": batch.sample_name,
                    "reference": ref.reference,
                    "match_method": match_method,
                    "error": f"{type(exc).__name__}: {exc}",
                })

        if limit is not None and len(rows) >= limit:
            break

    dataset = pd.DataFrame(rows)
    error_frame = pd.DataFrame(errors)
    manifest = {
        "git_revision": _git_revision(),
        "sealed_dates": sorted(SEALED_DATES),
        "rows": int(len(dataset)),
        "batches": sorted(dataset["batch_date"].astype(str).unique().tolist()) if not dataset.empty else [],
        "errors": int(len(error_frame)),
        "class_counts": {
            "error_gt_0_3": int(dataset.get("error_gt_0_3", pd.Series(dtype=int)).sum()),
            "error_gt_0_5": int(dataset.get("error_gt_0_5", pd.Series(dtype=int)).sum()),
        },
        "input_audit": audit,
    }
    return dataset, error_frame, manifest


def main() -> int:
    parser = argparse.ArgumentParser(description="Export reference-free integrator features for quality-model training.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE_PATH)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dates", nargs="*", default=None, help="Optional DDMMYYYY batch allow-list.")
    parser.add_argument("--limit", type=int, default=None, help="Smoke-test row limit; omit for the full dataset.")
    args = parser.parse_args()

    dates = {regression.normalize_date_token(value) for value in args.dates} if args.dates else None
    dataset, errors, manifest = build_dataset(args.data_dir, args.reference, dates=dates, limit=args.limit)
    if dataset.empty:
        raise SystemExit("No labeled samples were exported")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_csv(args.output, index=False)
    errors.to_csv(args.output.with_name(f"{args.output.stem}_errors.csv"), index=False)
    args.output.with_name(f"{args.output.stem}_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
