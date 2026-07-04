from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import hdbscan
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

import omega_chromatopy_clean as clean
from omega_regression import DEFAULT_DATA_DIR, DEFAULT_REFERENCE_PATH, load_excel_refs


C20_LEFT = 8.34
C20_RIGHT = 8.52
C20_GRID_SIZE = 240
C20_CODES = ["C20:4N6", "C20:5", "C20:3N8"]


def _sample_no(sample_name: str) -> int | None:
    match = re.search(r"O(\d+)", sample_name or "")
    return int(match.group(1)) if match else None


def _normalized_c20_trace(processed: pd.DataFrame) -> np.ndarray:
    mask = (processed["x"] >= C20_LEFT) & (processed["x"] <= C20_RIGHT)
    segment = processed.loc[mask, ["x", "y_corrected"]].copy()
    if len(segment) < 8:
        return np.zeros(C20_GRID_SIZE, dtype=float)
    x = segment["x"].to_numpy(dtype=float)
    y = np.clip(segment["y_corrected"].to_numpy(dtype=float), 0.0, None)
    grid = np.linspace(C20_LEFT, C20_RIGHT, C20_GRID_SIZE)
    interpolated = np.interp(grid, x, y)
    scale = float(np.nanmax(interpolated))
    if not np.isfinite(scale) or scale <= 0:
        return np.zeros(C20_GRID_SIZE, dtype=float)
    return interpolated / scale


def _area_by_code(matched: pd.DataFrame) -> tuple[dict[str, float], dict[str, float]]:
    area: dict[str, float] = {}
    percent: dict[str, float] = {}
    for code in C20_CODES:
        row = matched[matched["code"] == code]
        area[code] = float(row["area"].iloc[0]) if not row.empty and pd.notna(row["area"].iloc[0]) else np.nan
        percent[code] = (
            float(row["percent_area"].iloc[0])
            if not row.empty and pd.notna(row["percent_area"].iloc[0])
            else np.nan
        )
    return area, percent


def _c20_shape_features(processed: pd.DataFrame, matched: pd.DataFrame) -> dict[str, float]:
    mask = (processed["x"] >= C20_LEFT) & (processed["x"] <= C20_RIGHT)
    x = processed.loc[mask, "x"].to_numpy(dtype=float)
    y = np.clip(processed.loc[mask, "y_corrected"].to_numpy(dtype=float), 0.0, None)
    if len(x) < 8 or not np.any(y > 0):
        return {
            "c20_peak_count": 0,
            "c20_width_at_10": np.nan,
            "c20_width_at_50": np.nan,
            "c20_left_mass": np.nan,
            "c20_mid_mass": np.nan,
            "c20_right_mass": np.nan,
            "c20_centroid": np.nan,
            "epa_to_c20_3": np.nan,
            "epa_to_c20_4": np.nan,
            "c20_3_to_c20_4": np.nan,
            "epa_percent": np.nan,
        }

    y_norm = y / max(float(np.max(y)), 1e-9)
    peaks, _ = find_peaks(y_norm, prominence=0.03, distance=4)
    above_10 = x[y_norm >= 0.10]
    above_50 = x[y_norm >= 0.50]
    total = float(np.trapezoid(y, x))
    left_mass = float(np.trapezoid(y[x <= 8.402], x[x <= 8.402])) / total if total > 0 and np.any(x <= 8.402) else np.nan
    mid_mask = (x >= 8.398) & (x <= 8.440)
    mid_mass = float(np.trapezoid(y[mid_mask], x[mid_mask])) / total if total > 0 and np.any(mid_mask) else np.nan
    right_mass = float(np.trapezoid(y[x >= 8.435], x[x >= 8.435])) / total if total > 0 and np.any(x >= 8.435) else np.nan
    centroid = float(np.trapezoid(x * y, x) / total) if total > 0 else np.nan
    area, percent = _area_by_code(matched)

    return {
        "c20_peak_count": int(len(peaks)),
        "c20_width_at_10": float(above_10[-1] - above_10[0]) if len(above_10) > 1 else np.nan,
        "c20_width_at_50": float(above_50[-1] - above_50[0]) if len(above_50) > 1 else np.nan,
        "c20_left_mass": left_mass,
        "c20_mid_mass": mid_mass,
        "c20_right_mass": right_mass,
        "c20_centroid": centroid,
        "epa_to_c20_3": area["C20:5"] / area["C20:3N8"] if area["C20:3N8"] and np.isfinite(area["C20:3N8"]) else np.nan,
        "epa_to_c20_4": area["C20:5"] / area["C20:4N6"] if area["C20:4N6"] and np.isfinite(area["C20:4N6"]) else np.nan,
        "c20_3_to_c20_4": area["C20:3N8"] / area["C20:4N6"] if area["C20:4N6"] and np.isfinite(area["C20:4N6"]) else np.nan,
        "epa_percent": percent["C20:5"],
    }


def _cluster_features(features: pd.DataFrame, traces: list[np.ndarray]) -> pd.DataFrame:
    trace_matrix = np.vstack(traces)
    scalar_columns = [
        "c20_peak_count",
        "c20_width_at_10",
        "c20_width_at_50",
        "c20_left_mass",
        "c20_mid_mass",
        "c20_right_mass",
        "c20_centroid",
        "epa_to_c20_3",
        "epa_to_c20_4",
        "c20_3_to_c20_4",
        "epa_percent",
    ]
    scalar = features[scalar_columns].replace([np.inf, -np.inf], np.nan)
    scalar = scalar.fillna(scalar.median(numeric_only=True)).fillna(0.0).to_numpy(dtype=float)
    x_features = np.hstack([trace_matrix, scalar])
    scaled = StandardScaler().fit_transform(x_features)

    labels = hdbscan.HDBSCAN(min_cluster_size=6, min_samples=3).fit_predict(scaled)
    if len(set(labels)) <= 1 or np.mean(labels == -1) > 0.55:
        best_labels = None
        best_score = -np.inf
        for n_clusters in range(3, 7):
            candidate = KMeans(n_clusters=n_clusters, random_state=11, n_init=20).fit_predict(scaled)
            try:
                score = silhouette_score(scaled, candidate)
            except Exception:
                score = -np.inf
            if score > best_score:
                best_score = score
                best_labels = candidate
        labels = best_labels if best_labels is not None else labels

    out = features.copy()
    out["c20_cluster"] = labels
    return out


def build_dataset(data_dir: Path, reference_path: Path) -> pd.DataFrame:
    reference = clean.load_reference_targets(reference_path)
    config = clean.IntegrationConfig(use_chromatopy_fit=False)
    rows = []
    traces: list[np.ndarray] = []

    for xlsx_path in sorted(data_dir.glob("test_bigbatch_*.xlsx")):
        date = xlsx_path.stem.split("_")[-1]
        csv_path = data_dir / f"{date}.CSV"
        if not csv_path.exists():
            continue
        refs = dict(load_excel_refs(xlsx_path))
        batches = clean.load_batches(csv_path, cutoff_minutes=config.cutoff_minutes)
        for batch in batches:
            sample_no = _sample_no(batch["sample_name"])
            if sample_no not in refs:
                continue
            result = clean.integrate_batch(batch["dataframe"], reference, config)
            matched = result["matched_targets_df"]
            calc = float(result["omega3_trio"])
            ref = float(refs[sample_no])
            row = {
                "date": date,
                "sample_no": sample_no,
                "sample_name": batch["sample_name"],
                "reference": ref,
                "calculated": calc,
                "delta": calc - ref,
                "abs_delta": abs(calc - ref),
                "n_matched": int(matched["area"].notna().sum()),
                "omega_diagnostic": ";".join(sorted({
                    str(value)
                    for value in matched.get("omega_diagnostic", pd.Series(dtype=str)).dropna()
                    if str(value).strip()
                })),
            }
            row.update(_c20_shape_features(result["processed_df"], matched))
            rows.append(row)
            traces.append(_normalized_c20_trace(result["processed_df"]))

    features = pd.DataFrame(rows)
    if features.empty:
        return features
    return _cluster_features(features, traces)


def summarize_group(frame: pd.DataFrame) -> dict[str, float | int]:
    if frame.empty:
        return {"n": 0}
    delta = frame["delta"]
    abs_delta = frame["abs_delta"]
    return {
        "n": int(len(frame)),
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
        "outliers_gt_0_5": int((abs_delta > 0.5).sum()),
        "max_abs": float(abs_delta.max()),
    }


def summarize_by_cluster(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cluster_id, group in results.groupby("c20_cluster", sort=True):
        rows.append({
            "c20_cluster": cluster_id,
            **summarize_group(group),
            "mean_epa_to_c20_3": float(group["epa_to_c20_3"].mean()),
            "mean_epa_to_c20_4": float(group["epa_to_c20_4"].mean()),
            "mean_peak_count": float(group["c20_peak_count"].mean()),
            "mean_width_at_10": float(group["c20_width_at_10"].mean()),
        })
    return pd.DataFrame(rows).sort_values(["outliers_gt_0_5", "MAE"], ascending=[False, False])


def summarize_by_date(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for date, group in results.groupby("date", sort=True):
        rows.append({"date": date, **summarize_group(group)})
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="C20/EPA cluster experiment for clean ChromatoPy engine.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE_PATH)
    parser.add_argument("--out", type=Path, default=Path("omega_c20_cluster_experiment.xlsx"))
    args = parser.parse_args()

    results = build_dataset(args.data_dir, args.reference)
    if results.empty:
        raise SystemExit("No samples found.")
    overall = pd.DataFrame([{"scope": "ALL", **summarize_group(results)}])
    by_date = summarize_by_date(results)
    summary = summarize_by_cluster(results)
    outliers = results[results["abs_delta"] > 0.5].sort_values("abs_delta", ascending=False)

    print("Overall summary")
    print(overall.to_string(index=False))
    print("\nBy date")
    print(by_date.to_string(index=False))
    print("\nC20 cluster summary")
    print(summary.to_string(index=False))
    print("\nOutliers > 0.5 by C20 cluster")
    print(outliers[[
        "c20_cluster",
        "date",
        "sample_no",
        "sample_name",
        "reference",
        "calculated",
        "delta",
        "epa_to_c20_3",
        "epa_to_c20_4",
        "omega_diagnostic",
    ]].to_string(index=False))

    with pd.ExcelWriter(args.out) as writer:
        results.to_excel(writer, sheet_name="Samples", index=False)
        overall.to_excel(writer, sheet_name="OverallSummary", index=False)
        by_date.to_excel(writer, sheet_name="ByDate", index=False)
        summary.to_excel(writer, sheet_name="ClusterSummary", index=False)
        outliers.to_excel(writer, sheet_name="Outliers_gt_0_5", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
