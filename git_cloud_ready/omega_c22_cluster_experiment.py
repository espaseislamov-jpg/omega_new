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


C22_LEFT = 9.20
C22_RIGHT = 9.35
C22_GRID_SIZE = 220
C22_CODES = ["C22:6", "C22:5", "C22:4"]


def _sample_no(sample_name: str) -> int | None:
    match = re.search(r"O(\d+)", sample_name or "")
    return int(match.group(1)) if match else None


def _normalized_c22_trace(processed: pd.DataFrame) -> np.ndarray:
    mask = (processed["x"] >= C22_LEFT) & (processed["x"] <= C22_RIGHT)
    segment = processed.loc[mask, ["x", "y_corrected"]].copy()
    if len(segment) < 8:
        return np.zeros(C22_GRID_SIZE, dtype=float)
    x = segment["x"].to_numpy(dtype=float)
    y = np.clip(segment["y_corrected"].to_numpy(dtype=float), 0.0, None)
    grid = np.linspace(C22_LEFT, C22_RIGHT, C22_GRID_SIZE)
    interpolated = np.interp(grid, x, y)
    scale = float(np.nanmax(interpolated))
    if not np.isfinite(scale) or scale <= 0:
        return np.zeros(C22_GRID_SIZE, dtype=float)
    return interpolated / scale


def _c22_shape_features(processed: pd.DataFrame, matched: pd.DataFrame) -> dict[str, float]:
    mask = (processed["x"] >= C22_LEFT) & (processed["x"] <= C22_RIGHT)
    x = processed.loc[mask, "x"].to_numpy(dtype=float)
    y = np.clip(processed.loc[mask, "y_corrected"].to_numpy(dtype=float), 0.0, None)
    if len(x) < 8 or not np.any(y > 0):
        return {
            "c22_peak_count": 0,
            "c22_width_at_10": np.nan,
            "c22_width_at_50": np.nan,
            "c22_left_mass": np.nan,
            "c22_right_mass": np.nan,
            "c22_centroid": np.nan,
            "dpa_to_c22_4": np.nan,
            "dha_to_dpa": np.nan,
            "dpa_percent": np.nan,
        }

    y_norm = y / max(float(np.max(y)), 1e-9)
    peaks, _ = find_peaks(y_norm, prominence=0.03, distance=4)
    above_10 = x[y_norm >= 0.10]
    above_50 = x[y_norm >= 0.50]
    total = float(np.trapezoid(y, x))
    left_mass = float(np.trapezoid(y[x <= 9.268], x[x <= 9.268])) / total if total > 0 and np.any(x <= 9.268) else np.nan
    right_mass = float(np.trapezoid(y[x >= 9.298], x[x >= 9.298])) / total if total > 0 and np.any(x >= 9.298) else np.nan
    centroid = float(np.trapezoid(x * y, x) / total) if total > 0 else np.nan

    area = {}
    percent = {}
    for code in C22_CODES:
        row = matched[matched["code"] == code]
        area[code] = float(row["area"].iloc[0]) if not row.empty and pd.notna(row["area"].iloc[0]) else np.nan
        percent[code] = float(row["percent_area"].iloc[0]) if not row.empty and pd.notna(row["percent_area"].iloc[0]) else np.nan

    return {
        "c22_peak_count": int(len(peaks)),
        "c22_width_at_10": float(above_10[-1] - above_10[0]) if len(above_10) > 1 else np.nan,
        "c22_width_at_50": float(above_50[-1] - above_50[0]) if len(above_50) > 1 else np.nan,
        "c22_left_mass": left_mass,
        "c22_right_mass": right_mass,
        "c22_centroid": centroid,
        "dpa_to_c22_4": area["C22:5"] / area["C22:4"] if area["C22:4"] and np.isfinite(area["C22:4"]) else np.nan,
        "dha_to_dpa": area["C22:6"] / area["C22:5"] if area["C22:5"] and np.isfinite(area["C22:5"]) else np.nan,
        "dpa_percent": percent["C22:5"],
    }


def _cluster_features(features: pd.DataFrame, traces: list[np.ndarray]) -> pd.DataFrame:
    trace_matrix = np.vstack(traces)
    scalar_columns = [
        "c22_peak_count",
        "c22_width_at_10",
        "c22_width_at_50",
        "c22_left_mass",
        "c22_right_mass",
        "c22_centroid",
        "dpa_to_c22_4",
        "dha_to_dpa",
        "dpa_percent",
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
            candidate = KMeans(n_clusters=n_clusters, random_state=7, n_init=20).fit_predict(scaled)
            try:
                score = silhouette_score(scaled, candidate)
            except Exception:
                score = -np.inf
            if score > best_score:
                best_score = score
                best_labels = candidate
        labels = best_labels if best_labels is not None else labels

    out = features.copy()
    out["c22_cluster"] = labels
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
            row.update(_c22_shape_features(result["processed_df"], matched))
            rows.append(row)
            traces.append(_normalized_c22_trace(result["processed_df"]))

    features = pd.DataFrame(rows)
    if features.empty:
        return features
    return _cluster_features(features, traces)


def summarize_by_cluster(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cluster_id, group in results.groupby("c22_cluster", sort=True):
        rows.append({
            "c22_cluster": cluster_id,
            "n": int(len(group)),
            "MAE": float(group["abs_delta"].mean()),
            "RMSE": float(math.sqrt(float(np.mean(np.square(group["delta"]))))),
            "mean_delta": float(group["delta"].mean()),
            "median_abs": float(group["abs_delta"].median()),
            "within_0_5": int((group["abs_delta"] <= 0.5).sum()),
            "outliers_gt_0_5": int((group["abs_delta"] > 0.5).sum()),
            "max_abs": float(group["abs_delta"].max()),
            "mean_dpa_to_c22_4": float(group["dpa_to_c22_4"].mean()),
            "mean_dha_to_dpa": float(group["dha_to_dpa"].mean()),
            "mean_peak_count": float(group["c22_peak_count"].mean()),
            "mean_width_at_10": float(group["c22_width_at_10"].mean()),
        })
    return pd.DataFrame(rows).sort_values(["outliers_gt_0_5", "MAE"], ascending=[False, False])


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


def summarize_overall(results: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame([{"scope": "ALL", **summarize_group(results)}])


def summarize_by_date(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for date, group in results.groupby("date", sort=True):
        rows.append({"date": date, **summarize_group(group)})
    return pd.DataFrame(rows)


def summarize_by_date_and_cluster(results: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (date, cluster_id), group in results.groupby(["date", "c22_cluster"], sort=True):
        rows.append({"date": date, "c22_cluster": cluster_id, **summarize_group(group)})
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="C22-only cluster experiment for clean ChromatoPy engine.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE_PATH)
    parser.add_argument("--out", type=Path, default=Path("omega_c22_cluster_experiment.xlsx"))
    args = parser.parse_args()

    results = build_dataset(args.data_dir, args.reference)
    if results.empty:
        raise SystemExit("No samples found.")
    overall = summarize_overall(results)
    by_date = summarize_by_date(results)
    summary = summarize_by_cluster(results)
    by_date_cluster = summarize_by_date_and_cluster(results)
    outliers = results[results["abs_delta"] > 0.5].sort_values("abs_delta", ascending=False)

    print("Overall summary")
    print(overall.to_string(index=False))
    print("\nBy date")
    print(by_date.to_string(index=False))
    print()
    print("C22 cluster summary")
    print(summary.to_string(index=False))
    print("\nOutliers > 0.5 by C22 cluster")
    print(outliers[[
        "c22_cluster",
        "date",
        "sample_no",
        "sample_name",
        "reference",
        "calculated",
        "delta",
        "dpa_to_c22_4",
        "dha_to_dpa",
        "omega_diagnostic",
    ]].to_string(index=False))

    with pd.ExcelWriter(args.out) as writer:
        results.to_excel(writer, sheet_name="Samples", index=False)
        overall.to_excel(writer, sheet_name="OverallSummary", index=False)
        by_date.to_excel(writer, sheet_name="ByDate", index=False)
        summary.to_excel(writer, sheet_name="ClusterSummary", index=False)
        by_date_cluster.to_excel(writer, sheet_name="ByDateAndCluster", index=False)
        outliers.to_excel(writer, sheet_name="Outliers_gt_0_5", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
