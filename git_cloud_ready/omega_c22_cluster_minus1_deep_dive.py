from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import omega_chromatopy_clean as clean
from omega_c22_cluster_experiment import C22_LEFT, C22_RIGHT, build_dataset
from omega_regression import DEFAULT_DATA_DIR, DEFAULT_REFERENCE_PATH


C22_CODES = ["C22:6", "C22:5", "C22:4"]


def _sample_no(sample_name: str) -> int | None:
    match = re.search(r"O(\d+)", sample_name or "")
    return int(match.group(1)) if match else None


def _load_batch_lookup(data_dir: Path) -> dict[tuple[str, int], dict]:
    lookup = {}
    config = clean.IntegrationConfig(use_chromatopy_fit=False)
    for csv_path in sorted(data_dir.glob("*.CSV")):
        date = csv_path.stem
        if not re.fullmatch(r"\d{8}", date):
            continue
        try:
            batches = clean.load_batches(csv_path, cutoff_minutes=config.cutoff_minutes)
        except Exception:
            continue
        for batch in batches:
            sample_no = _sample_no(batch.get("sample_name", ""))
            if sample_no is not None:
                lookup[(date, sample_no)] = batch
    return lookup


def _c22_rows(matched: pd.DataFrame) -> pd.DataFrame:
    rows = matched[matched["code"].isin(C22_CODES)].copy()
    rows["area"] = pd.to_numeric(rows["area"], errors="coerce")
    rows["integration_start_x"] = pd.to_numeric(rows["integration_start_x"], errors="coerce")
    rows["integration_end_x"] = pd.to_numeric(rows["integration_end_x"], errors="coerce")
    rows["found_rt"] = pd.to_numeric(rows["found_rt"], errors="coerce")
    rows["width"] = rows["integration_end_x"] - rows["integration_start_x"]
    return rows


def _local_min_between(processed: pd.DataFrame, left_rt: float, right_rt: float) -> tuple[float, float]:
    mask = (processed["x"] >= left_rt) & (processed["x"] <= right_rt)
    segment = processed.loc[mask, ["x", "y_corrected"]]
    if segment.empty:
        return np.nan, np.nan
    idx = int(segment["y_corrected"].idxmin())
    return float(processed.at[idx, "x"]), float(processed.at[idx, "y_corrected"])


def _analyze_sample(row: pd.Series, batch: dict, reference: pd.DataFrame) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    config = clean.IntegrationConfig(use_chromatopy_fit=False)
    result = clean.integrate_batch(batch["dataframe"], reference, config)
    processed = result["processed_df"]
    matched = result["matched_targets_df"]
    c22 = _c22_rows(matched)
    area_by_code = c22.set_index("code")["area"].to_dict()
    width_by_code = c22.set_index("code")["width"].to_dict()
    rt_by_code = c22.set_index("code")["found_rt"].to_dict()
    boundary_by_code = c22.set_index("code")[["integration_start_x", "integration_end_x"]].to_dict("index")

    dha = float(area_by_code.get("C22:6", np.nan))
    dpa = float(area_by_code.get("C22:5", np.nan))
    c22_4 = float(area_by_code.get("C22:4", np.nan))
    v_dha_dpa_x, v_dha_dpa_y = _local_min_between(processed, float(rt_by_code.get("C22:6", 9.25)), float(rt_by_code.get("C22:5", 9.28)))
    v_dpa_224_x, v_dpa_224_y = _local_min_between(processed, float(rt_by_code.get("C22:5", 9.28)), float(rt_by_code.get("C22:4", 9.31)))

    detail = {
        "date": row["date"],
        "sample_no": row["sample_no"],
        "sample_name": row["sample_name"],
        "reference": row["reference"],
        "calculated": row["calculated"],
        "delta": row["delta"],
        "abs_delta": row["abs_delta"],
        "c22_cluster": row["c22_cluster"],
        "diagnostic": row.get("omega_diagnostic", ""),
        "dha_area": dha,
        "dpa_area": dpa,
        "c22_4_area": c22_4,
        "dpa_to_c22_4": dpa / c22_4 if c22_4 > 0 else np.nan,
        "dha_to_dpa": dha / dpa if dpa > 0 else np.nan,
        "dha_width": width_by_code.get("C22:6", np.nan),
        "dpa_width": width_by_code.get("C22:5", np.nan),
        "c22_4_width": width_by_code.get("C22:4", np.nan),
        "dha_rt": rt_by_code.get("C22:6", np.nan),
        "dpa_rt": rt_by_code.get("C22:5", np.nan),
        "c22_4_rt": rt_by_code.get("C22:4", np.nan),
        "valley_dha_dpa_x": v_dha_dpa_x,
        "valley_dha_dpa_y": v_dha_dpa_y,
        "valley_dpa_c22_4_x": v_dpa_224_x,
        "valley_dpa_c22_4_y": v_dpa_224_y,
        "dha_start": boundary_by_code.get("C22:6", {}).get("integration_start_x", np.nan),
        "dha_end": boundary_by_code.get("C22:6", {}).get("integration_end_x", np.nan),
        "dpa_start": boundary_by_code.get("C22:5", {}).get("integration_start_x", np.nan),
        "dpa_end": boundary_by_code.get("C22:5", {}).get("integration_end_x", np.nan),
        "c22_4_start": boundary_by_code.get("C22:4", {}).get("integration_start_x", np.nan),
        "c22_4_end": boundary_by_code.get("C22:4", {}).get("integration_end_x", np.nan),
    }
    return detail, processed, c22


def _plot_c22_panel(details: pd.DataFrame, processed_map: dict[tuple[str, int], pd.DataFrame], c22_map: dict[tuple[str, int], pd.DataFrame], out_path: Path) -> None:
    rows = details.head(18).copy()
    if rows.empty:
        return
    ncols = 3
    nrows = int(np.ceil(len(rows) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(16, 4.2 * nrows), squeeze=False)
    for axis in axes.ravel():
        axis.axis("off")

    for axis, (_, item) in zip(axes.ravel(), rows.iterrows()):
        key = (str(item["date"]), int(item["sample_no"]))
        processed = processed_map[key]
        c22 = c22_map[key]
        segment = processed[(processed["x"] >= C22_LEFT) & (processed["x"] <= C22_RIGHT)]
        axis.axis("on")
        axis.plot(segment["x"], segment["y_corrected"], color="#1f4e79", linewidth=1.0)
        for _, peak in c22.iterrows():
            sx = float(peak["integration_start_x"])
            ex = float(peak["integration_end_x"])
            mask = (segment["x"] >= sx) & (segment["x"] <= ex)
            axis.fill_between(segment.loc[mask, "x"], 0, segment.loc[mask, "y_corrected"], alpha=0.25)
            axis.axvline(sx, color="#cc7a00", linewidth=0.7, alpha=0.8)
            axis.axvline(ex, color="#cc7a00", linewidth=0.7, alpha=0.8)
            axis.text(float(peak["found_rt"]), float(segment["y_corrected"].max()) * 0.75, str(peak["code"]), rotation=90, fontsize=8)
        axis.set_title(
            f"{item['date']} O{int(item['sample_no'])} d={item['delta']:+.2f} "
            f"DPA/C22:4={item['dpa_to_c22_4']:.2f}",
            fontsize=10,
        )
        axis.set_xlim(C22_LEFT, C22_RIGHT)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Deep dive for C22 cluster -1.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE_PATH)
    parser.add_argument("--cluster-report", type=Path, default=Path("omega_chromatopy_clean_clustered_stats.xlsx"))
    parser.add_argument("--out", type=Path, default=Path("omega_c22_cluster_minus1_deep_dive.xlsx"))
    parser.add_argument("--plot", type=Path, default=Path("omega_c22_cluster_minus1_panels.png"))
    args = parser.parse_args()

    if args.cluster_report.exists():
        samples = pd.read_excel(args.cluster_report, sheet_name="Samples")
    else:
        samples = build_dataset(args.data_dir, args.reference)
    minus1 = samples[samples["c22_cluster"] == -1].copy().sort_values("abs_delta", ascending=False)
    reference = clean.load_reference_targets(args.reference)
    lookup = _load_batch_lookup(args.data_dir)

    details = []
    processed_map = {}
    c22_map = {}
    for _, row in minus1.iterrows():
        key = (str(row["date"]), int(row["sample_no"]))
        batch = lookup.get(key)
        if batch is None:
            continue
        detail, processed, c22 = _analyze_sample(row, batch, reference)
        details.append(detail)
        processed_map[key] = processed
        c22_map[key] = c22

    detail_df = pd.DataFrame(details)
    high_dpa = detail_df[detail_df["dpa_to_c22_4"] > 1.15].copy().sort_values("dpa_to_c22_4", ascending=False)
    low_dpa = detail_df[detail_df["dpa_to_c22_4"] <= 1.15].copy().sort_values("dpa_to_c22_4", ascending=True)
    outliers = detail_df[detail_df["abs_delta"] > 0.5].copy().sort_values("abs_delta", ascending=False)
    aggregate = detail_df.describe(include="all")

    print("Cluster -1 detail")
    print(detail_df[[
        "date", "sample_no", "reference", "calculated", "delta", "dpa_to_c22_4", "dha_to_dpa",
        "dha_width", "dpa_width", "c22_4_width", "diagnostic",
    ]].to_string(index=False))
    print("\nHigh DPA/C22:4")
    print(high_dpa[["date", "sample_no", "delta", "dpa_to_c22_4", "dha_to_dpa", "diagnostic"]].to_string(index=False))

    with pd.ExcelWriter(args.out) as writer:
        detail_df.to_excel(writer, sheet_name="ClusterMinus1", index=False)
        outliers.to_excel(writer, sheet_name="Outliers", index=False)
        high_dpa.to_excel(writer, sheet_name="HighDPA", index=False)
        low_dpa.to_excel(writer, sheet_name="LowDPA", index=False)
        aggregate.to_excel(writer, sheet_name="Describe")

    _plot_c22_panel(outliers, processed_map, c22_map, args.plot)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
