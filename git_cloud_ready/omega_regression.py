from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

import omega_core


PROJECT_DIR = Path(__file__).resolve().parent
DEFAULT_DATA_DIR = Path(r"C:\Users\marat\Desktop\CSV_Omega")
DEFAULT_REFERENCE_PATH = PROJECT_DIR / "reference_targets_reverted_c22fixed.json"
OLD45_DATES = {"13032026", "14012026", "20032026"}


def load_excel_refs(xlsx_path: Path) -> list[tuple[int, float]]:
    raw = pd.read_excel(xlsx_path, header=None)
    refs: list[tuple[int, float]] = []
    for row in raw.itertuples(index=False):
        sample_no = None
        ref_value = None
        for value in row:
            if pd.isna(value):
                continue
            if sample_no is None:
                text = str(value).strip()
                if text.isdigit():
                    sample_no = int(text)
                    continue
            if ref_value is None:
                try:
                    ref_value = float(str(value).strip().replace(",", "."))
                except ValueError:
                    pass
        if sample_no is not None and ref_value is not None:
            refs.append((sample_no, ref_value))
    return refs


def run_current_engine(
    data_dir: Path = DEFAULT_DATA_DIR,
    reference_path: Path = DEFAULT_REFERENCE_PATH,
) -> pd.DataFrame:
    reference_targets = omega_core.load_reference_targets(reference_path)
    rows = []

    for xlsx_path in sorted(data_dir.glob("test_bigbatch_*.xlsx")):
        date = xlsx_path.stem.split("_")[-1]
        csv_path = data_dir / f"{date}.CSV"
        if not csv_path.exists():
            continue

        refs = load_excel_refs(xlsx_path)
        batches = omega_core.load_batches(csv_path, cutoff_minutes=4.0)
        for sample_no, reference in refs:
            if sample_no < 1 or sample_no > len(batches):
                continue
            batch = batches[sample_no - 1]
            result = omega_core.process_batch(batch["dataframe"], reference_targets)
            calculated = float(result["omega_report"])
            delta = calculated - float(reference)
            confidence = result.get("confidence", {})
            confidence_score = confidence.get("score", np.nan) if isinstance(confidence, dict) else np.nan
            judge_decisions = result.get("judge_decisions_df")
            if isinstance(judge_decisions, pd.DataFrame) and not judge_decisions.empty:
                accepted_count = int((judge_decisions.get("decision") == "accepted").sum())
                rejected_count = int((judge_decisions.get("decision") == "rejected").sum())
                judge_codes = ",".join(sorted({str(code) for code in judge_decisions.get("code", pd.Series(dtype=str)).dropna()}))
                judge_reasons = ",".join(sorted({str(reason) for reason in judge_decisions.get("reason", pd.Series(dtype=str)).dropna()}))
            else:
                accepted_count = 0
                rejected_count = 0
                judge_codes = ""
                judge_reasons = ""
            rows.append({
                "date": date,
                "sample_no": sample_no,
                "sample_name": batch.get("sample_name", ""),
                "reference": float(reference),
                "calculated": calculated,
                "delta": delta,
                "abs_delta": abs(delta),
                "confidence": confidence_score,
                "baseline_mode": result.get("baseline_mode", ""),
                "cluster_quality_score": result.get("cluster_quality_score", np.nan),
                "judge_accepted": accepted_count,
                "judge_rejected": rejected_count,
                "judge_codes": judge_codes,
                "judge_reasons": judge_reasons,
            })

    return pd.DataFrame(rows)


def summarize(results: pd.DataFrame) -> dict[str, float | int]:
    if results.empty:
        return {"n": 0}
    abs_delta = results["abs_delta"]
    delta = results["delta"]
    return {
        "n": int(len(results)),
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
    old45 = results[results["date"].isin(OLD45_DATES)]
    rows.append({"scope": "OLD45", **summarize(old45)})
    for date, group in results.groupby("date", sort=True):
        rows.append({"scope": date, **summarize(group)})
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description="Regression harness for the current Omega engine.")
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE_PATH)
    parser.add_argument("--out", type=Path, default=PROJECT_DIR / "omega_regression_current.xlsx")
    args = parser.parse_args()

    results = run_current_engine(args.data_dir, args.reference)
    summary = build_summary_table(results)
    outliers = results[results["abs_delta"] > 0.5].sort_values("abs_delta", ascending=False)

    print(summary.to_string(index=False))
    if not outliers.empty:
        print("\nOutliers > 0.5")
        print(outliers[["date", "sample_no", "sample_name", "reference", "calculated", "delta", "confidence"]].to_string(index=False))

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with pd.ExcelWriter(args.out) as writer:
            results.to_excel(writer, sheet_name="Results", index=False)
            summary.to_excel(writer, sheet_name="Summary", index=False)
            outliers.to_excel(writer, sheet_name="Outliers_gt_0_5", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
