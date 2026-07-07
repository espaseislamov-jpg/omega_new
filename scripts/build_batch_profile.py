#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

import pandas as pd

from omega_core.batch_profile import blend_profiles, build_profile, profile_to_frame, save_profile


def _markdown_table(df: pd.DataFrame) -> str:
    if df.empty:
        return "_empty_"
    def fmt(value):
        if pd.isna(value):
            return ""
        if isinstance(value, float):
            return f"{value:.6f}"
        return str(value)
    lines = ["| " + " | ".join(map(str, df.columns)) + " |", "| " + " | ".join(["---"] * len(df.columns)) + " |"]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(fmt(row[col]) for col in df.columns) + " |")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build robust RT/boundary profile from regression peak diagnostics.")
    parser.add_argument("--peak-diagnostics", type=Path, default=Path("regression_outputs/omega_regression_peak_diagnostics.csv"))
    parser.add_argument("--sample-diagnostics", type=Path, default=Path("regression_outputs/omega_regression_sample_diagnostics.csv"))
    parser.add_argument("--dates", nargs="*", default=["02072026", "03072026"])
    parser.add_argument("--out-json", type=Path, default=Path("regression_outputs/batch_profiles/latest_rt_boundary_profile.json"))
    parser.add_argument("--out-md", type=Path, default=Path("regression_outputs/batch_profiles/latest_rt_boundary_profile.md"))
    parser.add_argument("--previous-json", type=Path, default=None, help="Optional previous profile JSON to blend as warm-start.")
    parser.add_argument("--blend-alpha", type=float, default=0.35, help="Current-batch weight when --previous-json is used.")
    args = parser.parse_args()

    peaks = pd.read_csv(args.peak_diagnostics)
    peaks["date"] = peaks["date"].astype(str).str.replace(r"\.0$", "", regex=True).str.zfill(8)
    selected = peaks[peaks["date"].isin(args.dates)].copy() if args.dates else peaks
    profile = build_profile(selected)
    if args.previous_json is not None and args.previous_json.exists():
        previous = json.loads(args.previous_json.read_text(encoding="utf-8"))
        profile = blend_profiles(previous, profile, alpha=args.blend_alpha)
    save_profile(profile, args.out_json)
    profile_frame = profile_to_frame(profile)

    report = [
        "# Latest batch RT/boundary profile",
        "",
        f"Source peak diagnostics: `{args.peak_diagnostics}`",
        f"Dates: `{', '.join(args.dates) if args.dates else 'ALL'}`",
        f"Rows: `{profile.get('global', {}).get('n_rows')}`; samples: `{profile.get('global', {}).get('n_samples')}`",
        f"Blending: `{profile.get('blending', 'none')}`",
        "",
        "## Target profile summary",
        "",
        _markdown_table(profile_frame),
    ]

    if args.sample_diagnostics.exists():
        samples = pd.read_csv(args.sample_diagnostics)
        samples["date"] = samples["date"].astype(str).str.replace(r"\.0$", "", regex=True).str.zfill(8)
        samples = samples[samples["date"].isin(args.dates)].copy() if args.dates else samples
        high = samples.sort_values("abs_delta", ascending=False).head(16)
        report.extend([
            "",
            "## Highest-error samples in profiled dates",
            "",
            _markdown_table(high[[col for col in ["date", "sample_name", "reference", "calculated", "delta", "abs_delta", "safety_judge_band", "diagnostic_bucket"] if col in high.columns]]),
        ])

    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.write_text("\n".join(report), encoding="utf-8")
    print(f"Wrote {args.out_json}")
    print(f"Wrote {args.out_md}")
    print(profile_frame.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
