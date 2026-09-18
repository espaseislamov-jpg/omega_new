"""Headless, portable calculation report for comparing computers."""
import argparse
import json
import math
from pathlib import Path

from . import io, pipeline, runtime, instrument_profiles
from .manual_edits import file_digest


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if hasattr(value, "item"):
        return clean(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify-csv", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--samples", help="Comma-separated sample numbers, starting at 1; omit for all")
    parser.add_argument('--profile-file', type=Path, help='Exported instrument profile; default is the current 118 method')
    args = parser.parse_args(argv)
    selected = {int(v)-1 for v in args.samples.split(",")} if args.samples else None
    runtime.verify_versions()
    targets = io.load_reference_targets()
    if args.profile_file:
        raw = json.loads(args.profile_file.read_text(encoding='utf-8-sig'))
        profile = instrument_profiles.validate_profile(raw.get('profile',raw),targets.code)
    else:
        profile = instrument_profiles.legacy_profile(targets)
    targets = instrument_profiles.apply_profile_to_targets(targets,profile)
    report = {"input_sha256": file_digest(args.verify_csv), "runtime": runtime.snapshot(targets),
              "profile": profile, "samples": [], "errors": []}
    try:
        batches = io.load_batches(args.verify_csv)
        if selected is not None and (not selected or min(selected) < 0 or max(selected) >= len(batches)):
            raise ValueError("Sample number outside the file")
        for index, batch in enumerate(batches):
            if selected is not None and index not in selected:
                continue
            result = pipeline.process_batch(batch["dataframe"], targets)
            columns = [c for c in ("code", "area", "found_rt", "integration_start_x", "integration_end_x", "status",
                                   'geometry_method', 'left_kind', 'right_kind', 'valley_fraction',
                                   'boundary_uncertainty', 'assignment_margin')
                       if c in result["matched_targets_df"]]
            report["samples"].append({"index": index, "sample": batch["sample_name"],
                                      "omega_report": result["omega_report"],
                                      "omega_estimate": result.get("omega_estimate"),
                                      "report_status": result["report_status"], "reasons": result["report_reasons"],
                                      "separation": result.get("separation_decisions", []),
                                      'baseline_mode': result['baseline_mode'],
                                      'geometry_baseline_mode': result.get('geometry_baseline_mode'),
                                      'area_baseline_mode': result.get('area_baseline_mode'),
                                      'baseline_decisions': result.get('baseline_decisions', []),
                                      'confidence': result.get('confidence', {}),
                                      "peaks": result["matched_targets_df"][columns].to_dict("records")})
    except Exception as exc:
        report["errors"].append(f"{type(exc).__name__}: {exc}")
        raise
    finally:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(clean(report), ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
