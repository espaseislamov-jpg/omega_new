"""Single final calculation path shared by batch processing and manual edits."""
from . import metrics, instrument_profiles
from .peak_evidence import assess_final_evidence


def finalize_result(result, profile):
    out = dict(result)
    out.pop("omega_estimate", None)
    out["matched_targets_df"] = metrics.annotate_peak_heights(out["processed_df"], out["matched_targets_df"])
    out["omega"] = metrics.compute_omega(out["matched_targets_df"])
    out["omega_report"] = out["omega"]["omega3_trio"]
    out = metrics.annotate_result(out, out.get("baseline_mode", "chebyshev"))
    out["confidence"] = assess_final_evidence(out["matched_targets_df"], out["confidence"])
    out = instrument_profiles.apply_result_multiplier(out, profile)
    return metrics.apply_reporting_gate(out)
