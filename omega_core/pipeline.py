"""One production path: baseline, geometry, ordered identity, evidence, result."""
from pathlib import Path
import pandas as pd
from . import instrument_profiles, io, matching, metrics, rt_profile, signal, runtime
from .peak_geometry import apply_geometry
from .peak_evidence import measurement_states
from .results import finalize_result
from .measurement import apply_area_baseline
from .epa_recovery import recover_missing_epa


def process_from_baseline(processed, reference_targets, strict_matching=False):
    processed, window = signal.add_smoothing_and_derivatives(processed)
    peaks = signal.detect_peak_candidates(processed)
    matched, shift = matching.match_targets_to_peaks(reference_targets, peaks, strict=strict_matching)
    matched = apply_geometry(matched, peaks)
    history = list(matched.attrs.get('boundary_history', []))
    matched = rt_profile.annotate_rt_profile(matched)
    matched = metrics.annotate_peak_heights(processed, matched)
    omega = metrics.compute_omega(matched)
    return dict(processed_df=processed, best_window=window, peaks_df=peaks,
                matched_targets_df=matched, rt_shift=shift, omega=omega,
                omega_report=omega['omega3_trio'], boundary_history=history,
                judge_decisions_df=pd.DataFrame(), separation_decisions=[])


def assignment_quality(result):
    """Baseline comparison uses evidence only, never an expected Omega value."""
    frame = result['matched_targets_df']
    states = measurement_states(frame)
    measured = states.isin(['observed', 'manual'])
    missing = int((~measured).sum())
    distances = pd.to_numeric(frame.loc[measured, 'match_score'], errors='coerce')
    return missing, float(distances.fillna(.035).sum())


def process_batch(dataframe: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    runtime.verify_versions()
    profile = instrument_profiles.upgrade_legacy_profile(instrument_profiles.profile_from_targets(reference_targets))
    reference_targets = instrument_profiles.apply_profile_to_targets(reference_targets, profile)
    # ASLS is the documented primary baseline. A single alternate baseline is
    # allowed when assignments are missing, replacing Omega-driven retries.
    if signal.ENABLE_ASLS_SHAPE_FALLBACK:
        if signal.Baseline is None:
            raise RuntimeError('PyBaselines не загрузился: расчёт остановлен.')
        processed, mode = signal.add_asls_baseline(dataframe), 'asls'
    else:
        processed, mode = signal.add_baseline(dataframe, **signal.BASELINE_KWARGS), 'chebyshev'
    result = process_from_baseline(processed, reference_targets)
    result['baseline_mode'] = mode
    decisions = []
    quality = assignment_quality(result)
    if quality[0] and mode == 'asls':
        alternative = process_from_baseline(signal.add_baseline(dataframe, **signal.BASELINE_KWARGS), reference_targets)
        alt_quality = assignment_quality(alternative)
        accepted = alt_quality < quality
        decisions.append(dict(stage='baseline_evidence', current=list(quality), candidate=list(alt_quality), accepted=accepted))
        if accepted:
            result = alternative
            result['baseline_mode'] = 'chebyshev_fallback'
    result['baseline_decisions'] = decisions
    result = apply_area_baseline(result, dataframe, profile.get('integration_baseline', 'geometry'))
    result = recover_missing_epa(result)
    result['analysis_environment'] = runtime.snapshot(reference_targets)
    return finalize_result(result, profile)


def process_file(file_path: Path, reference_path: Path = io.DEFAULT_REFERENCE_PATH,
                 cutoff_minutes: float = 4.0) -> list[dict]:
    targets = io.load_reference_targets(reference_path)
    return [{**batch, **process_batch(batch['dataframe'], targets)}
            for batch in io.load_batches(file_path, cutoff_minutes=cutoff_minutes)]
