"""Ordered one-to-one RT matching, independent of peak integration and Omega."""
import numpy as np
import pandas as pd

from . import rt_profile

RT_WINDOW_MINUTES = .035
RT_SCALE_MINUTES = .012
MISSING_COST = 6.0


def match_ordered(targets, peaks):
    from .matching import estimate_rt_shift, estimate_profile_rt_calibration, _assign_peak
    out = targets.sort_values('order_index').reset_index(drop=True).copy()
    custom = rt_profile.uses_custom_instrument_profile(out)
    if not custom:
        # One effective RT table replaces coarse JSON + three cluster-specific
        # fallback tables. Preserve unknown codes' explicitly supplied RTs.
        known = out.code.map(rt_profile.MANUAL_TABLE_RTS)
        out['expected_rt'] = known.where(known.notna(), out.expected_rt)
    expected = out.expected_rt.to_numpy(dtype=float)
    peaks = peaks.sort_values('apex_x').reset_index(drop=True)
    apex = peaks.apex_x.to_numpy(dtype=float)
    selected = set()
    prominence = peaks.get('prominence', peaks.area).to_numpy(dtype=float)
    for rt in expected[np.isfinite(expected)]:
        nearby = np.flatnonzero(abs(apex-rt) <= .20)
        if len(nearby):
            selected.update(nearby[prominence[nearby] >= max(prominence[nearby])*.02].tolist())
    anchors = apex[sorted(selected)] if selected else apex
    finite = np.isfinite(expected)
    shift = estimate_rt_shift(expected[finite], anchors)
    slope, intercept = 1., shift
    if custom:
        _, slope, intercept = estimate_profile_rt_calibration(expected[finite], anchors)
    out['corrected_target_rt'] = expected*slope+intercept
    out['rt_calibration_slope'] = slope
    out['rt_calibration_intercept'] = intercept
    out['rt_local_shift'] = out.corrected_target_rt-out.expected_rt
    for name in ('found_rt','area','percent_area','matched_peak_id','match_score',
                 'integration_start_x','integration_end_x'):
        out[name] = np.nan
    out['status'] = 'not_found'
    # Dynamic programming allows skipping a missing target without stealing its
    # neighbour. Extra measured peaks are free to remain unidentified.
    n, m = len(out), len(peaks)
    costs = np.full((n+1, m+1), np.inf)
    actions = np.zeros((n+1, m+1), dtype=np.uint8)
    costs[0, :] = 0
    costs[:, 0] = np.arange(n+1)*MISSING_COST
    actions[1:, 0] = 1
    for i in range(1, n+1):
        rt = float(out.at[i-1, 'corrected_target_rt'])
        distances = abs(apex-rt)
        for j in range(1, m+1):
            options = [costs[i, j-1], costs[i-1, j]+MISSING_COST, np.inf]
            if np.isfinite(rt) and distances[j-1] <= RT_WINDOW_MINUTES:
                options[2] = costs[i-1, j-1]+(distances[j-1]/RT_SCALE_MINUTES)**2
            action = int(np.argmin(options))
            costs[i, j], actions[i, j] = options[action], action
    i, j = n, m
    while i > 0:
        action = actions[i, j]
        if action == 2:
            distance = abs(float(out.at[i-1, 'corrected_target_rt'])-apex[j-1])
            _assign_peak(out, i-1, peaks.iloc[j-1], 'matched_ordered_rt', distance)
            i, j = i-1, j-1
        elif action == 1 or j == 0:
            i -= 1
        else:
            j -= 1
    # A similarly plausible unused maximum must be visible as uncertainty.
    used = set(out.matched_peak_id.dropna())
    for index, row in out.loc[out.found_rt.notna()].iterrows():
        alternatives = peaks.loc[~peaks.peak_id.isin(used)]
        lower = out.loc[:index-1, 'found_rt'].max() if index else -np.inf
        upper = out.loc[index+1:, 'found_rt'].min() if index+1 < len(out) else np.inf
        lower = lower if np.isfinite(lower) else -np.inf
        upper = upper if np.isfinite(upper) else np.inf
        alternatives = alternatives.loc[(alternatives.apex_x > lower) & (alternatives.apex_x < upper)]
        if len(alternatives):
            distances = abs(alternatives.apex_x-float(row.corrected_target_rt))
            margin = float((distances.min()/RT_SCALE_MINUTES)**2-(row.match_score/RT_SCALE_MINUTES)**2)
            out.at[index, 'assignment_margin'] = margin
            if distances.min() <= RT_WINDOW_MINUTES and margin < .25:
                out.at[index, 'status'] = 'assignment_unresolved'
    return out, float(np.nanmedian(out.rt_local_shift)) if finite.any() else 0.
