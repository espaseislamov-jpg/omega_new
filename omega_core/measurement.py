"""Integrate displayed signal in the detector's fixed intervals."""
import numpy as np

from . import signal
from .integration import integrate_interval


def integrate_fixed_intervals(frame, processed):
    out = frame.copy(deep=True)
    x, y = processed.x_corrected.to_numpy(), processed.y_corrected.to_numpy()
    for idx, row in out.iterrows():
        a, apex, b = row.integration_start_x, row.found_rt, row.integration_end_x
        if np.isfinite([a, apex, b]).all() and x[0] <= a < apex < b <= x[-1]:
            out.at[idx, 'area'] = integrate_interval(x, y, a, b).area
        else:
            out.at[idx, 'area'] = np.nan
    total = out.area.sum()
    out['percent_area'] = 100 * out.area / total if total > 0 else np.nan
    return out


def apply_area_baseline(result, dataframe, method):
    """Geometry is already selected; area baseline cannot change assignments."""
    out = dict(result)
    out['geometry_baseline_mode'] = result['baseline_mode']
    if method == 'geometry':
        out['area_baseline_mode'] = result['baseline_mode']
        return out
    if method != 'chebyshev':
        raise ValueError(f'Неизвестная базовая линия интегрирования: {method}')
    processed = signal.add_baseline(dataframe, **signal.BASELINE_KWARGS)
    processed, window = signal.add_smoothing_and_derivatives(processed)
    out['processed_df'], out['best_window'] = processed, window
    out['matched_targets_df'] = integrate_fixed_intervals(result['matched_targets_df'], processed)
    peaks = result['peaks_df'].copy(deep=True)
    x, y = processed.x_corrected.to_numpy(), processed.y_corrected.to_numpy()
    for idx, row in peaks.iterrows():
        peaks.at[idx, 'area'] = integrate_interval(x, y, row.start_x, row.end_x).area
    total = peaks.area.sum()
    peaks['percent_area'] = 100 * peaks.area / total if total > 0 else np.nan
    out['peaks_df'] = peaks
    out['baseline_mode'] = 'chebyshev_fixed_g2_geometry'
    out['area_baseline_mode'] = 'chebyshev'
    out['calculation_method'] = 'fixed_geometry_direct_areas_g3'
    return out
