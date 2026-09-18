"""Compound-independent detection and integration of observed chromatographic peaks.

Derivatives locate maxima and valleys. Inflection points are diagnostics, not
integration limits. Original samples supply the integrals; a uniform time grid
is used only for geometry. Every neighbouring pair shares one valley decision.
"""
import numpy as np
import pandas as pd
from scipy.signal import find_peaks, savgol_filter
from scipy.ndimage import median_filter, percentile_filter

from .integration import validate_signal, integrate_interval, BoundaryProposal

SMOOTHING_MINUTES = 0.003
PERSISTENCE_MINUTES = 0.002
NOISE_WINDOW_MINUTES = 0.05
PROMINENCE_SIGMA = 5.0
BASELINE_SIGMA = 3.0
BASELINE_HEIGHT_FRACTION = 0.001
MAX_HALF_WIDTH_MINUTES = 0.25
GEOMETRY_METHOD = "derivative_valley_v1"


def derivative_signal(x, y, duration=SMOOTHING_MINUTES):
    x, y = validate_signal(x, y)
    step = float(np.median(np.diff(x)))
    count = int(round((x[-1]-x[0])/step)) + 1
    # An isolated dense run must not produce an unbounded interpolation grid.
    count = min(max(5, count), max(5, 4*len(x)))
    gx = np.linspace(x[0], x[-1], count)
    gy = np.interp(gx, x, y)
    step = float(gx[1]-gx[0])
    window = max(5, int(round(duration/step)) | 1)
    window = min(window, len(gx) if len(gx) % 2 else len(gx)-1)
    smooth = savgol_filter(gy, window, 3)
    dy = savgol_filter(gy, window, 3, deriv=1, delta=step)
    d2y = savgol_filter(gy, window, 3, deriv=2, delta=step)
    return gx, gy, smooth, dy, d2y, window


def detect_geometry(x, y):
    x, y = validate_signal(x, y)
    gx, gy, smooth, dy, d2y, window = derivative_signal(x, y)
    step = gx[1]-gx[0]
    # Second differences reject drift. A local MAD prevents a tall peak or a
    # quiet part of the run from determining the noise floor everywhere else.
    differences = np.r_[0., np.diff(gy, n=2), 0.]
    noise_window = max(7, int(round(NOISE_WINDOW_MINUTES/step)) | 1)
    noise_window = min(noise_window, len(gx) if len(gx) % 2 else len(gx)-1)
    center = median_filter(differences, size=noise_window, mode="reflect")
    noise = 1.4826*median_filter(abs(differences-center), size=noise_window, mode="reflect")/np.sqrt(6)
    residual_window = min(window*3, len(gx) if len(gx) % 2 else len(gx)-1)
    residual = gy-savgol_filter(gy, residual_window, 3)
    residual_center = median_filter(residual, size=noise_window, mode='reflect')
    residual_noise = 1.4826*median_filter(abs(residual-residual_center), size=noise_window, mode='reflect')
    noise = np.maximum(noise, residual_noise)
    epsilon = max(np.max(abs(gy))*1e-10, np.finfo(float).tiny)
    noise = np.maximum(noise, epsilon)
    # A densely occupied peak window contains curvature in the residual. Use
    # surrounding quiet-window estimates for baseline support, rather than
    # mistaking that curvature for noise and cutting the central peak's tails.
    baseline_noise = percentile_filter(noise, percentile=20, size=noise_window*3, mode='reflect')
    peaks, props = find_peaks(smooth, prominence=PROMINENCE_SIGMA*noise,
                             height=BASELINE_SIGMA*noise, width=2)
    # Require a maximum to survive a second smoothing scale. No fixed 0.03 min
    # distance rejects physically distinct narrow peaks or small neighbours.
    coarse_window = min((window*2-1), len(gx) if len(gx) % 2 else len(gx)-1)
    coarse = savgol_filter(gy, coarse_window, 3)
    coarse_peaks, _ = find_peaks(coarse, prominence=3*noise)
    supported = {}
    for pos, apex in enumerate(peaks):
        if not len(coarse_peaks):
            continue
        nearest = coarse_peaks[np.argmin(abs(coarse_peaks-apex))]
        if abs(gx[nearest]-gx[apex]) <= max(SMOOTHING_MINUTES, 2*step):
            previous = supported.get(int(nearest))
            if previous is None or props['prominences'][pos] > props['prominences'][previous]:
                supported[int(nearest)] = pos
    valid = sorted(supported.values())
    peaks = peaks[valid]
    prominences = props['prominences'][valid]
    valleys = [int(a+np.argmin(smooth[a:b+1])) for a, b in zip(peaks, peaks[1:])]
    valley_spreads = [abs(float(gx[v]-gx[a+np.argmin(coarse[a:b+1])]))
                      for a, b, v in zip(peaks, peaks[1:], valleys)]
    hold = max(2, int(np.ceil(PERSISTENCE_MINUTES/step)))
    records = []

    for position, apex in enumerate(peaks):
        height = float(smooth[apex])
        floor = max(BASELINE_SIGMA*baseline_noise[apex], BASELINE_HEIGHT_FRACTION*height)

        def edge(direction):
            neighbour = position-1 if direction < 0 else position+1
            if 0 <= neighbour < len(peaks):
                limit = valleys[min(position, neighbour)]
                kind = 'valley'
            else:
                limit = 0 if direction < 0 else len(gx)-1
                kind = 'edge'
            reach = max(2, int(round(MAX_HALF_WIDTH_MINUTES/step)))
            if abs(limit-apex) > reach:
                limit, kind = int(apex+direction*reach), 'limit'
            indices = np.arange(apex, limit+direction, direction)
            low = smooth[indices] <= floor
            sustained = np.convolve(low.astype(int), np.ones(hold, dtype=int), mode='valid') if len(low) >= hold else []
            hits = np.flatnonzero(np.asarray(sustained) == hold)
            if len(hits):
                return int(indices[hits[0]]), 'baseline'
            return limit, kind

        left, left_kind = edge(-1)
        right, right_kind = edge(1)
        if not left < apex < right:
            continue
        integral = integrate_interval(x, y, float(gx[left]), float(gx[right]))
        # Inflections bound the central lobe only; retain them for diagnostics.
        crossings = np.flatnonzero(np.diff(np.signbit(d2y[left:right+1])))+left
        left_inflections = crossings[crossings < apex]
        right_inflections = crossings[crossings >= apex]
        records.append(dict(
            peak_id=len(records)+1, start_idx=int(np.searchsorted(x, gx[left])),
            apex_idx=int(np.argmin(abs(x-gx[apex]))), end_idx=int(np.searchsorted(x, gx[right])),
            start_x=float(gx[left]), apex_x=float(gx[apex]), end_x=float(gx[right]),
            height=height, prominence=float(prominences[position]), width_points=float(right-left),
            area=integral.area, left_kind=left_kind, right_kind=right_kind,
            noise=float(noise[apex]), geometry_method=GEOMETRY_METHOD,
            valley_fraction=max(float(smooth[left]/height) if left_kind == 'valley' else 0.,
                                float(smooth[right]/height) if right_kind == 'valley' else 0.),
            boundary_uncertainty=max(valley_spreads[position-1] if position and left_kind == 'valley' else 0.,
                                     valley_spreads[position] if position < len(valleys) and right_kind == 'valley' else 0.),
            left_inflection=float(gx[left_inflections[-1]]) if len(left_inflections) else np.nan,
            right_inflection=float(gx[right_inflections[0]]) if len(right_inflections) else np.nan,
        ))
    columns = ['peak_id','start_idx','apex_idx','end_idx','start_x','apex_x','end_x',
               'height','prominence','width_points','area','percent_area','left_kind','right_kind',
               'noise','geometry_method','left_inflection','right_inflection','valley_fraction','boundary_uncertainty']
    out = pd.DataFrame(records, columns=columns)
    if len(out):
        out['percent_area'] = 100*out.area/out.area.sum()
    return out


def apply_geometry(matched, peaks):
    """One authoritative integration per observed assignment; manual/model stay explicit."""
    from .boundary_audit import snapshot, record_changes
    out = matched.copy(deep=True)
    out['geometry_method'] = GEOMETRY_METHOD
    before = snapshot(out)
    lookup = peaks.set_index('peak_id')
    for index, row in out.iterrows():
        status = str(row.get('status', ''))
        if 'manual' in status or 'fit' in status or 'estimated' in status:
            continue
        peak_id = row.get('matched_peak_id')
        if peak_id not in lookup.index:
            continue
        peak = lookup.loc[peak_id]
        out.loc[index, ['found_rt','integration_start_x','integration_end_x','area']] = [
            peak.apex_x, peak.start_x, peak.end_x, peak.area]
        for key in ('left_kind','right_kind','noise','geometry_method','left_inflection','right_inflection',
                    'valley_fraction','prominence','boundary_uncertainty'):
            out.at[index, key] = peak[key]
        if peak.left_kind not in ('baseline','valley') or peak.right_kind not in ('baseline','valley'):
            if 'unresolved' not in status:
                out.at[index, 'status'] = status+'_boundary_unresolved'
    out['percent_area'] = 100*out.area/out.area.fillna(0).sum()
    history = list(matched.attrs.get('boundary_history', []))
    record_changes(before, out, GEOMETRY_METHOD, history)
    out.attrs['boundary_history'] = history
    return out


def proposal_for_apex(peaks, apex):
    """The GUI proposes the same detector interval, with no second integrator."""
    if peaks.empty or not np.isfinite(apex):
        raise ValueError('Нет наблюдаемой вершины. Укажите границы на графике вручную.')
    peak = peaks.iloc[int(np.argmin(abs(peaks.apex_x-float(apex))))]
    if abs(peak.apex_x-apex) > .01 or not peak.start_x < apex < peak.end_x:
        raise ValueError('Рядом нет подтверждённой вершины. Укажите границы вручную.')
    return BoundaryProposal(float(peak.start_x), float(peak.end_x), float(peak.apex_x),
                            peak.left_kind, peak.right_kind, float(peak.noise))
