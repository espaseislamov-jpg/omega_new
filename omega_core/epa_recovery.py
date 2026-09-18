"""Second pass for missing EPA supported by observed maxima and flanking acids.

Applied after the integration baseline is selected. Existing assignments are
untouched. Recovery requires review and never estimates a hidden component.
"""
import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.signal import find_peaks, savgol_filter

from . import peak_geometry as geometry
from .integration import integrate_interval, boundary_issues
from .boundary_audit import snapshot, record_changes
from .peak_evidence import measurement_states

RECOVERED_STATUS = 'recovered_epa_local_review'


def _nearest_positions(sorted_values, values):
    """Nearest entries in an ordered, nonempty array; ties use the first."""
    right = np.searchsorted(sorted_values, values)
    left = np.clip(right - 1, 0, len(sorted_values) - 1)
    right = np.minimum(right, len(sorted_values) - 1)
    return np.where(abs(sorted_values[right] - values) < abs(sorted_values[left] - values), right, left)


def _detection_trace(processed):
    gx,gy,smooth,dy,d2,window=geometry.derivative_signal(processed.x_corrected.to_numpy(),processed.y_corrected.to_numpy())
    step=gx[1]-gx[0];nw=max(7,int(round(.05/step))|1)
    nw=min(nw,len(gx) if len(gx)%2 else len(gx)-1)
    differences=np.r_[0.,np.diff(gy,n=2),0.]
    center=median_filter(differences,size=nw,mode='reflect')
    diff_noise=1.4826*median_filter(abs(differences-center),size=nw,mode='reflect')/np.sqrt(6)
    residual=gy-savgol_filter(gy,min(window*3,len(gx) if len(gx)%2 else len(gx)-1),3)
    rc=median_filter(residual,size=nw,mode='reflect')
    rn=1.4826*median_filter(abs(residual-rc),size=nw,mode='reflect')
    noise=np.maximum(np.maximum(diff_noise,rn),max(np.max(abs(gy))*1e-10,np.finfo(float).tiny))
    coarse=savgol_filter(gy,min(window*2-1,len(gx) if len(gx)%2 else len(gx)-1),3)
    cp,cprops=find_peaks(coarse,prominence=0)
    peaks,props=find_peaks(smooth,prominence=0,width=0)
    good_cp=cp[cprops['prominences']>=3*noise[cp]]
    # Reuse nearest coarse peaks instead of scanning every coarse peak several
    # times for every fine peak. Use time coordinates separately for distances
    # so floating-point grid rounding remains identical to the original search.
    nearest = _nearest_positions(cp, peaks) if len(cp) else None
    nearest_time = _nearest_positions(gx[cp], gx[peaks]) if len(cp) else None
    nearest_good = _nearest_positions(gx[good_cp], gx[peaks]) if len(good_cp) else None
    table=pd.DataFrame(dict(index=peaks,rt=gx[peaks],height=smooth[peaks],prominence=props['prominences'],
        noise=noise[peaks],diff_noise=diff_noise[peaks],residual_noise=rn[peaks],
        snr=props['prominences']/noise[peaks],width=props['widths'],
        coarse_snr=cprops['prominences'][nearest]/noise[cp[nearest]] if len(cp) else [0 for p in peaks],
        coarse_diff_snr=cprops['prominences'][nearest]/np.maximum(diff_noise[cp[nearest]],1e-10) if len(cp) else [0 for p in peaks],
        coarse_rt=gx[cp[nearest]] if len(cp) else np.full(len(peaks),np.nan),
        coarse_distance=abs(gx[good_cp[nearest_good]]-gx[peaks]) if len(good_cp) else np.full(len(peaks),np.nan),
        coarse_any_distance=abs(gx[cp[nearest_time]]-gx[peaks]) if len(cp) else np.full(len(peaks),np.nan)))
    return gx,gy,smooth,coarse,noise,table


def propose_missing_epa(processed, frame):
    f=frame.set_index('code')
    if 'C20:5' not in f.index:return dict(decision='missing_target')
    if not all(code in f.index for code in ('C20:4N6', 'C20:3N8')):
        return dict(decision='missing_anchor')
    if not f.index.is_unique:
        return dict(decision='ambiguous_targets')
    epa=f.loc['C20:5']
    if 'manual' in str(epa.status):
        return dict(decision='manual_target')
    if np.isfinite(epa.found_rt):return dict(decision='already_present')
    left,right=f.loc['C20:4N6'],f.loc['C20:3N8']
    if not np.isfinite([left.found_rt,right.found_rt]).all():return dict(decision='missing_anchor')
    if not measurement_states(frame.loc[frame.code.isin(['C20:4N6', 'C20:3N8'])]).eq('observed').all():
        return dict(decision='unconfirmed_anchor')
    anchors = np.array([left.expected_rt, epa.expected_rt, right.expected_rt,
                        left.found_rt, right.found_rt], dtype=float)
    if not np.isfinite(anchors).all() or not left.expected_rt < epa.expected_rt < right.expected_rt or not left.found_rt < right.found_rt:
        return dict(decision='invalid_order')
    # Correct local drift using flanking assigned acids, not the final Omega.
    fraction=(epa.expected_rt-left.expected_rt)/(right.expected_rt-left.expected_rt)
    predicted=float(left.found_rt+fraction*(right.found_rt-left.found_rt))
    if not 0<fraction<1:return dict(decision='invalid_order')
    gx,gy,fine,coarse,noise,table=_detection_trace(processed)
    step=float(gx[1]-gx[0])
    near=table[(abs(table.rt-predicted)<=.012)&(table.rt>left.found_rt+.005)&(table.rt<right.found_rt-.005)].copy()
    near['high_frequency_snr']=near.prominence/near.diff_noise.clip(lower=1e-8)
    viable=near[(near.high_frequency_snr>=5)&(near.coarse_diff_snr>=3)&
                (near.coarse_any_distance<=.003)&(near.width*step>=.0015)&
                (near.width*step<=.020)&(near.height>3*near.diff_noise)]
    # Fine maxima pointing to the same broad-scale maximum represent one lobe.
    viable=viable.sort_values('prominence',ascending=False).drop_duplicates('coarse_rt')
    if len(viable)!=1:
        return dict(decision='no_supported_peak' if not len(viable) else 'ambiguous_peaks',
                    predicted_rt=predicted,candidates=len(viable))
    candidate=viable.iloc[0];apex=int(candidate['index'])
    li=int(np.argmin(abs(gx-left.found_rt)));ri=int(np.argmin(abs(gx-right.found_rt)))
    a=li+int(np.argmin(coarse[li:apex+1]));b=apex+int(np.argmin(coarse[apex:ri+1]))
    if not li<a<apex<b<ri:return dict(decision='unresolved_boundary',predicted_rt=predicted)
    # Use shared valleys. The area is the displayed, unsmoothed corrected signal.
    # Do not use a model subtraction or a learned numerical area correction.
    integral=integrate_interval(processed.x_corrected.to_numpy(),processed.y_corrected.to_numpy(),gx[a],gx[b])
    if integral.area<=0:return dict(decision='nonpositive_area')
    return dict(decision='proposed_review',requires_review=True,predicted_rt=predicted,
                apex=float(gx[apex]),start=float(gx[a]),end=float(gx[b]),area=integral.area,
                high_frequency_snr=float(candidate.high_frequency_snr),
                conservative_snr=float(candidate.snr),coarse_snr=float(candidate.coarse_snr),
                coarse_high_frequency_snr=float(candidate.coarse_diff_snr),
                noise=float(candidate.diff_noise),height=float(candidate.height),
                prominence=float(candidate.prominence),coarse_rt=float(candidate.coarse_rt))


def recover_missing_epa(result):
    """Atomically insert a supported peak and trim overlapping flanking intervals."""
    proposal = propose_missing_epa(result['processed_df'], result['matched_targets_df'])
    if proposal['decision'] != 'proposed_review':
        return result
    frame = result['matched_targets_df'].copy(deep=True)
    before = snapshot(frame)
    processed = result['processed_df']
    x, y = processed.x_corrected.to_numpy(), processed.y_corrected.to_numpy()
    peaks = result['peaks_df'].copy(deep=True)
    ids = pd.concat([frame.matched_peak_id, peaks.peak_id]).dropna()
    peak_id = float(ids.max()) + 1 if len(ids) else 1.
    index = frame.index[frame.code == 'C20:5'][0]
    values = dict(found_rt=proposal['apex'], area=proposal['area'],
                  integration_start_x=proposal['start'], integration_end_x=proposal['end'],
                  matched_peak_id=peak_id, match_score=abs(proposal['apex']-proposal['predicted_rt']),
                  status=RECOVERED_STATUS, noise=proposal['noise'], prominence=proposal['prominence'],
                  left_kind='valley', right_kind='valley', geometry_method=geometry.GEOMETRY_METHOD)
    for key, value in values.items():
        frame.at[index, key] = value
    changed = ['C20:5']
    for code, key, bound in [('C20:4N6', 'integration_end_x', proposal['start']),
                             ('C20:3N8', 'integration_start_x', proposal['end'])]:
        j = frame.index[frame.code == code][0]
        overlaps = frame.at[j, key] > bound if key == 'integration_end_x' else frame.at[j, key] < bound
        if overlaps:
            frame.at[j, key] = bound
            a, b = frame.at[j, 'integration_start_x'], frame.at[j, 'integration_end_x']
            if not np.isfinite([a, b]).all() or not a < frame.at[j, 'found_rt'] < b:
                return result
            frame.at[j, 'area'] = integrate_interval(x, y, a, b).area
            frame.at[j, 'right_kind' if key == 'integration_end_x' else 'left_kind'] = 'valley'
            changed.append(code)
    # No new overlap with any measured acid is allowed, including unusual profiles.
    if any(issue['code'] in changed or issue.get('other_code') in changed
           for issue in boundary_issues(frame)):
        return result
    for code in changed:
        row = frame.loc[frame.code == code].iloc[0]
        record = dict(peak_id=row.matched_peak_id, start_x=row.integration_start_x,
                      apex_x=row.found_rt, end_x=row.integration_end_x, area=row.area,
                      start_idx=int(np.searchsorted(x, row.integration_start_x)),
                      apex_idx=int(np.argmin(abs(x-row.found_rt))),
                      end_idx=int(np.searchsorted(x, row.integration_end_x)),
                      height=float(np.interp(row.found_rt, x, y)),
                      width_points=float(np.searchsorted(x, row.integration_end_x)-np.searchsorted(x, row.integration_start_x)))
        for key in ('noise', 'prominence', 'left_kind', 'right_kind', 'geometry_method'):
            record[key] = row[key]
        existing = peaks.index[peaks.peak_id == row.matched_peak_id]
        if len(existing):
            for key, value in record.items():
                peaks.at[existing[0], key] = value
        else:
            peaks = pd.concat([peaks, pd.DataFrame([record])], ignore_index=True)
    frame['percent_area'] = 100*frame.area/frame.area.sum()
    peaks['percent_area'] = 100*peaks.area/peaks.area.sum()
    history = list(result.get('boundary_history', []))
    record_changes(before, frame, 'epa_local_recovery', history)
    frame.attrs['boundary_history'] = history
    return {**result, 'matched_targets_df': frame, 'peaks_df': peaks,
            'boundary_history': history, 'epa_recovery': proposal}

