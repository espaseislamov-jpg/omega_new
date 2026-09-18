"""Separate overlapping observed peaks at a reproducible shared valley."""
import numpy as np
from scipy.signal import savgol_filter, find_peaks

from .integration import integrate_interval, validate_signal, propose_boundaries
from . import boundary_audit


def recover_missing_oleic(processed, matched, peaks=None):
    """Recover an observed oleic maximum between two assigned C18 neighbours."""
    codes = ("C18:2N6C", "C18:1N9C", "C18:3N3")
    indices = [matched.index[matched.code == code].tolist() for code in codes]
    if any(len(i) != 1 for i in indices):
        return matched
    li, oi, ri = [i[0] for i in indices]
    oleic = matched.loc[oi]
    if "manual" in str(oleic.status) or (np.isfinite(oleic.found_rt) and np.isfinite(oleic.area) and oleic.area > 0
                                           and str(oleic.status) != "matched_c18_oleic_observed"):
        return matched
    left, right = matched.loc[li], matched.loc[ri]
    if any("manual" in str(row.status) or "fit" in str(row.status) for row in (left, right)):
        return matched
    a, b = float(left.found_rt), float(right.found_rt)
    if not np.isfinite([a, b]).all() or not .035 <= b-a <= .100:
        return matched
    x = processed["x_corrected" if "x_corrected" in processed else "x"].to_numpy()
    y = processed.y_corrected.to_numpy()
    validate_signal(x, y)
    step = float(np.median(np.diff(x)))
    mask = (x >= a-.015) & (x <= b+.015)
    sx, sy = x[mask], y[mask]
    window = max(5, int(round(.003/step)) | 1)
    if len(sx) <= window:
        return matched
    smooth = savgol_filter(sy, window, 2)
    residual = sy-smooth
    noise = 1.4826*np.median(abs(residual-np.median(residual)))
    candidates, _ = find_peaks(smooth, prominence=max(5*noise, .015*max(smooth)))
    candidates = [p for p in candidates if a+.008 < sx[p] < b-.008]
    if len(candidates) != 1:
        return matched
    apex = float(sx[candidates[0]])
    left_geometry = shared_valley(x, y, a, apex, details=True)
    right_geometry = shared_valley(x, y, apex, b, details=True)
    start = left_geometry["valley"] if left_geometry else None
    end = right_geometry["valley"] if right_geometry else None
    if start is None:
        return matched
    if end is None:
        proposal = propose_boundaries(x, y, apex, [a, b])
        if proposal.right_kind not in {"baseline", "valley"}:
            return matched
        end = proposal.end
    if not a < start < apex < end < b:
        return matched
    out = matched.copy(deep=True)
    integral = integrate_interval(x, y, start, end)
    out.loc[oi, ["found_rt", "integration_start_x", "integration_end_x", "area"]] = [apex, start, end, integral.area]
    out.at[oi, "status"] = "matched_c18_oleic_observed"
    out.at[oi, "matched_peak_id"] = np.nan
    out.at[oi, "match_score"] = np.nan
    out.at[li, "found_rt"] = left_geometry["left_apex"]
    if right_geometry:
        out.at[ri, "found_rt"] = right_geometry["right_apex"]
    if peaks is not None and not peaks.empty:
        near = peaks[abs(peaks.apex_x-apex) <= .004]
        used = set(out.loc[out.index != oi, "matched_peak_id"].dropna())
        near = near[~near.peak_id.isin(used)]
        if not near.empty:
            out.at[oi, "matched_peak_id"] = near.iloc[np.argmin(abs(near.apex_x-apex))].peak_id
    # Shared boundaries also recover tails cut off before the true valley.
    for index, edge, value in ((li, "integration_end_x", start), (ri, "integration_start_x", end)):
        row = out.loc[index]
        out.at[index, edge] = value
        out.at[index, "area"] = integrate_interval(x, y, out.at[index, "integration_start_x"], out.at[index, "integration_end_x"]).area
    out["percent_area"] = 100*out.area/out.area.fillna(0).sum()
    return out


def shared_valley(x, y, left_apex, right_apex, details=False):
    x, y = validate_signal(x, y)
    if not x[0] < left_apex < right_apex < x[-1]:
        return None
    step = float(np.median(np.diff(x)))
    padding = .012
    first = int(np.ceil((max(x[0], left_apex-padding)-x[0])/step))
    last = int(np.floor((min(x[-1], right_apex+padding)-x[0])/step))
    gx = x[0]+np.arange(first, last+1)*step
    gy = np.interp(gx, x, y)
    lo, hi = [int(np.argmin(abs(gx-v))) for v in (left_apex, right_apex)]
    if hi-lo < 5:
        return None
    candidates, apex_pairs = [], []
    for duration in (.0015, .003, .006):
        window = max(5, int(round(duration/step)) | 1)
        if window >= len(gx):
            return None
        smooth = savgol_filter(gy, window, 2)
        residual = gy-smooth
        noise = 1.4826*np.median(abs(residual-np.median(residual)))
        # RT matching can put an apex on a flank. Locate an observed maximum
        # near each assignment before asking whether a valley separates them.
        maxima, _ = find_peaks(smooth)
        # Disjoint neighbourhoods preserve peak order while tolerating an
        # assigned RT on a flank (e.g. O1's left apex is displaced by 0.008 min).
        radius = min(.015, .45*(right_apex-left_apex))
        refined = []
        for expected in (left_apex, right_apex):
            nearby = maxima[abs(gx[maxima]-expected) <= radius]
            if not len(nearby):
                return None
            refined.append(int(nearby[np.argmax(smooth[nearby])]))
        lo, hi = refined
        if hi-lo < 5:
            return None
        valley = lo+int(np.argmin(smooth[lo:hi+1]))
        height = min(smooth[lo], smooth[hi])
        if not lo+1 < valley < hi-1 or height <= 8*noise:
            return None
        if height-smooth[valley] <= max(4*noise, .02*height):
            return None
        # A shoulder along a monotone slope is not a second observed apex.
        reach = max(2, int(np.ceil(.0015/step)))
        for apex in (lo, hi):
            a, b = max(0, apex-reach), min(len(gx)-1, apex+reach)
            if smooth[apex] <= max(smooth[a], smooth[b]):
                return None
        candidates.append(float(gx[valley]))
        apex_pairs.append([float(gx[lo]), float(gx[hi])])
    if np.ptp(candidates) > max(2*step, .0015):
        return None
    valley = float(np.median(candidates))
    if details:
        apices = np.median(apex_pairs, axis=0)
        return {"valley": valley, "left_apex": float(apices[0]), "right_apex": float(apices[1])}
    return valley


def separate_overlaps(processed, matched):
    """Preserve fitted/manual components; partition observed overlaps once."""
    out = matched.copy(deep=True)
    before = boundary_audit.snapshot(out)
    x = processed["x_corrected" if "x_corrected" in processed else "x"].to_numpy()
    y = processed["y_corrected"].to_numpy()
    validate_signal(x, y)
    records = []
    for index, row in out.iterrows():
        values = np.asarray([row.get(k, np.nan) for k in
                             ("integration_start_x", "found_rt", "integration_end_x", "area")], dtype=float)
        if not np.isfinite(values).all() or values[3] <= 0 or not values[0] < values[1] < values[2]:
            continue
        records.append((values[1], index))
    records.sort()
    decisions, changed = [], set()
    for (_, left), (_, right) in zip(records, records[1:]):
        a, b = out.loc[left], out.loc[right]
        if a.integration_end_x <= b.integration_start_x+1e-10:
            continue
        statuses = [str(a.status), str(b.status)]
        if any("fit" in s or "manual" in s or "not_found" in s for s in statuses):
            continue
        geometry = shared_valley(x, y, float(a.found_rt), float(b.found_rt), details=True)
        valley = geometry["valley"] if geometry is not None else None
        decision = {"left": str(a.code), "right": str(b.code), "valley": valley,
                    "status": "separated" if valley is not None else "unresolved"}
        decisions.append(decision)
        if valley is None:
            for index in (left, right):
                if "unresolved" not in str(out.at[index, "status"]):
                    out.at[index, "status"] = str(out.at[index, "status"])+"_overlap_unresolved"
            continue
        out.at[left, "integration_end_x"] = valley
        out.at[right, "integration_start_x"] = valley
        out.at[left, "found_rt"] = geometry["left_apex"]
        out.at[right, "found_rt"] = geometry["right_apex"]
        changed.update((left, right))
    for index in changed:
        row = out.loc[index]
        integral = integrate_interval(x, y, float(row.integration_start_x), float(row.integration_end_x))
        out.at[index, "area"] = integral.area
        out.at[index, "status"] = str(row.status)+"_valley_separated"
    total = out.area.fillna(0).sum()
    if changed:
        out["percent_area"] = 100*out.area/total if total > 0 else np.nan
    history = list(matched.attrs.get("boundary_history", []))
    boundary_audit.record_changes(before, out, "observed_peak_separation", history)
    out.attrs["boundary_history"] = history
    return out, decisions
