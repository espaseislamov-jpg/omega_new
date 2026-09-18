"""Signal geometry in minutes, independent of compound names and Omega values."""
from dataclasses import dataclass

import numpy as np
from scipy.signal import savgol_filter


def validate_signal(x, y):
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    if x.ndim != 1 or y.shape != x.shape or len(x) < 2:
        raise ValueError("Сигнал должен содержать минимум две точки времени и интенсивности.")
    if not np.isfinite(x).all() or not np.isfinite(y).all() or np.any(np.diff(x) <= 0):
        raise ValueError("Время должно строго возрастать; NaN и бесконечности недопустимы.")
    return x, y


def interval_points(x, y, start, end):
    """Include exact fractional endpoints and zero crossings of linear segments."""
    x, y = validate_signal(x, y)
    if not np.isfinite([start, end]).all() or not x[0] <= start < end <= x[-1]:
        raise ValueError("Границы должны лежать внутри сигнала; начало меньше конца.")
    a, b = np.searchsorted(x, [start, end], side="right")
    interior = x[a:b]
    interior = interior[interior < end]
    sx = np.r_[start, interior, end]
    sy = np.interp(sx, x, y)
    crossing = np.flatnonzero(sy[:-1] * sy[1:] < 0)
    if len(crossing):
        zeros = sx[crossing] - sy[crossing] * np.diff(sx)[crossing] / np.diff(sy)[crossing]
        sx = np.sort(np.r_[sx, zeros])
        sy = np.interp(sx, x, y)
    return sx, sy


@dataclass(frozen=True)
class Integral:
    start: float
    end: float
    apex: float
    area: float
    baseline: str


def integrate_interval(x, y, start, end, baseline="zero"):
    """Integrate the positive signal with a zero or endpoint-to-endpoint baseline."""
    x, y = validate_signal(x, y)
    x, y = interval_points(x, y, start, end)
    if baseline == "linear":
        levels = np.interp([start, end], x, y)
        y = y - np.interp(x, [start, end], levels)
    elif baseline != "zero":
        raise ValueError("Неизвестный способ построения базовой линии.")
    sx, sy = interval_points(x, y, start, end)
    return Integral(float(start), float(end), float(sx[np.argmax(sy)]),
                    float(np.trapezoid(np.maximum(sy, 0), sx)), baseline)


@dataclass(frozen=True)
class BoundaryProposal:
    start: float
    end: float
    apex: float
    left_kind: str
    right_kind: str
    noise: float

    @property
    def resolved(self):
        return self.left_kind in {"baseline", "valley"} and self.right_kind in {"baseline", "valley"}


def propose_boundaries(x, y, apex, neighbors=(), max_half_width=0.12,
                       smoothing_minutes=0.003, persistence_minutes=0.002):
    """Walk from an observed apex to sustained baseline support or a shared valley.

    Work on a uniform local grid so acquisition density does not set the search
    extent. A missing valley between supplied apices is unresolved, not a midpoint.
    """
    x, y = validate_signal(x, y)
    if not x[0] < apex < x[-1]:
        raise ValueError("Вершина должна находиться внутри сигнала.")
    left = max(float(x[0]), apex - max_half_width)
    right = min(float(x[-1]), apex + max_half_width)
    dx = float(np.median(np.diff(x)))
    first = int(np.ceil((left-x[0])/dx))
    last = int(np.floor((right-x[0])/dx))
    gx = x[0] + np.arange(first, last+1)*dx
    if len(gx) < 5:
        raise ValueError("Недостаточно точек для определения формы пика.")
    gy = np.interp(gx, x, y)
    step = float(gx[1]-gx[0])
    window = max(5, int(round(smoothing_minutes/step)) | 1)
    window = min(window, len(gx) if len(gx) % 2 else len(gx)-1)
    smooth = savgol_filter(gy, window, 2)
    residual = gy-smooth
    noise = float(1.4826*np.median(np.abs(residual-np.median(residual))))
    center = int(np.argmin(abs(gx-apex)))
    height = float(smooth[center])
    if height <= max(8*noise, 0):
        return BoundaryProposal(apex, apex, apex, "noise", "noise", noise)
    reach = max(2, int(np.ceil(persistence_minutes/step)))
    lo, hi = max(0, center-reach), min(len(gx)-1, center+reach)
    if center == lo or center == hi or min(smooth[center]-smooth[lo], smooth[center]-smooth[hi]) <= 0:
        return BoundaryProposal(apex, apex, apex, "not_apex", "not_apex", noise)
    floor = max(3*noise, height*0.005)
    hold = max(2, int(np.ceil(persistence_minutes/step)))
    neighbors = sorted(float(v) for v in neighbors if np.isfinite(v) and abs(v-apex) > dx*2)

    def side(direction):
        candidates = [v for v in neighbors if direction*(v-apex) > 0]
        neighbor = min(candidates, key=lambda v: abs(v-apex)) if candidates else None
        limit = 0 if direction < 0 else len(gx)-1
        kind = "limit"
        if neighbor is not None and left <= neighbor <= right:
            n = int(np.argmin(abs(gx-neighbor)))
            lo, hi = sorted((center, n))
            valley = lo+int(np.argmin(smooth[lo:hi+1]))
            # A real interior minimum must rise above noise toward both apices.
            if lo < valley < hi and min(smooth[lo], smooth[hi])-smooth[valley] > 3*noise:
                limit, kind = valley, "valley"
            else:
                limit, kind = n, "unresolved"
        indices = np.arange(center, limit+direction, direction, dtype=int)
        for offset in range(1, max(1, len(indices)-hold+1)):
            segment = indices[offset:offset+hold]
            if len(segment) == hold and np.all(smooth[segment] <= floor):
                return float(gx[segment[0]]), "baseline"
        return float(gx[limit]), kind

    start, left_kind = side(-1)
    end, right_kind = side(1)
    return BoundaryProposal(start, end, float(apex), left_kind, right_kind, noise)


def boundary_issues(frame):
    """Geometry checks, not an accuracy claim or a fitted-component area check."""
    issues = []
    intervals = []
    for row in frame.to_dict("records"):
        code = str(row["code"])
        start, apex, end = (row.get(key, np.nan) for key in
                            ("integration_start_x", "found_rt", "integration_end_x"))
        if not np.isfinite([start, apex, end]).all():
            continue
        if not start <= apex <= end or end <= start:
            issues.append({"code": code, "reason": "apex_outside_bounds"})
        # Fitted components may legitimately have overlapping support intervals.
        if "fit" not in str(row.get("status", "")):
            intervals.append((float(start), float(end), code))
    intervals.sort()
    for i, (start, end, code) in enumerate(intervals):
        for other_start, other_end, other_code in intervals[i+1:]:
            if other_start >= end-1e-9:
                break
            issues.append({"code": code, "other_code": other_code, "reason": "overlapping_intervals"})
    return issues
