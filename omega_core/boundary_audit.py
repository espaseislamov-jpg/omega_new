"""Small, JSON-friendly audit trail for integration changes; no numeric decisions."""
import math

FIELDS = ("found_rt", "integration_start_x", "integration_end_x", "area", "status")
ATTR = "boundary_history"


def snapshot(frame):
    return {str(row["code"]): {field: row.get(field) for field in FIELDS}
            for row in frame.to_dict("records")}


def _clean(value):
    if value is None or isinstance(value, str):
        return value
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return str(value)


def record_changes(before, after, stage, history):
    current = snapshot(after)
    for code in sorted(set(before) | set(current)):
        old = {key: _clean(value) for key, value in before.get(code, {}).items()}
        new = {key: _clean(value) for key, value in current.get(code, {}).items()}
        if old != new:
            history.append({"stage": stage, "code": code, "before": old, "after": new})
    return current
