"""Reviewable manual intervals, undo history, and source-bound JSON persistence."""
from copy import deepcopy
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.signal import find_peaks

from .integration import integrate_interval, interval_points


def manual_integral(batch, start, end):
    df = batch["processed_df"]
    x = df["x_corrected" if "x_corrected" in df else "x"].to_numpy()
    y = df["y_corrected"].to_numpy()
    integral = integrate_interval(x, y, start, end)
    if integral.area > 0 and not start < integral.apex < end:
        # An adjacent peak's tail can be higher at the selected edge than the
        # small peak the operator is integrating. Locate an actual interior
        # maximum without changing the chosen interval, signal or area.
        local_x, local_y = interval_points(x, y, start, end)
        peaks, properties = find_peaks(local_y, prominence=0)
        positive = local_y[peaks] > 0
        peaks = peaks[positive]
        if len(peaks):
            apex = peaks[np.argmax(properties['prominences'][positive])]
            integral = replace(integral, apex=float(local_x[apex]))
    if integral.area <= 0 or not start < integral.apex < end:
        raise ValueError("Внутри выбранных границ нет положительной локальной вершины. Уточните интервал на графике.")
    return integral


def _push_undo(batch):
    undo = batch.setdefault('manual_undo', [])
    undo.append((batch['matched_targets_df'].copy(deep=True), deepcopy(batch.get('manual_overrides', {})),
                 deepcopy(batch.get('manual_extra_peaks', []))))
    del undo[:-100]


def add_unassigned_peak(batch, start, end, name='Новый пик'):
    integral = manual_integral(batch, start, end)
    _push_undo(batch)
    peak = dict(id=hashlib.sha256(f'{start}:{end}:{name}'.encode()).hexdigest()[:16],
                name=str(name).strip() or 'Новый пик', start=float(start), end=float(end),
                apex=integral.apex, area=integral.area)
    if any(p['id']==peak['id'] for p in batch.get('manual_extra_peaks', [])):
        batch['manual_undo'].pop()
        raise ValueError('Этот пик уже добавлен.')
    batch.setdefault('manual_extra_peaks', []).append(peak)
    return peak


def apply_edit(batch, code, start, end, source="manual"):
    frame = batch["matched_targets_df"].copy(deep=True)
    indices = frame.index[frame["code"] == code]
    if len(indices) != 1:
        raise ValueError(f"Не найден единственный пик {code}.")
    integral = manual_integral(batch, start, end)
    index = indices[0]
    old = {key: frame.at[index, key] for key in
           ("integration_start_x", "integration_end_x", "found_rt", "area", "status")}
    _push_undo(batch)
    for key, value in {"integration_start_x": start, "integration_end_x": end,
                       "found_rt": integral.apex, "area": integral.area,
                       "matched_peak_id": np.nan, "match_score": np.nan,
                       "status": "manual_" + source}.items():
        frame.at[index, key] = value
    total = frame["area"].fillna(0).sum()
    frame["percent_area"] = 100*frame["area"]/total if total > 0 else np.nan
    batch["matched_targets_df"] = frame
    if source == 'insert':
        batch['manual_extra_peaks'] = [p for p in batch.get('manual_extra_peaks', [])
                                       if abs(p['start']-start)>1e-7 or abs(p['end']-end)>1e-7]
    batch.setdefault("manual_overrides", {})[str(code)] = {"start": float(start), "end": float(end), "source": source}
    batch.setdefault("manual_history", []).append({
        "time": datetime.now(timezone.utc).isoformat(), "action": "integrate",
        "code": code, "source": source, "before": old,
        "after": {"start": start, "end": end, "apex": integral.apex, "area": integral.area},
    })
    return integral


def undo_edit(batch):
    history = batch.get("manual_undo", [])
    if not history:
        return False
    previous = history.pop()
    frame, overrides = previous[:2]
    batch['manual_extra_peaks'] = previous[2] if len(previous)>2 else []
    batch["matched_targets_df"] = frame
    batch["manual_overrides"] = overrides
    batch.setdefault("manual_history", []).append({"time": datetime.now(timezone.utc).isoformat(), "action": "undo"})
    return True


def file_digest(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def _profile_digest(profile):
    return hashlib.sha256(json.dumps(profile, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def export_edits(path, source_path, profile, batches):
    payload = {"schema": 1, "source_sha256": file_digest(source_path),
               "profile_sha256": _profile_digest(profile), "batches": []}
    for index, batch in enumerate(batches):
        if batch.get("manual_overrides") or batch.get('manual_extra_peaks'):
            payload["batches"].append({"index": index, "sample": batch["sample_name"],
                                       "baseline": batch.get("baseline_mode"),
                                       "overrides": batch.get("manual_overrides", {}),
                                       'extra_peaks': batch.get('manual_extra_peaks', [])})
    Path(path).write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def import_edits(path, source_path, profile, batches):
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    if payload.get("schema") != 1 or payload.get("source_sha256") != file_digest(source_path):
        raise ValueError("Правки относятся к другому CSV или неизвестному формату.")
    if payload.get("profile_sha256") != _profile_digest(profile):
        raise ValueError("Правки рассчитаны для другого профиля прибора.")
    # Validate and calculate all edits in independent copies before committing any.
    staged = {}
    for entry in payload["batches"]:
        index = entry["index"]
        if type(index) is not int or not 0 <= index < len(batches) or index in staged:
            raise ValueError("Некорректный или повторный номер пробы.")
        original = batches[index]
        if entry["sample"] != original["sample_name"] or entry["baseline"] != original.get("baseline_mode"):
            raise ValueError("Проба или базовая линия не соответствует сохранённой правке.")
        copy = dict(original)
        for key in ("manual_undo", "manual_overrides", "manual_history", 'manual_extra_peaks'):
            copy[key] = deepcopy(original.get(key, {} if key == "manual_overrides" else []))
        for code, values in entry["overrides"].items():
            apply_edit(copy, code, float(values["start"]), float(values["end"]), values.get("source", "manual"))
        for peak in entry.get('extra_peaks', []):
            if not any(p['id']==peak['id'] for p in copy.get('manual_extra_peaks', [])):
                add_unassigned_peak(copy, float(peak['start']), float(peak['end']), peak['name'])
        staged[index] = copy
    for index, value in staged.items():
        batches[index].update(value)
    return len(staged)
