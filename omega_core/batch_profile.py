from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

PROFILE_TARGET_CODES = (
    "C18:2N6C",
    "C18:1N9C",
    "C18:3N3",
    "C20:4N6",
    "C20:5",
    "C20:3N8",
    "C22:6",
    "C22:5",
    "C22:4",
)


def _safe_float(value, default=np.nan) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float(default)
    return out if np.isfinite(out) else float(default)


def _robust_scale(values: np.ndarray, floor: float) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size <= 1:
        return float(floor)
    median = float(np.median(arr))
    mad = float(np.median(np.abs(arr - median)))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale < floor:
        scale = max(float(np.std(arr)), float(floor))
    return float(max(scale, floor))


def _summary(values: Iterable[float], scale_floor: float) -> dict:
    arr = np.asarray(list(values), dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return {"n": 0, "median": None, "mean": None, "std": None, "q25": None, "q75": None, "robust_scale": None}
    return {
        "n": int(arr.size),
        "median": float(np.median(arr)),
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=0)),
        "q25": float(np.quantile(arr, 0.25)),
        "q75": float(np.quantile(arr, 0.75)),
        "robust_scale": _robust_scale(arr, scale_floor),
    }


def build_profile(
    peak_diagnostics: pd.DataFrame,
    target_codes: Iterable[str] = PROFILE_TARGET_CODES,
    rt_floor: float = 0.0015,
    width_floor: float = 0.0020,
) -> dict:
    """Build robust RT/boundary profile from per-sample peak diagnostics.

    The profile is intentionally read-only: it learns stable medians and robust
    scales that can later be used as priors by a global integrator.
    """
    if peak_diagnostics is None or peak_diagnostics.empty:
        return {"targets": {}, "global": {"n_rows": 0}}

    frame = peak_diagnostics.copy()
    for column in ["found_rt", "rt_error", "width", "left_width", "right_width", "area", "asymmetry"]:
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    targets: dict[str, dict] = {}
    for code in target_codes:
        target = frame[frame["target_code"] == code].copy()
        if target.empty:
            continue
        rt_summary = _summary(target["found_rt"].dropna(), rt_floor)
        width_summary = _summary(target["width"].dropna(), width_floor)
        left_summary = _summary(target["left_width"].dropna(), width_floor)
        right_summary = _summary(target["right_width"].dropna(), width_floor)
        rt_center = rt_summary.get("median")
        rt_scale = rt_summary.get("robust_scale") or rt_floor
        width_center = width_summary.get("median")
        width_scale = width_summary.get("robust_scale") or width_floor
        rt_values = target["found_rt"].to_numpy(dtype=float)
        width_values = target["width"].to_numpy(dtype=float)
        rt_z = (rt_values - float(rt_center)) / float(rt_scale) if rt_center is not None else np.asarray([])
        width_z = (width_values - float(width_center)) / float(width_scale) if width_center is not None else np.asarray([])
        rt_z = rt_z[np.isfinite(rt_z)]
        width_z = width_z[np.isfinite(width_z)]
        targets[code] = {
            "rt": rt_summary,
            "width": width_summary,
            "left_width": left_summary,
            "right_width": right_summary,
            "chi2_like_rt_mean": float(np.mean(np.square(rt_z))) if rt_z.size else None,
            "chi2_like_width_mean": float(np.mean(np.square(width_z))) if width_z.size else None,
            "status_top": target["status"].fillna("").astype(str).value_counts().head(5).to_dict() if "status" in target else {},
        }
    return {
        "schema_version": 1,
        "target_codes": list(target_codes),
        "global": {
            "n_rows": int(len(frame)),
            "n_samples": int(frame[["date", "sample_name"]].drop_duplicates().shape[0]) if {"date", "sample_name"}.issubset(frame.columns) else None,
        },
        "targets": targets,
    }



def blend_profiles(previous: dict | None, current: dict, alpha: float = 0.35) -> dict:
    """Blend a previous persisted profile with a freshly measured batch profile.

    `alpha` is the current-batch weight. The function only blends numeric summary
    values that exist in both profiles; status counters and chi-square-like values
    are kept from the current run because they describe current convergence.
    """
    if not previous:
        return current
    alpha = float(min(1.0, max(0.0, alpha)))
    beta = 1.0 - alpha
    blended = json.loads(json.dumps(current))
    previous_targets = previous.get("targets", {}) if isinstance(previous, dict) else {}
    for code, current_target in current.get("targets", {}).items():
        previous_target = previous_targets.get(code)
        if not isinstance(previous_target, dict):
            continue
        out_target = blended["targets"].get(code, {})
        for section in ["rt", "width", "left_width", "right_width"]:
            current_section = current_target.get(section, {})
            previous_section = previous_target.get(section, {})
            out_section = out_target.get(section, {})
            for key in ["median", "mean", "std", "q25", "q75", "robust_scale"]:
                c_val = current_section.get(key)
                p_val = previous_section.get(key)
                if c_val is None or p_val is None:
                    continue
                try:
                    c_float = float(c_val)
                    p_float = float(p_val)
                except (TypeError, ValueError):
                    continue
                if np.isfinite(c_float) and np.isfinite(p_float):
                    out_section[key] = float(beta * p_float + alpha * c_float)
            out_section["n"] = int(current_section.get("n") or 0)
            out_target[section] = out_section
        blended["targets"][code] = out_target
    blended["blending"] = {
        "alpha_current": alpha,
        "previous_schema_version": previous.get("schema_version") if isinstance(previous, dict) else None,
    }
    return blended

def profile_to_frame(profile: dict) -> pd.DataFrame:
    rows = []
    for code, item in profile.get("targets", {}).items():
        rows.append({
            "target_code": code,
            "n": item.get("rt", {}).get("n"),
            "rt_median": item.get("rt", {}).get("median"),
            "rt_robust_scale": item.get("rt", {}).get("robust_scale"),
            "width_median": item.get("width", {}).get("median"),
            "left_width_median": item.get("left_width", {}).get("median"),
            "right_width_median": item.get("right_width", {}).get("median"),
            "chi2_like_rt_mean": item.get("chi2_like_rt_mean"),
            "chi2_like_width_mean": item.get("chi2_like_width_mean"),
        })
    return pd.DataFrame(rows)


def save_profile(profile: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(profile, ensure_ascii=False, indent=2), encoding="utf-8")
