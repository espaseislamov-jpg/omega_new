from __future__ import annotations

import json
import math
import os
import tempfile
import uuid
from copy import deepcopy
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from . import io, rt_profile


SCHEMA_VERSION = 1
STORE_FILE_NAME = "omega_instrument_profiles.json"
LEGACY_PROFILE_ID = "omega-v26-compatible"
PROFILE_FILE_SUFFIX = ".omega-profile.json"
LEGACY_CALCULATION_MODE = "legacy"
DIRECT_COMPONENTS_CALCULATION_MODE = "direct_components"
DEFAULT_OMEGA_COMPONENT_CODES = ("C20:5", "C22:6", "C22:5")
CALCULATION_MODES = frozenset({LEGACY_CALCULATION_MODE, DIRECT_COMPONENTS_CALCULATION_MODE})
LEGACY_JUDGE_MANUAL_SAMPLES = 411
LEGACY_JUDGE_ERROR_SAMPLES = 28
LEGACY_JUDGE_VALIDATION_BATCHES = 17
LEGACY_C22_HEIGHT_PRIORITY_THRESHOLD = 1.3635
EMERGENCY_JUDGE_TARGET_ABS_ERROR = 0.8


def _ordered_codes(reference_targets: pd.DataFrame | None = None) -> list[str]:
    if reference_targets is not None and not reference_targets.empty:
        ordered = reference_targets.sort_values("order_index")
        return ordered["code"].astype(str).tolist()
    return list(rt_profile.MANUAL_TABLE_RTS)


def legacy_profile(reference_targets: pd.DataFrame | None = None) -> dict[str, Any]:
    return {
        "id": LEGACY_PROFILE_ID,
        "name": "118 прибор",
        "custom_rt": False,
        "calculation_mode": LEGACY_CALCULATION_MODE,
        "omega_component_codes": list(DEFAULT_OMEGA_COMPONENT_CODES),
        "result_multiplier": 1.0,
        "judge_calibrated": True,
        "judge_emergency_enabled": False,
        "judge_target_abs_error": 0.5,
        "judge_manual_samples": LEGACY_JUDGE_MANUAL_SAMPLES,
        "judge_error_samples": LEGACY_JUDGE_ERROR_SAMPLES,
        "judge_validation_batches": LEGACY_JUDGE_VALIDATION_BATCHES,
        "judge_c22_height_priority_threshold": LEGACY_C22_HEIGHT_PRIORITY_THRESHOLD,
        "retention_times": {
            code: float(rt_profile.MANUAL_TABLE_RTS[code])
            for code in _ordered_codes(reference_targets)
            if code in rt_profile.MANUAL_TABLE_RTS
        },
    }


def new_profile(name: str, reference_targets: pd.DataFrame) -> dict[str, Any]:
    profile = legacy_profile(reference_targets)
    profile.update({
        "id": str(uuid.uuid4()),
        "name": str(name).strip(),
        "custom_rt": True,
        "judge_calibrated": False,
        "judge_emergency_enabled": False,
        "judge_target_abs_error": None,
        "judge_manual_samples": 0,
        "judge_error_samples": 0,
        "judge_validation_batches": 0,
        "judge_c22_height_priority_threshold": None,
    })
    return profile


def validate_profile(profile: dict[str, Any], required_codes: Iterable[str]) -> dict[str, Any]:
    if not isinstance(profile, dict):
        raise ValueError("Профиль должен быть объектом JSON.")
    out = deepcopy(profile)
    name = str(out.get("name", "")).strip()
    if not name:
        raise ValueError("Укажите название профиля.")
    multiplier = float(out.get("result_multiplier", 1.0))
    if not math.isfinite(multiplier) or not 0.01 <= multiplier <= 100.0:
        raise ValueError("Множитель должен быть положительным числом от 0.01 до 100.")

    calculation_mode = str(out.get("calculation_mode", LEGACY_CALCULATION_MODE)).strip()
    if calculation_mode not in CALCULATION_MODES:
        raise ValueError("Неизвестная схема расчёта Omega-3.")

    custom_rt = bool(out.get("custom_rt", True))
    judge_calibrated = bool(out.get("judge_calibrated", not custom_rt))
    judge_emergency_enabled = bool(out.get("judge_emergency_enabled", False))
    raw_target_abs_error = out.get("judge_target_abs_error")
    if raw_target_abs_error in (None, ""):
        judge_target_abs_error = None
    else:
        try:
            judge_target_abs_error = float(raw_target_abs_error)
        except (TypeError, ValueError):
            raise ValueError("Целевой порог экстренного судьи должен быть положительным числом.") from None
        if not math.isfinite(judge_target_abs_error) or judge_target_abs_error <= 0:
            raise ValueError("Целевой порог экстренного судьи должен быть положительным числом.")
    if judge_emergency_enabled and judge_target_abs_error is None:
        judge_target_abs_error = EMERGENCY_JUDGE_TARGET_ABS_ERROR
    try:
        judge_manual_samples = max(0, int(out.get("judge_manual_samples", 0)))
        judge_error_samples = max(0, int(out.get("judge_error_samples", 0)))
        judge_validation_batches = max(0, int(out.get("judge_validation_batches", 0)))
    except (TypeError, ValueError):
        raise ValueError("Сведения о калибровке судьи должны быть целыми неотрицательными числами.") from None
    if judge_error_samples > judge_manual_samples:
        raise ValueError("Число контрольных ошибок судьи не может превышать число ручных проб.")
    raw_height_threshold = out.get("judge_c22_height_priority_threshold")
    if raw_height_threshold in (None, ""):
        judge_c22_height_priority_threshold = None
    else:
        try:
            judge_c22_height_priority_threshold = float(raw_height_threshold)
        except (TypeError, ValueError):
            raise ValueError("Порог отношения высот C22 должен быть положительным числом.") from None
        if not math.isfinite(judge_c22_height_priority_threshold) or judge_c22_height_priority_threshold <= 0:
            raise ValueError("Порог отношения высот C22 должен быть положительным числом.")
    raw_rts = out.get("retention_times")
    if not isinstance(raw_rts, dict):
        raise ValueError("В профиле отсутствует таблица времён выхода.")
    codes = [str(code) for code in required_codes]
    raw_component_codes = out.get("omega_component_codes", DEFAULT_OMEGA_COMPONENT_CODES)
    if not isinstance(raw_component_codes, (list, tuple)):
        raise ValueError("Пики числителя Omega-3 должны быть списком.")
    omega_component_codes = list(dict.fromkeys(str(code).strip() for code in raw_component_codes))
    if not 2 <= len(omega_component_codes) <= 5:
        raise ValueError("В числителе Omega-3 должно быть от 2 до 5 разных пиков.")
    unknown_components = [code for code in omega_component_codes if code not in codes]
    if unknown_components:
        raise ValueError(f"Неизвестные пики числителя: {', '.join(unknown_components)}.")
    retention_times: dict[str, float] = {}
    for code in codes:
        try:
            value = float(raw_rts[code])
        except (KeyError, TypeError, ValueError):
            raise ValueError(f"Не задано корректное время выхода для {code}.") from None
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"Время выхода для {code} должно быть положительным числом.")
        retention_times[code] = value
    if custom_rt:
        values = [retention_times[code] for code in codes]
        if any(right <= left for left, right in zip(values, values[1:])):
            raise ValueError("Времена выхода должны строго возрастать в порядке пиков.")

    out.update({
        "id": str(out.get("id") or uuid.uuid4()),
        "name": name,
        "custom_rt": custom_rt,
        "calculation_mode": calculation_mode,
        "omega_component_codes": omega_component_codes,
        "result_multiplier": multiplier,
        "judge_calibrated": judge_calibrated,
        "judge_emergency_enabled": judge_emergency_enabled,
        "judge_target_abs_error": judge_target_abs_error,
        "judge_manual_samples": judge_manual_samples,
        "judge_error_samples": judge_error_samples,
        "judge_validation_batches": judge_validation_batches,
        "judge_c22_height_priority_threshold": judge_c22_height_priority_threshold,
        "retention_times": retention_times,
    })
    return out


def default_store(reference_targets: pd.DataFrame) -> dict[str, Any]:
    profile = legacy_profile(reference_targets)
    return {
        "schema_version": SCHEMA_VERSION,
        "active_profile_id": profile["id"],
        "profiles": [profile],
    }


def store_path() -> Path:
    return io.get_runtime_app_dir() / STORE_FILE_NAME


def _atomic_json_write(path: Path, payload: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    handle, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            json.dump(payload, stream, ensure_ascii=False, indent=2)
            stream.write("\n")
        os.replace(temp_name, path)
    except Exception:
        try:
            os.unlink(temp_name)
        except OSError:
            pass
        raise


def normalize_store(raw: Any, reference_targets: pd.DataFrame) -> dict[str, Any]:
    codes = _ordered_codes(reference_targets)
    legacy = legacy_profile(reference_targets)
    profiles: list[dict[str, Any]] = [legacy]
    seen = {LEGACY_PROFILE_ID}
    if isinstance(raw, dict):
        for item in raw.get("profiles", []):
            if not isinstance(item, dict) or str(item.get("id")) == LEGACY_PROFILE_ID:
                continue
            try:
                profile = validate_profile(item, codes)
            except (ValueError, TypeError):
                continue
            if profile["id"] in seen:
                profile["id"] = str(uuid.uuid4())
            profiles.append(profile)
            seen.add(profile["id"])
    active_id = str(raw.get("active_profile_id", LEGACY_PROFILE_ID)) if isinstance(raw, dict) else LEGACY_PROFILE_ID
    if active_id not in seen:
        active_id = LEGACY_PROFILE_ID
    return {"schema_version": SCHEMA_VERSION, "active_profile_id": active_id, "profiles": profiles}


def load_store(reference_targets: pd.DataFrame, path: Path | None = None) -> dict[str, Any]:
    path = Path(path) if path is not None else store_path()
    raw: Any = None
    if path.exists():
        try:
            with path.open("r", encoding="utf-8") as stream:
                raw = json.load(stream)
        except (OSError, json.JSONDecodeError):
            raw = None
    store = normalize_store(raw, reference_targets)
    if raw != store:
        _atomic_json_write(path, store)
    return store


def save_store(store: dict[str, Any], reference_targets: pd.DataFrame, path: Path | None = None) -> dict[str, Any]:
    normalized = normalize_store(store, reference_targets)
    _atomic_json_write(Path(path) if path is not None else store_path(), normalized)
    return normalized


def active_profile(store: dict[str, Any]) -> dict[str, Any]:
    active_id = str(store.get("active_profile_id", LEGACY_PROFILE_ID))
    for profile in store.get("profiles", []):
        if str(profile.get("id")) == active_id:
            return profile
    return store["profiles"][0]


def unique_profile_name(store: dict[str, Any], desired: str) -> str:
    base = str(desired).strip() or "Импортированный профиль"
    existing = {str(item.get("name", "")).casefold() for item in store.get("profiles", [])}
    if base.casefold() not in existing:
        return base
    number = 2
    while f"{base} ({number})".casefold() in existing:
        number += 1
    return f"{base} ({number})"


def export_profile(profile: dict[str, Any], path: Path, required_codes: Iterable[str]) -> None:
    clean = validate_profile(profile, required_codes)
    payload = {"schema_version": SCHEMA_VERSION, "profile": clean}
    _atomic_json_write(Path(path), payload)


def import_profile(path: Path, store: dict[str, Any], reference_targets: pd.DataFrame) -> dict[str, Any]:
    with Path(path).open("r", encoding="utf-8") as stream:
        raw = json.load(stream)
    candidate = raw.get("profile") if isinstance(raw, dict) and "profile" in raw else raw
    profile = validate_profile(candidate, _ordered_codes(reference_targets))
    profile["id"] = str(uuid.uuid4())
    profile["custom_rt"] = True
    profile["name"] = unique_profile_name(store, profile["name"])
    return profile


def apply_profile_to_targets(reference_targets: pd.DataFrame, profile: dict[str, Any]) -> pd.DataFrame:
    out = reference_targets.copy()
    custom_rt = bool(profile.get("custom_rt", False))
    if custom_rt:
        rts = profile["retention_times"]
        out["expected_rt"] = out["code"].map(rts).astype(float)
        out["rt_reliable"] = True
    out["instrument_profile_id"] = str(profile.get("id", LEGACY_PROFILE_ID))
    out["instrument_profile_name"] = str(profile.get("name", ""))
    out["instrument_profile_custom_rt"] = custom_rt
    out["instrument_profile_calculation_mode"] = str(
        profile.get("calculation_mode", LEGACY_CALCULATION_MODE)
    )
    component_codes = set(profile.get("omega_component_codes", DEFAULT_OMEGA_COMPONENT_CODES))
    out["instrument_profile_omega_component"] = out["code"].astype(str).isin(component_codes)
    out["result_multiplier"] = float(profile.get("result_multiplier", 1.0))
    out["instrument_profile_judge_calibrated"] = bool(profile.get("judge_calibrated", not custom_rt))
    out["instrument_profile_judge_emergency_enabled"] = bool(
        profile.get("judge_emergency_enabled", False)
    )
    target_abs_error = profile.get("judge_target_abs_error")
    out["instrument_profile_judge_target_abs_error"] = (
        float(target_abs_error) if target_abs_error not in (None, "") else np.nan
    )
    out["instrument_profile_judge_manual_samples"] = int(profile.get("judge_manual_samples", 0))
    out["instrument_profile_judge_error_samples"] = int(profile.get("judge_error_samples", 0))
    out["instrument_profile_judge_validation_batches"] = int(profile.get("judge_validation_batches", 0))
    threshold = profile.get("judge_c22_height_priority_threshold")
    out["instrument_profile_c22_height_priority_threshold"] = (
        float(threshold) if threshold not in (None, "") else np.nan
    )
    return out


def profile_from_targets(reference_targets: pd.DataFrame) -> dict[str, Any]:
    if reference_targets is None or reference_targets.empty:
        return legacy_profile(reference_targets)
    row = reference_targets.iloc[0]
    component_codes = (
        reference_targets.loc[
            reference_targets["instrument_profile_omega_component"].fillna(False).astype(bool),
            "code",
        ].astype(str).tolist()
        if "instrument_profile_omega_component" in reference_targets
        else list(DEFAULT_OMEGA_COMPONENT_CODES)
    )
    return {
        "id": str(row.get("instrument_profile_id", LEGACY_PROFILE_ID)),
        "name": str(row.get("instrument_profile_name", "118 прибор")),
        "custom_rt": bool(row.get("instrument_profile_custom_rt", False)),
        "calculation_mode": str(
            row.get("instrument_profile_calculation_mode", LEGACY_CALCULATION_MODE)
        ),
        "omega_component_codes": component_codes,
        "result_multiplier": float(row.get("result_multiplier", 1.0)),
        "judge_calibrated": bool(
            row.get(
                "instrument_profile_judge_calibrated",
                not bool(row.get("instrument_profile_custom_rt", False)),
            )
        ),
        "judge_emergency_enabled": bool(
            row.get("instrument_profile_judge_emergency_enabled", False)
        ),
        "judge_target_abs_error": (
            float(row.get("instrument_profile_judge_target_abs_error"))
            if pd.notna(row.get("instrument_profile_judge_target_abs_error"))
            else None
        ),
        "judge_manual_samples": int(row.get("instrument_profile_judge_manual_samples", 0)),
        "judge_error_samples": int(row.get("instrument_profile_judge_error_samples", 0)),
        "judge_validation_batches": int(row.get("instrument_profile_judge_validation_batches", 0)),
        "judge_c22_height_priority_threshold": (
            float(row.get("instrument_profile_c22_height_priority_threshold"))
            if pd.notna(row.get("instrument_profile_c22_height_priority_threshold"))
            else None
        ),
    }


def judge_calibration_label(profile: dict[str, Any]) -> str:
    if bool(profile.get("judge_emergency_enabled", False)):
        threshold = float(
            profile.get("judge_target_abs_error") or EMERGENCY_JUDGE_TARGET_ABS_ERROR
        )
        samples = max(0, int(profile.get("judge_manual_samples", 0)))
        suffix = f", {samples} проб" if samples else ""
        return f"Экстренный ±{threshold:.1f}{suffix}"
    if not bool(profile.get("judge_calibrated", False)):
        return "Не настроен"
    samples = max(0, int(profile.get("judge_manual_samples", 0)))
    errors = max(0, int(profile.get("judge_error_samples", 0)))
    if samples:
        return f"Готов: {samples}/{errors}"
    return "Готов"


def apply_result_multiplier(result: dict[str, Any], profile: dict[str, Any]) -> dict[str, Any]:
    out = dict(result)
    multiplier = float(profile.get("result_multiplier", 1.0))
    omega = dict(out.get("omega", {}))
    for key in ("omega3_trio", "omega3_trio_strict", "omega3_trio_corrected"):
        value = float(omega.get(key, np.nan))
        omega[f"{key}_unscaled"] = value
        if np.isfinite(value):
            omega[key] = value * multiplier
    out["omega"] = omega
    out["omega_report_unscaled"] = float(result.get("omega_report", np.nan))
    out["omega_report"] = omega.get("omega3_trio", np.nan)
    out["instrument_profile_id"] = str(profile.get("id", LEGACY_PROFILE_ID))
    out["instrument_profile_name"] = str(profile.get("name", ""))
    out["result_multiplier"] = multiplier
    return out
