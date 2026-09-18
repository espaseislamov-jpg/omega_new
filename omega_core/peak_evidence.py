"""Measurement evidence is independent of detector IDs and the Omega value."""
import numpy as np
import pandas as pd
from .integration import boundary_issues

CORE_DENOMINATOR_CODES = ("C16:0", "C18:2N6C", "C18:1N9C", "C18:0")


def measurement_states(frame):
    states = []
    for row in frame.to_dict("records"):
        status = str(row.get("status", ""))
        values = pd.to_numeric(pd.Series([row.get(k, np.nan) for k in
            ("area", "found_rt", "integration_start_x", "integration_end_x")]), errors="coerce").to_numpy()
        area, apex, start, end = values
        if not np.isfinite(values).all() or area <= 0 or not start < apex < end:
            states.append("missing" if not np.isfinite(area) or area <= 0 else "invalid")
        elif "manual" in status:
            states.append("manual")
        elif "unresolved" in status or "not_found" in status or "estimated" in status:
            states.append("unresolved")
        elif "fit" in status or "recovered_c20_missing_epa_component" in status:
            states.append("model")
        else:
            states.append("observed")
    return pd.Series(states, index=frame.index, dtype=str)


def uses_configured_base_judge(frame):
    from .rt_profile import uses_custom_instrument_profile
    return (
        not frame.empty and not uses_custom_instrument_profile(frame)
        and 'instrument_profile_judge_calibrated' in frame
        and frame['instrument_profile_judge_calibrated'].fillna(False).astype(bool).all()
    )


def assess_final_evidence(frame, confidence):
    out = _assess_final_evidence(frame, confidence)
    recovered = frame.status.eq('recovered_epa_local_review')
    if recovered.any():
        from copy import deepcopy
        out = deepcopy(out)
        codes = frame.loc[recovered, 'code'].astype(str).tolist()
        reason = 'ЭПК найдена дополнительным поиском: проверьте вершину и границы интегрирования'
        out['reasons'] = [reason] + [r for r in out.get('reasons', []) if r != reason]
        risk = out.setdefault('high_error_risk', {})
        risk['score'] = max(85, risk.get('score', 0))
        risk['peak_codes'] = sorted(set(risk.get('peak_codes', [])) | set(codes))
        risk['reason_codes'] = sorted(set(risk.get('reason_codes', [])) | {'epa_local_recovery'})
        out['score'] = min(out.get('score', 100), 54)
        if risk['score'] < 95 and not str(out.get('button_text', '')).startswith('СТОП'):
            out['button_text'] = 'ПРОВЕРИТЬ — ' + ', '.join(risk['peak_codes'])
        out['label'] = out['level'] = out['button_text']
    return out


def _assess_final_evidence(frame, confidence):
    """Read-only final judge; never participates in baseline/candidate selection."""
    if ('geometry_method' in frame and frame.geometry_method.eq('derivative_valley_v1').any()
            and not uses_configured_base_judge(frame)):
        return assess_geometry_evidence(frame)
    from copy import deepcopy
    out = deepcopy(confidence)
    states = measurement_states(frame)
    confirmed = states.isin(["observed", "manual", "model"])
    out["evidence"] = {
        "confirmed_count": int(confirmed.sum()), "total_count": len(frame),
        "states": dict(zip(frame.code.astype(str), states)),
        "missing_codes": frame.loc[~confirmed, "code"].astype(str).tolist(),
    }
    # A manual/local integral need not have a native detector ID.
    old_missing = int(frame.matched_peak_id.isna().sum())
    actual_missing = int((~confirmed).sum())
    old_score = float(out.get("geometry_score", out.get("score", 0)))
    out["selection_geometry_score"] = old_score
    out["geometry_score"] = float(np.clip(old_score + min(18, 4*old_missing) - min(18, 4*actual_missing), 0, 100))
    out["reasons"] = [r for r in out.get("reasons", []) if "Некоторые пики не удалось определить" not in r]
    if actual_missing:
        out["reasons"].insert(0, "Не подтверждены: " + ", ".join(out["evidence"]["missing_codes"]))
    measured = frame.loc[states.isin(["observed", "manual"])].copy()
    overlaps = [i for i in boundary_issues(measured) if i["reason"] == "overlapping_intervals"]
    model_codes = frame.loc[states == "model", "code"].astype(str).tolist()
    warnings, affected = [], set()
    if overlaps:
        affected.update(i["code"] for i in overlaps)
        affected.update(i["other_code"] for i in overlaps)
        warnings.append("Перекрываются интервалы измеренных пиков: " + ", ".join(sorted(affected)))
    if model_codes:
        affected.update(model_codes)
        warnings.append("Площади получены моделью: " + ", ".join(model_codes))
    if warnings:
        risk = out.setdefault("high_error_risk", {})
        risk["score"] = max(85, risk.get("score", 0))
        risk["peak_codes"] = sorted(affected | set(risk.get("peak_codes", [])))
        out["reasons"] = warnings + out["reasons"]
        out["button_text"] = "ПРОВЕРИТЬ — " + ", ".join(sorted(affected))
    risk_score = out.get("high_error_risk", {}).get("score", 0)
    geometry = out["geometry_score"]
    out["score"] = min(geometry, 35 if risk_score >= 95 else 54 if risk_score >= 85 else 100)
    if risk_score >= 95:
        out["button_text"] = "СТОП — ручная проверка пиков"
    elif risk_score >= 85:
        out["button_text"] = "ПРОВЕРИТЬ — " + ", ".join(out.get("high_error_risk", {}).get("peak_codes", []))
    elif geometry < 55:
        out["button_text"] = "СТОП — ручная проверка пиков"
    elif geometry < 85:
        out["button_text"] = "ПРОВЕРИТЬ — геометрию пиков"
    else:
        out["button_text"] = "ПРОВЕРКИ ПРОЙДЕНЫ"
    out["metrics"] = [m for m in out.get("metrics", []) if not m.startswith("Matched peaks:")]
    out["metrics"].append(f"Подтверждённые измерения: {int(confirmed.sum())}/{len(frame)}")
    return out


def assess_geometry_evidence(frame):
    """Evaluate physical support, not empirical area ratios or the final percentage."""
    states = measurement_states(frame)
    measured = states.isin(['observed', 'manual'])
    reasons, affected = [], set(frame.loc[~measured, 'code'].astype(str))
    if affected:
        reasons.append('Не подтверждены измерением: ' + ', '.join(sorted(affected)))
    for _, row in frame.loc[states == 'observed'].iterrows():
        noise = float(row.get('noise', np.nan))
        prominence = float(row.get('prominence', np.nan))
        if noise > 0 and prominence/noise < 10:
            affected.add(str(row.code))
            reasons.append(f'{row.code}: слабое превышение над локальным шумом')
        if float(row.get('valley_fraction', 0)) > .5:
            affected.add(str(row.code))
            reasons.append(f'{row.code}: высокая впадина между пиками; площади чувствительны к разделению')
        if float(row.get('boundary_uncertainty', 0)) > .002:
            affected.add(str(row.code))
            reasons.append(f'{row.code}: положение впадины меняется при изменении масштаба поиска')
    issues = boundary_issues(frame.loc[measured])
    for issue in issues:
        affected.add(issue['code'])
        if 'other_code' in issue:
            affected.add(issue['other_code'])
        reasons.append(f"{issue['code']}: {issue['reason']}")
    ids = frame.loc[measured, 'matched_peak_id'].dropna()
    duplicate = bool(ids.duplicated().any())
    if duplicate:
        reasons.append('Один измеренный пик назначен нескольким кислотам')
        affected.update(frame.loc[frame.matched_peak_id.isin(ids[ids.duplicated()]), 'code'])
    stop = bool(issues or duplicate or states.isin(['invalid','unresolved']).any())
    score = 35 if stop else 70 if affected else 100
    label = 'СТОП — проверка пиков' if stop else ('ПРОВЕРИТЬ — ' + ', '.join(sorted(affected))) if affected else 'ПРОВЕРКИ ПРОЙДЕНЫ'
    return dict(score=score, geometry_score=score, level=label, label=label, button_text=label,
                reasons=reasons, metrics=[f'Измеренные пики: {int(measured.sum())}/{len(frame)}',
                                        'Оценка геометрии; не вероятность аналитической точности'],
                high_error_risk=dict(score=100 if stop else 85 if affected else 0,
                                     peak_codes=sorted(affected), reason_codes=[]),
                evidence=dict(confirmed_count=int(measured.sum()), total_count=len(frame),
                              states=dict(zip(frame.code.astype(str), states)),
                              missing_codes=frame.loc[~measured, 'code'].astype(str).tolist()))
