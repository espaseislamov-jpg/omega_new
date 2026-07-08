from __future__ import annotations

import os

import numpy as np
import pandas as pd

# Manual workbook (ИНДЕКС РАСЧ.xlsx) retention-time landmarks.  These are the
# effective tabular RTs used by the spreadsheet algorithm, not the older
# reference JSON values where several C20/C22 components are intentionally
# duplicated or absent.
MANUAL_TABLE_RTS: dict[str, float] = {
    "C16:1N7": 6.594,
    "C16:0": 6.708,
    "C18:3N6": 7.497,
    "C18:2N6C": 7.582,
    "C18:1N9C": 7.615,
    "C18:3N3": 7.643,
    "C18:0": 7.743,
    "C20:4N6": 8.373,
    "C20:5": 8.402,
    "C20:3N8": 8.460,
    "C22:6": 9.239,
    "C22:5": 9.272,
    "C22:4": 9.301,
    "C24:1N9": 10.291,
    "C24:0": 10.383,
}

# Across the old+new batches these two targets had the most stable coefficient
# behaviour and are separated enough to avoid relying on the problematic C20:5
# or C22 cluster itself.
ANCHOR_CODES: tuple[str, ...] = ("C18:3N3", "C20:3N8")
ANCHOR_COEFFICIENT_FALLBACK = 0.99952
ENABLE_MANUAL_RT_PROFILE_TARGETING = os.environ.get("OMEGA_MANUAL_RT_PROFILE_TARGETING", "0").strip() == "1"


def estimate_anchor_coefficient(matched_targets: pd.DataFrame) -> float:
    """Return table_rt / observed_rt for the current chromatogram.

    The coefficient is a multiplicative RT stretch/compression estimate.  It is
    deliberately based only on stable anchors; if anchors are missing or clearly
    invalid, fall back to the corpus median instead of the older additive shift.
    """
    if matched_targets is None or matched_targets.empty or "code" not in matched_targets:
        return float(ANCHOR_COEFFICIENT_FALLBACK)

    values: list[float] = []
    for code in ANCHOR_CODES:
        row = matched_targets[matched_targets["code"] == code]
        if row.empty:
            continue
        observed = pd.to_numeric(row["found_rt"], errors="coerce").iloc[0]
        table_rt = MANUAL_TABLE_RTS.get(code)
        if table_rt is None or not np.isfinite(observed) or observed <= 0:
            continue
        coef = float(table_rt) / float(observed)
        if 0.985 <= coef <= 1.015:
            values.append(coef)

    if not values:
        return float(ANCHOR_COEFFICIENT_FALLBACK)
    return float(np.median(values))


def expected_rt(code: str, anchor_coefficient: float) -> float:
    table_rt = MANUAL_TABLE_RTS[code]
    coefficient = float(anchor_coefficient)
    if not np.isfinite(coefficient) or coefficient <= 0:
        coefficient = ANCHOR_COEFFICIENT_FALLBACK
    return float(table_rt) / coefficient


def expected_rts(codes: list[str] | tuple[str, ...], anchor_coefficient: float) -> list[float]:
    return [expected_rt(code, anchor_coefficient) for code in codes]


def choose_expected_rts(codes: list[str] | tuple[str, ...], anchor_coefficient: float, fallback_rts) -> list[float]:
    """Return manual-profile targets only when the guarded experiment is enabled."""
    if ENABLE_MANUAL_RT_PROFILE_TARGETING:
        return expected_rts(codes, anchor_coefficient)
    return [float(value) for value in fallback_rts]


def annotate_rt_profile(matched_targets: pd.DataFrame) -> pd.DataFrame:
    """Attach manual-table RT coefficients for diagnostics/regression reports."""
    out = matched_targets.copy()
    if out.empty or "code" not in out:
        return out
    anchor = estimate_anchor_coefficient(out)
    out["rt_profile_anchor_coefficient"] = anchor
    out["manual_table_rt"] = out["code"].map(MANUAL_TABLE_RTS)
    found = pd.to_numeric(out.get("found_rt"), errors="coerce")
    table = pd.to_numeric(out["manual_table_rt"], errors="coerce")
    out["rt_profile_coefficient"] = table / found
    out.loc[~np.isfinite(out["rt_profile_coefficient"]), "rt_profile_coefficient"] = np.nan
    out["rt_profile_expected_rt"] = table / anchor if np.isfinite(anchor) and anchor > 0 else np.nan
    out["rt_profile_delta_to_anchor"] = out["rt_profile_coefficient"] - anchor
    return out
