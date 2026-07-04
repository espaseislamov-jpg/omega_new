from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class BatchInput:
    sample_name: str
    dataframe: pd.DataFrame
    source_path: Path | None = None


@dataclass(frozen=True)
class OmegaResult:
    sample_name: str
    omega_report: float
    matched_targets: pd.DataFrame
    processed_chromatogram: pd.DataFrame
    confidence: dict
