from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import io, pipeline


DEFAULT_REFERENCE_PATH = io.DEFAULT_REFERENCE_PATH


def load_reference_targets(reference_path: Path = DEFAULT_REFERENCE_PATH) -> pd.DataFrame:
    return io.load_reference_targets(reference_path)


def load_batches(file_path: Path, cutoff_minutes: float = 4.0) -> list[dict]:
    return io.load_batches(file_path, cutoff_minutes=cutoff_minutes)


def process_batch(dataframe: pd.DataFrame, reference_targets: pd.DataFrame) -> dict:
    return pipeline.process_batch(dataframe, reference_targets)


def process_file(
    file_path: Path,
    reference_path: Path = DEFAULT_REFERENCE_PATH,
    cutoff_minutes: float = 4.0,
) -> list[dict]:
    return pipeline.process_file(
        file_path=file_path,
        reference_path=reference_path,
        cutoff_minutes=cutoff_minutes,
    )
