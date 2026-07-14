from __future__ import annotations

import csv
import json
import re
import shutil
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from omega_batch_order import sort_batches_by_acquisition


PROJECT_DIR = Path(__file__).resolve().parents[1]
DEFAULT_REFERENCE_PATH = PROJECT_DIR / "reference_targets_reverted_c22fixed.json"

_SAMPLE_NAME_RE = re.compile(r"\bO\d+_[A-Za-z0-9._-]+\b")


def get_runtime_resource_dir() -> Path:
    if getattr(sys, "frozen", False):
        meipass = getattr(sys, "_MEIPASS", None)
        if meipass:
            return Path(meipass)
        return Path(sys.executable).resolve().parent
    return PROJECT_DIR


def get_runtime_app_dir() -> Path:
    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent
    return PROJECT_DIR


def ensure_runtime_file(file_name: str) -> Path:
    app_dir = get_runtime_app_dir()
    app_dir.mkdir(parents=True, exist_ok=True)
    target_path = app_dir / file_name
    if target_path.exists():
        return target_path

    source_path = get_runtime_resource_dir() / file_name
    if source_path.exists():
        try:
            shutil.copy2(source_path, target_path)
        except OSError:
            pass
    return target_path


def load_reference_targets(reference_path: Path = DEFAULT_REFERENCE_PATH) -> pd.DataFrame:
    reference_path = Path(reference_path)
    if not reference_path.exists() and reference_path.name == DEFAULT_REFERENCE_PATH.name:
        reference_path = ensure_runtime_file(reference_path.name)
    if not reference_path.exists():
        raise FileNotFoundError(f"Reference JSON not found: {reference_path}")

    with reference_path.open("r", encoding="utf-8") as f:
        raw_targets = json.load(f)

    if not isinstance(raw_targets, list):
        raise ValueError("Reference JSON must contain a list of targets.")

    records: list[dict[str, Any]] = []
    for idx, item in enumerate(raw_targets, start=1):
        if not isinstance(item, dict):
            continue
        component = item.get("component", item.get("code", f"target_{idx}"))
        code = item.get("code", component)
        records.append({
            "component": component,
            "code": code,
            "display_name": item.get("display_name", component),
            "order_index": item.get("order_index", idx),
            "expected_rt": item.get("expected_rt", item.get("target_rt")),
            "rt_reliable": bool(item.get("rt_reliable", False)),
            "historical_area": item.get("historical_area", 0.0),
            "historical_percent": item.get("historical_percent", 0.0),
            "notes": item.get("notes", ""),
        })

    df = pd.DataFrame(records)
    if df.empty:
        raise ValueError("Reference target list is empty.")

    fallback_order = pd.Series(np.arange(1, len(df) + 1), index=df.index, dtype=float)
    order_index = pd.to_numeric(df["order_index"], errors="coerce")
    df["order_index"] = order_index.where(order_index.notna(), fallback_order).astype(int)
    df["expected_rt"] = pd.to_numeric(df["expected_rt"], errors="coerce")
    df["historical_area"] = pd.to_numeric(df["historical_area"], errors="coerce").fillna(0.0)
    df["historical_percent"] = pd.to_numeric(df["historical_percent"], errors="coerce").fillna(0.0)
    df["rt_reliable"] = df["rt_reliable"].fillna(False).astype(bool)
    return df.sort_values("order_index").reset_index(drop=True)


def extract_sample_name_from_header(file_path: Path) -> str:
    try:
        with Path(file_path).open("r", encoding="utf-8", errors="ignore") as f:
            header_lines = [next(f, "") for _ in range(3)]
    except OSError:
        return Path(file_path).stem

    match = _SAMPLE_NAME_RE.search(" ".join(header_lines))
    return match.group(0) if match else Path(file_path).stem


def extract_sample_name_from_text(text: str, fallback: str) -> str:
    match = _SAMPLE_NAME_RE.search(text or "")
    return match.group(0) if match else fallback


def finalize_chromatogram_dataframe(df: pd.DataFrame, cutoff_minutes: float = 4.0) -> pd.DataFrame:
    df = df.copy()
    df["x"] = pd.to_numeric(df["x"], errors="coerce")
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    df = df.dropna(subset=["x", "y"]).reset_index(drop=True)

    unique_x = np.sort(df["x"].unique())
    if len(unique_x) < 2:
        raise ValueError("Not enough unique x values to estimate chromatogram step.")

    step = float(np.median(np.diff(unique_x)))
    corrected_x: list[float] = []
    x_values = df["x"].to_numpy(dtype=float)
    i = 0
    n = len(df)
    while i < n:
        current_x = x_values[i]
        j = i + 1
        while j < n and x_values[j] == current_x:
            j += 1
        group_size = j - i
        offsets = np.linspace(0.0, step * (group_size - 1) / max(group_size, 1), group_size)
        corrected_x.extend((current_x + offsets).tolist())
        i = j

    df["x_corrected"] = corrected_x
    return df[df["x_corrected"] >= cutoff_minutes].reset_index(drop=True)


def load_chromatogram(file_path: Path, cutoff_minutes: float = 4.0) -> pd.DataFrame:
    df = pd.read_csv(file_path, skiprows=3, header=None, names=["x", "y"])
    return finalize_chromatogram_dataframe(df, cutoff_minutes=cutoff_minutes)


def is_chromtab_file(file_path: Path) -> bool:
    try:
        with Path(file_path).open("r", encoding="utf-8", errors="ignore") as f:
            first_line = next(f, "").strip()
            next(f, "")
            third_line = next(f, "").strip()
    except OSError:
        return False
    return first_line.startswith('"Path","File","Date Acquired"') and third_line.startswith('"Signal: ')


def load_chromtab_batches(file_path: Path, cutoff_minutes: float = 4.0) -> list[dict[str, Any]]:
    batches: list[dict[str, Any]] = []
    current_meta: dict[str, str] | None = None
    current_rows: list[tuple[float, float]] = []

    def flush_current() -> None:
        nonlocal current_meta, current_rows
        if current_meta is None or not current_rows:
            current_rows = []
            return

        raw_df = pd.DataFrame(current_rows, columns=["x", "y"])
        batch_df = finalize_chromatogram_dataframe(raw_df, cutoff_minutes=cutoff_minutes)
        file_name = current_meta.get("file_name", "")
        sample_name = extract_sample_name_from_text(
            " ".join(filter(None, [current_meta.get("signal_name", ""), file_name])),
            fallback=Path(file_name).stem if file_name else f"batch_{len(batches) + 1}",
        )
        batches.append({
            "file_name": file_name,
            "signal_name": current_meta.get("signal_name", ""),
            "acquired_at": current_meta.get("acquired_at", ""),
            "source_path": current_meta.get("source_path", ""),
            "sample_name": sample_name,
            "dataframe": batch_df,
        })
        current_meta = None
        current_rows = []

    with Path(file_path).open("r", encoding="utf-8", errors="ignore") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith('"Path","File","Date Acquired"'):
                flush_current()
                current_meta = {"source_path": "", "file_name": "", "acquired_at": "", "signal_name": ""}
                continue

            if current_meta is not None and not current_meta["file_name"] and line.startswith('"'):
                parsed = next(csv.reader([line]))
                current_meta["source_path"] = parsed[0] if len(parsed) > 0 else ""
                current_meta["file_name"] = parsed[1] if len(parsed) > 1 else ""
                current_meta["acquired_at"] = parsed[2] if len(parsed) > 2 else ""
                continue

            if current_meta is not None and line.startswith('"Signal: '):
                current_meta["signal_name"] = line.strip('"')
                continue

            if current_meta is None:
                continue

            parsed = next(csv.reader([line]), [])
            if len(parsed) < 2:
                continue
            try:
                current_rows.append((float(parsed[0]), float(parsed[1])))
            except ValueError:
                continue

    flush_current()
    if not batches:
        raise ValueError(f"No chromatogram batches parsed from {file_path}")
    return sort_batches_by_acquisition(batches)


def load_batches(file_path: Path, cutoff_minutes: float = 4.0) -> list[dict[str, Any]]:
    file_path = Path(file_path)
    if is_chromtab_file(file_path):
        return load_chromtab_batches(file_path, cutoff_minutes=cutoff_minutes)

    return [{
        "file_name": file_path.name,
        "signal_name": "",
        "acquired_at": "",
        "source_path": str(file_path),
        "sample_name": extract_sample_name_from_header(file_path),
        "dataframe": load_chromatogram(file_path, cutoff_minutes=cutoff_minutes),
    }]
