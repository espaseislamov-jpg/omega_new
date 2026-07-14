from __future__ import annotations

from datetime import datetime
import re
from typing import Any, Iterable


_ENGLISH_MONTHS = {
    "jan": 1,
    "feb": 2,
    "mar": 3,
    "apr": 4,
    "may": 5,
    "jun": 6,
    "jul": 7,
    "aug": 8,
    "sep": 9,
    "oct": 10,
    "nov": 11,
    "dec": 12,
}


def sample_number_sort_key(name: str) -> tuple[int, int | str, str]:
    text = str(name or "")
    match = re.search(r"\bO(\d+)\b", text) or re.search(r"\bO(\d+)", text)
    if match:
        return (0, int(match.group(1)), text)
    return (1, text, text)


def parse_acquired_at(value: Any) -> datetime | None:
    """Parse Agilent CHROMTAB's locale-independent Date Acquired value."""
    text = " ".join(str(value or "").strip().split())
    if not text:
        return None

    english = re.fullmatch(
        r"(\d{1,2})\s+([A-Za-z]{3,9})\s+(\d{4})\s+(\d{1,2}):(\d{2})(?::(\d{2}))?",
        text,
    )
    if english:
        day, month_text, year, hour, minute, second = english.groups()
        month = _ENGLISH_MONTHS.get(month_text[:3].lower())
        if month is not None:
            try:
                return datetime(
                    int(year),
                    month,
                    int(day),
                    int(hour),
                    int(minute),
                    int(second or 0),
                )
            except ValueError:
                return None

    for format_text in ("%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M", "%d.%m.%Y %H:%M:%S", "%d.%m.%Y %H:%M"):
        try:
            return datetime.strptime(text, format_text)
        except ValueError:
            continue
    return None


def sort_batches_by_acquisition(batches: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    """Prefer acquisition chronology and retain the old O-number fallback.

    CHROMTAB exports can contain several consecutive runs and may arrange their
    sections lexicographically (O1, O1, O10, O10, ...).  Date Acquired is the
    only reliable way to keep the complete early run before a later run whose
    sample numbering starts again at O1.
    """
    indexed = list(enumerate(batches))
    if not indexed:
        return []

    parsed = [(index, batch, parse_acquired_at(batch.get("acquired_at"))) for index, batch in indexed]
    if any(acquired_at is not None for _, _, acquired_at in parsed):
        valid = [(index, batch, acquired_at) for index, batch, acquired_at in parsed if acquired_at is not None]
        missing = [(index, batch) for index, batch, acquired_at in parsed if acquired_at is None]
        valid.sort(key=lambda item: (item[2], item[0]))
        return [batch for _, batch, _ in valid] + [batch for _, batch in missing]

    return [
        batch
        for _, batch in sorted(
            indexed,
            key=lambda item: (sample_number_sort_key(item[1].get("sample_name", item[1].get("file_name", ""))), item[0]),
        )
    ]
