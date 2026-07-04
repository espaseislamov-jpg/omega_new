from omega_path_compat import configure_windows_path_compat

configure_windows_path_compat()

from .current_engine import (
    load_batches,
    load_reference_targets,
    process_batch,
    process_file,
)
from .pipeline import process_from_baseline

__all__ = [
    "load_batches",
    "load_reference_targets",
    "process_batch",
    "process_file",
    "process_from_baseline",
]
