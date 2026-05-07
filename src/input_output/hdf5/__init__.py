"""Public HDF5 read/write helpers."""

from .core import (
    initialize_output_h5,
    normalize_h5_path,
    open_h5,
    resolve_dataset_target,
    set_attr_safe,
    write_value_dataset,
)

__all__ = [
    "initialize_output_h5",
    "normalize_h5_path",
    "open_h5",
    "resolve_dataset_target",
    "set_attr_safe",
    "write_value_dataset",
]
