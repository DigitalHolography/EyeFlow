"""Dataset payload contract shared by producers and HDF5 storage."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass
class DatasetValue:
    """One dataset's data, metadata, and HDF5 creation options."""

    data: Any
    attrs: dict[str, Any] | None = None
    h5_options: dict[str, Any] | None = None


def dataset_parts(value: Any) -> tuple[Any, dict[str, Any] | None, dict[str, Any] | None]:
    """Interpret every supported dataset payload in one place."""
    if isinstance(value, DatasetValue):
        return value.data, value.attrs, value.h5_options
    if isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], dict):
        return value[0], value[1], None
    if hasattr(value, "data") and hasattr(value, "attrs"):
        return value.data, value.attrs, getattr(value, "h5_options", None)
    return value, None, None


def payload_data(value: Any) -> Any:
    return dataset_parts(value)[0]


def payload_array(value: Any) -> np.ndarray:
    return np.asarray(payload_data(value))
