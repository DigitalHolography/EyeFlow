"""Source metadata and typed access helpers for EyeFlow HDF5 adapters."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from input_output.readers.nodes import is_dataset


def normalize_h5_path(path: object) -> str:
    """Use logical HDF5 separators independent of the host filesystem."""
    return str(path).replace("\\", "/").strip("/")


@dataclass(frozen=True)
class SourceFileLayout:
    """Companion application file layout metadata used for run discovery."""

    label: str
    companion_suffix: str
    h5_folder_name: str
    h5_filename_template: str
    config_dir_name: str | None = None
    config_filename: str | None = None


class TypedSource:
    """Thin typed facade over one raw HDF5 reader and its sidecar config."""

    INPUTS = {}

    def __init__(self, reader, config: dict[str, object] | None = None) -> None:
        self._reader = reader
        self._config = dict(config or {})

    @property
    def filename(self) -> str | None:
        return self._reader.filename

    def input_datasets(self, name: str) -> tuple:
        """Resolve one declarative input through the adapter's lazy resolver table."""
        try:
            resolve = self.INPUTS[name]
        except KeyError:
            raise KeyError(f"Unknown {self._reader.label} input declaration: {name}") from None
        value = resolve(self)
        nodes = value if isinstance(value, tuple) else (value,)
        return tuple(node for node in nodes if is_dataset(node))


def scalar_from_value(value):
    """Return the first scalar from a scalar-like HDF5/JSON value."""

    if value is None:
        return None
    array = np.asarray(value).reshape(-1)
    if array.size == 0:
        return None
    scalar = array[0]
    if isinstance(scalar, bytes):
        return scalar.decode("utf-8")
    return scalar.item() if hasattr(scalar, "item") else scalar
