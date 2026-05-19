"""Small source metadata and reader helpers for EyeFlow HDF5 adapters."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np


MISSING = object()


@dataclass(frozen=True)
class SourceFileLayout:
    """Companion application file layout metadata used for run discovery."""

    label: str
    companion_suffix: str
    h5_folder_name: str
    h5_filename_template: str
    config_dir_name: str | None = None
    config_filename: str | None = None


@dataclass(frozen=True)
class HolodopplerTiming:
    """Sampling metadata exported by Holodoppler."""

    sampling_freq: float
    batch_stride: float

    @property
    def dt_seconds(self) -> float:
        return self.batch_stride / self.sampling_freq


class TypedSource:
    """Thin typed facade over one raw HDF5 reader and its sidecar config."""

    layout: SourceFileLayout

    def __init__(self, reader, config: dict[str, object] | None = None) -> None:
        self._reader = reader
        self._config = dict(config or {})

    @property
    def filename(self) -> str | None:
        return self._reader.filename

    def require(self) -> None:
        self._reader.require()

    def _array(
        self,
        path: str,
        *,
        dtype=None,
        default: Any = MISSING,
    ) -> Any:
        return self._reader.array(path, dtype=dtype, default=default)

    def _value(self, path: str, *, default: Any = MISSING):
        return self._reader.value(path, default=default)

    def _scalar_h5_or_config(self, h5_path: str, config_key: str):
        value = scalar_from_value(self._value(h5_path, default=None))
        if value is not None:
            return value
        return scalar_from_value(self._config.get(config_key))

    def _config_value(self, section: str, key: str, default):
        source = self._config.get(section, {})
        if not isinstance(source, dict):
            return default
        return source.get(key, default)

    def config_value(self, section: str, key: str, default=None):
        return self._config_value(section, key, default)


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
