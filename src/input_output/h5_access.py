"""HDF5 source readers and output access for pipeline runs."""

from collections.abc import Iterator, Mapping
from typing import Any

import h5py
import numpy as np

from .writers.h5 import normalize_h5_path, set_attr_safe, write_value_dataset

_MISSING = object()


class MergedAttrs(Mapping[str, object]):
    """Read ordered HDF5 and config attributes as one mapping."""

    def __init__(self, *sources: h5py.File | Mapping[str, object] | None) -> None:
        self._sources = [
            source.attrs if isinstance(source, h5py.File) else source
            for source in sources
            if source is not None
        ]

    def __getitem__(self, key: str) -> object:
        sentinel = object()
        value = self.get(key, sentinel)
        if value is sentinel:
            raise KeyError(key)
        return value

    def __iter__(self) -> Iterator[str]:
        seen: set[str] = set()
        for source in self._sources:
            for key in source.keys():
                if key not in seen:
                    seen.add(key)
                    yield str(key)

    def __len__(self) -> int:
        return sum(1 for _ in self.__iter__())

    def get(self, key: str, default=None):
        for source in self._sources:
            if key in source:
                return source[key]
        return default


class PipelineInputSource:
    """One explicit-path HDF5 reader with its parsed sidecar configuration."""

    def __init__(
        self,
        *,
        h5file: h5py.File | None,
        label: str,
        config: Mapping[str, object] | None = None,
    ) -> None:
        self.h5file = h5file
        self.label = label
        self.config = dict(config or {})

    @property
    def filename(self) -> str | None:
        if self.h5file is None:
            return None
        return self.h5file.filename

    @property
    def available(self) -> bool:
        return self.h5file is not None

    def keys(self):
        if self.h5file is None:
            return ()
        return self.h5file.keys()

    def get(self, path: str, default=None):
        if self.h5file is None:
            return default
        found = self.h5file.get(normalize_h5_path(path))
        return default if found is None else found

    def __contains__(self, path: object) -> bool:
        return isinstance(path, str) and self.get(path) is not None

    def dataset(self, path: str) -> h5py.Dataset:
        found = self.get(path)
        if not isinstance(found, h5py.Dataset):
            raise KeyError(
                f"Missing {self.label} dataset at path '{normalize_h5_path(path)}'."
            )
        return found

    def value(self, path: str, default: Any = _MISSING):
        try:
            return self.dataset(path)[()]
        except KeyError:
            if default is not _MISSING:
                return default
            raise

    def array(
        self,
        path: str,
        *,
        dtype=None,
        default: Any = _MISSING,
    ) -> Any:
        return _optional_dataset_array(
            self.get(path),
            dtype=dtype,
            default=default,
            missing_message=(
                f"Missing {self.label} dataset at path '{normalize_h5_path(path)}'."
            ),
        )

    def as_holodoppler(self):
        from .schema import HolodopplerSource

        return HolodopplerSource(self, self.config)

    def as_dopplerview(self):
        from .schema import DopplerViewSource

        return DopplerViewSource(self, self.config)


def _read_dataset_array(dataset: h5py.Dataset, *, dtype=None) -> np.ndarray:
    """Read numeric data in the requested dtype without a second full-size copy."""
    if dtype is None:
        return dataset[()]

    requested = np.dtype(dtype)
    if requested == np.dtype(bool):
        return np.asarray(dataset[()], dtype=requested)
    if requested == dataset.dtype:
        return dataset[()]
    if requested.kind in "iufc" and dataset.dtype.kind in "iufc":
        return dataset.astype(requested)[()]
    return np.asarray(dataset[()], dtype=requested)


def _optional_dataset_array(
    found,
    *,
    dtype,
    default: Any,
    missing_message: str,
) -> Any:
    if not isinstance(found, h5py.Dataset):
        if default is _MISSING:
            raise KeyError(missing_message)
        return None if default is None else np.asarray(default, dtype=dtype)
    return _read_dataset_array(found, dtype=dtype)


class PipelineH5Output:
    """Read and write the EyeFlow work/output HDF5 file."""

    def __init__(self, work_h5: h5py.File) -> None:
        self.file = work_h5

    @property
    def filename(self) -> str | None:
        return self.file.filename

    def get(self, path: str, default=None):
        found = self.file.get(normalize_h5_path(path))
        return default if found is None else found

    def read(self, path: str, default: Any = _MISSING):
        found = self.get(path)
        if not isinstance(found, h5py.Dataset):
            if default is not _MISSING:
                return default
            raise KeyError(path)
        return found[()]

    def array(
        self,
        path: str,
        *,
        dtype=None,
        default: Any = _MISSING,
    ) -> Any:
        return _optional_dataset_array(
            self.get(path), dtype=dtype, default=default, missing_message=path
        )

    def write_many(self, metrics: Mapping[str, Any]) -> None:
        for path, value in metrics.items():
            write_value_dataset(self.file, path, value)

    def set_attr(self, key: str, value: Any) -> None:
        if key == "pipeline":
            return
        set_attr_safe(self.file, key, value)

    def set_attrs(self, attrs: Mapping[str, Any] | None) -> None:
        for key, value in (attrs or {}).items():
            self.set_attr(str(key), value)

    def flush(self) -> None:
        self.file.flush()
