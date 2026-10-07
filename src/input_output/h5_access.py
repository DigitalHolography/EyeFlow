"""HDF5 source readers and output access for pipeline runs."""

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import h5py
import numpy as np

from .writers.h5 import normalize_h5_path, set_attr_safe, write_value_dataset

_MISSING = object()


class RawH5SourceReader:
    """Explicit-path reader for one locked HDF5 source."""

    def __init__(self, *, h5file: h5py.File | None, label: str) -> None:
        self.h5file = h5file
        self.label = label

    @property
    def filename(self) -> str | None:
        if self.h5file is None:
            return None
        return self.h5file.filename

    @property
    def available(self) -> bool:
        return self.h5file is not None

    def require(self) -> None:
        if self.h5file is None:
            raise ValueError(f"{self.label} HDF5 input is required.")

    def keys(self):
        if self.h5file is None:
            return ()
        return self.h5file.keys()

    def get(self, path: str, default=None):
        if self.h5file is None:
            return default
        found = self.h5file.get(normalize_h5_path(path))
        return default if found is None else found

    def __getitem__(self, path: str):
        found = self.get(path)
        if found is None:
            raise KeyError(path)
        return found

    def __contains__(self, path: object) -> bool:
        return isinstance(path, str) and self.get(path) is not None

    def dataset(self, path: str) -> h5py.Dataset:
        found = self.get(path)
        if not isinstance(found, h5py.Dataset):
            raise KeyError(f"Missing {self.label} dataset at path '{normalize_h5_path(path)}'.")
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
        flatten: bool = False,
        default: Any = _MISSING,
    ) -> Any:
        try:
            array = _read_dataset_array(self.dataset(path), dtype=dtype)
        except KeyError:
            if default is _MISSING:
                raise
            if default is None:
                return None
            return np.asarray(default, dtype=dtype)
        return np.ravel(array) if flatten else array


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


@dataclass(frozen=True)
class PipelineInputSource:
    """One input HDF5 reader and its parsed sidecar configuration."""

    h5: RawH5SourceReader
    config: dict[str, object]

    @property
    def filename(self) -> str | None:
        return self.h5.filename

    @property
    def available(self) -> bool:
        return self.h5.available

    def require(self) -> None:
        self.h5.require()

    def keys(self):
        return self.h5.keys()

    def get(self, path: str, default=None):
        return self.h5.get(path, default)

    def dataset(self, path: str) -> h5py.Dataset:
        return self.h5.dataset(path)

    def value(self, path: str, default: Any = _MISSING):
        return self.h5.value(path, default)

    def array(
        self,
        path: str,
        *,
        dtype=None,
        flatten: bool = False,
        default: Any = _MISSING,
    ) -> Any:
        return self.h5.array(path, dtype=dtype, flatten=flatten, default=default)

    def as_holodoppler(self):
        from .schema import HolodopplerSource

        return HolodopplerSource(self.h5, self.config)

    def as_dopplerview(self):
        from .schema import DopplerViewSource

        return DopplerViewSource(self.h5, self.config)


class PipelineH5Output:
    """Read and write the EyeFlow work/output HDF5 file."""

    def __init__(
        self,
        work_h5: h5py.File,
        *,
        processing_root: str | None = None,
        provenance: Mapping[str, Any] | None = None,
    ) -> None:
        self.file = work_h5
        self.processing_root = processing_root
        self.provenance = dict(provenance or {})

    def _path(self, path: str) -> str:
        from .schema.eyeflow_output import processing_path

        if self.processing_root is None:
            return normalize_h5_path(path)
        return processing_path(path, self.processing_root)

    @property
    def attrs(self):
        if self.processing_root is None:
            return self.file.attrs
        return self.file[self.processing_root].attrs

    @property
    def filename(self) -> str | None:
        return self.file.filename

    def get(self, path: str, default=None):
        found = self.file.get(self._path(path))
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
        flatten: bool = False,
        default: Any = _MISSING,
    ) -> Any:
        try:
            array = np.asarray(self.read(path), dtype=dtype)
        except KeyError:
            if default is _MISSING:
                raise
            if default is None:
                return None
            return np.asarray(default, dtype=dtype)
        return np.ravel(array) if flatten else array

    def write(self, path: str, value: Any, **attrs: Any) -> None:
        payload = (value, attrs) if attrs else value
        self.write_many({path: payload})

    def write_many(self, metrics: Mapping[str, Any]) -> None:
        for path, value in metrics.items():
            target = self._path(path)
            write_value_dataset(self.file, target, value)
            if self.processing_root and target.startswith(self.processing_root + "/"):
                dataset = self.file[target]
                for key, attr in self.provenance.items():
                    set_attr_safe(dataset, key, attr)
                # References in dataset attributes use the same workflow namespace.
                for key, attr in list(dataset.attrs.items()):
                    if isinstance(attr, str) and attr.lstrip("/").startswith("Processing/"):
                        set_attr_safe(dataset, key, self._path(attr))

    def set_attr(self, key: str, value: Any) -> None:
        if key == "pipeline":
            return
        target = (
            self.file
            if self.processing_root is None
            else self.file.require_group(self.processing_root)
        )
        set_attr_safe(target, key, value)

    def set_attrs(self, attrs: Mapping[str, Any] | None) -> None:
        for key, value in (attrs or {}).items():
            self.set_attr(str(key), value)

    def flush(self) -> None:
        self.file.flush()
