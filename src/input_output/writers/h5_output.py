"""Pipeline-facing HDF5 output access and workflow writes."""

from collections.abc import Mapping
from typing import Any

import h5py
import numpy as np

from ..schema.base import normalize_h5_path
from ..schema.provenance import processing_summary_keys, read_dataset_attrs
from .h5 import set_attr_safe, write_value_dataset

_MISSING = object()


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
        self._write_guard = None

    def _path(self, path: str) -> str:
        from ..schema.eyeflow_output import processing_path

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

    def read_attrs(self, path: str) -> dict[str, object]:
        """Read dataset attributes including its workflow's shared provenance."""
        found = self.get(path)
        if not isinstance(found, h5py.Dataset):
            raise KeyError(path)
        return read_dataset_attrs(found)

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
        if self._write_guard is not None:
            self._write_guard()
        if self.processing_root:
            self.set_attrs(self.provenance)
        for path, value in metrics.items():
            target = self._path(path)
            in_workflow = bool(
                self.processing_root and target.startswith(self.processing_root + "/")
            )
            write_value_dataset(self.file, target, value)
            if in_workflow:
                dataset = self.file[target]
                # References in dataset attributes use the same workflow namespace.
                for key, attr in list(dataset.attrs.items()):
                    if isinstance(attr, str) and attr.lstrip("/").startswith("ProcessingMoments/"):
                        set_attr_safe(dataset, key, self._path(attr))

    def set_attr(self, key: str, value: Any) -> None:
        if self._write_guard is not None:
            self._write_guard()
        if key == "pipeline" or (
            self.processing_root and key not in processing_summary_keys(self.processing_root)
        ):
            return
        target = self.file.require_group(self.processing_root) if self.processing_root else self.file
        set_attr_safe(target, key, value)

    def set_attrs(self, attrs: Mapping[str, Any] | None) -> None:
        for key, value in (attrs or {}).items():
            self.set_attr(str(key), value)

    def flush(self) -> None:
        self.file.flush()
