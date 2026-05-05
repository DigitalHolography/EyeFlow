from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import h5py
import numpy as np

from input_output.hdf5 import normalize_h5_path, set_attr_safe, write_value_dataset
from input_output.inputs import EyeFlowView, MergedAttrs
from input_output.schema import DOPPLER_VIEW_SCHEMA, HOLODOPPLER_SCHEMA, H5SourceSchema

from .base import DatasetValue, ProcessResult

_MISSING = object()


class H5SourceReader:
    """Schema-aware reader for one locked upstream HDF5 source."""

    def __init__(
        self,
        *,
        h5file: h5py.File | None,
        schema: H5SourceSchema,
        config: Mapping[str, object] | None = None,
    ) -> None:
        self.h5file = h5file
        self.schema = schema
        self.config_values = dict(config or {})

    @property
    def filename(self) -> str | None:
        if self.h5file is None:
            return None
        return self.h5file.filename

    def keys(self):
        if self.h5file is None:
            return ()
        return self.h5file.keys()

    def get(self, key_or_path: str, default=None):
        if self.h5file is None:
            return default
        path = self.path(key_or_path)
        found = self.h5file.get(path)
        return default if found is None else found

    def __getitem__(self, key_or_path: str):
        found = self.get(key_or_path)
        if found is None:
            raise KeyError(key_or_path)
        return found

    def __contains__(self, key_or_path: object) -> bool:
        return isinstance(key_or_path, str) and self.get(key_or_path) is not None

    def path(self, key_or_path: str) -> str:
        if key_or_path in self.schema.datasets:
            return self.schema.dataset_path(key_or_path)
        return normalize_h5_path(key_or_path)

    def dataset(self, key_or_path: str) -> h5py.Dataset:
        found = self.get(key_or_path)
        if not isinstance(found, h5py.Dataset):
            raise KeyError(
                f"Missing {self.schema.label} dataset '{key_or_path}'. "
                f"Expected: {self.path(key_or_path)}"
            )
        return found

    def value(self, key_or_path: str, default: Any = _MISSING):
        try:
            return self.dataset(key_or_path)[()]
        except KeyError:
            if default is not _MISSING:
                return default
            raise

    def array(
        self,
        key_or_path: str,
        *,
        dtype=None,
        flatten: bool = False,
        default: Any = _MISSING,
    ) -> np.ndarray:
        try:
            array = np.asarray(self.dataset(key_or_path)[()], dtype=dtype)
        except KeyError:
            if default is not _MISSING:
                return np.asarray(default, dtype=dtype)
            raise
        return np.ravel(array) if flatten else array

    def config(self, key: str, default: Any = _MISSING):
        spec = self.schema.config_value(key)
        value = self._h5_config_value(spec.h5_path)
        if value is not None:
            return value
        value = spec.read_json_config(dict(self.config_values))
        if value is not None:
            return value
        if default is not _MISSING:
            return default
        return spec.default

    def _h5_config_value(self, path: str | None):
        if self.h5file is None or path is None:
            return None
        found = self.h5file.get(path)
        if not isinstance(found, h5py.Dataset):
            return None
        array = np.asarray(found[()]).reshape(-1)
        if array.size == 0:
            return None
        value = array[0]
        if isinstance(value, bytes):
            return value.decode("utf-8")
        return value.item() if hasattr(value, "item") else value


class PipelineContext:
    """Runtime object passed to every pipeline."""

    def __init__(
        self,
        *,
        work_h5: h5py.File,
        holodoppler_h5: h5py.File | None,
        doppler_vision_h5: h5py.File | None,
        holodoppler_config: Mapping[str, object] | None = None,
        doppler_vision_config: Mapping[str, object] | None = None,
        preferred_input: str = "both",
        pipeline_name: str = "",
    ) -> None:
        self.work_h5 = work_h5
        self.work = work_h5
        self.ef = EyeFlowView(work_h5)
        self.hd = H5SourceReader(
            h5file=holodoppler_h5,
            schema=HOLODOPPLER_SCHEMA,
            config=holodoppler_config,
        )
        self.dv = H5SourceReader(
            h5file=doppler_vision_h5,
            schema=DOPPLER_VIEW_SCHEMA,
            config=doppler_vision_config,
        )
        self.hd_h5 = holodoppler_h5
        self.dv_h5 = doppler_vision_h5
        self.hd_config = dict(holodoppler_config or {})
        self.dv_config = dict(doppler_vision_config or {})
        self.preferred_input = preferred_input
        self.pipeline_name = pipeline_name
        self.attrs = MergedAttrs(
            self.work_h5,
            self._preferred_raw_source(),
            self._secondary_raw_source(),
            self.hd_config,
            self.dv_config,
        )

    @property
    def filename(self) -> str:
        primary = self._preferred_raw_source()
        if primary is not None and primary.filename is not None:
            return str(primary.filename)
        if self.work_h5.filename is not None:
            return str(self.work_h5.filename)
        return ""

    def _preferred_raw_source(self) -> h5py.File | None:
        if self.preferred_input == "dv":
            return self.dv_h5 or self.hd_h5
        return self.hd_h5 or self.dv_h5

    def _secondary_raw_source(self) -> h5py.File | None:
        preferred = self._preferred_raw_source()
        if preferred is self.hd_h5:
            return self.dv_h5
        if preferred is self.dv_h5:
            return self.hd_h5
        return None

    def get(self, path: str, default=None):
        normalized_path = normalize_h5_path(path)
        if not normalized_path:
            return default
        for source in (
            self.work_h5,
            self._preferred_raw_source(),
            self._secondary_raw_source(),
        ):
            found = source.get(normalized_path) if source is not None else None
            if found is not None:
                return found
        return default

    def __getitem__(self, path: str):
        found = self.get(path)
        if found is None:
            raise KeyError(path)
        return found

    def __contains__(self, path: object) -> bool:
        return isinstance(path, str) and self.get(path) is not None

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
    ) -> np.ndarray:
        try:
            array = np.asarray(self.read(path), dtype=dtype)
        except KeyError:
            if default is not _MISSING:
                return np.asarray(default, dtype=dtype)
            raise
        return np.ravel(array) if flatten else array

    def write(self, path: str, value: Any, **attrs: Any) -> None:
        payload = DatasetValue(value, attrs) if attrs else value
        write_value_dataset(self.work_h5, path, payload)

    def write_many(self, metrics: Mapping[str, Any]) -> None:
        for path, value in metrics.items():
            write_value_dataset(self.work_h5, path, value)

    def set_attr(self, key: str, value: Any) -> None:
        if key == "pipeline":
            return
        set_attr_safe(self.work_h5, key, value)

    def set_attrs(self, attrs: Mapping[str, Any] | None) -> None:
        for key, value in (attrs or {}).items():
            self.set_attr(str(key), value)

    def apply_result(self, result: ProcessResult | Mapping[str, Any] | None) -> None:
        if result is None:
            return
        if isinstance(result, ProcessResult):
            self.set_attrs(result.attrs)
            self.write_many(result.metrics)
            return
        if isinstance(result, Mapping):
            self.write_many(result)
            return
        raise TypeError(
            "Pipeline must return None, a metrics dict, or ProcessResult. "
            f"Got: {type(result).__name__}"
        )

    def finish_pipeline(self, pipeline_name: str) -> None:
        self.work_h5.attrs["last_pipeline"] = pipeline_name
        self.work_h5.flush()
