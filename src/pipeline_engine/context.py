from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import h5py
import numpy as np

from input_output.inputs import EyeFlowView, MergedAttrs
from input_output.output_manager import OutputManager
from input_output.writers.h5 import (
    normalize_h5_path,
    set_attr_safe,
    write_value_dataset,
)

from .base import DatasetValue, ProcessResult

_MISSING = object()


class RawH5SourceReader:
    """Explicit-path reader for one locked HDF5 source."""

    def __init__(
        self,
        *,
        h5file: h5py.File | None,
        label: str,
    ) -> None:
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
        flatten: bool = False,
        default: Any = _MISSING,
    ) -> Any:
        try:
            array = np.asarray(self.dataset(path)[()], dtype=dtype)
        except KeyError:
            if default is not _MISSING:
                if default is None:
                    return None
                return np.asarray(default, dtype=dtype)
            raise
        return np.ravel(array) if flatten else array


@dataclass(frozen=True)
class PipelineSources:
    """Raw explicit-path source readers available on PipelineContext."""

    work: RawH5SourceReader
    hd: RawH5SourceReader
    dv: RawH5SourceReader


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
        variables: dict[str, Any] | None = None,
        output_manager: OutputManager | None = None,
    ) -> None:
        self.work_h5 = work_h5
        self.work = work_h5
        self.ef = EyeFlowView(work_h5)
        self.sources = PipelineSources(
            work=RawH5SourceReader(h5file=work_h5, label="work"),
            hd=RawH5SourceReader(h5file=holodoppler_h5, label="HD"),
            dv=RawH5SourceReader(h5file=doppler_vision_h5, label="DV"),
        )
        self.hd_h5 = holodoppler_h5
        self.dv_h5 = doppler_vision_h5
        self.hd_config = dict(holodoppler_config or {})
        self.dv_config = dict(doppler_vision_config or {})
        self.preferred_input = preferred_input
        self.output = output_manager
        self.pipeline_name = pipeline_name
        # Shared in-memory state for this run. It is not persisted unless a
        # pipeline explicitly writes a value to the work H5.
        self.vars = variables if variables is not None else {}
        self.state = self.vars
        self.attrs = MergedAttrs(
            self.work_h5,
            self._preferred_raw_source(),
            self._secondary_raw_source(),
            self.hd_config,
            self.dv_config,
        )

    def require_inputs(self, *inputs: str) -> None:
        requested = {name.lower() for name in inputs} or {"hd", "dv"}
        missing: list[str] = []
        if "hd" in requested and not self.sources.hd.available:
            missing.append("HD")
        if "dv" in requested and not self.sources.dv.available:
            missing.append("DV")
        if missing:
            raise ValueError(f"Missing required input(s): {', '.join(missing)}.")

    def set_var(self, key: str, value: Any) -> None:
        self.vars[str(key)] = value

    def get_var(self, key: str, default: Any = None) -> Any:
        return self.vars.get(str(key), default)

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
