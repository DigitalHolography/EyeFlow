"""Pipeline context namespaces for inputs, runtime state, outputs, and run logs."""

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any

import h5py
import numpy as np

from input_output.inputs import MergedAttrs
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
        from input_output.schema import HolodopplerSource

        return HolodopplerSource(self.h5, self.config)

    def as_dopplerview(self):
        from input_output.schema import DopplerViewSource

        return DopplerViewSource(self.h5, self.config)


@dataclass(frozen=True)
class PipelineInputs:
    """External application inputs consumed by EyeFlow."""

    hd: PipelineInputSource
    dv: PipelineInputSource


class PipelineState:
    """Shared in-memory values for the current pipeline run."""

    def __init__(self, values: dict[str, Any] | None = None) -> None:
        self._values = values if values is not None else {}

    def set(self, key: str, value: Any) -> None:
        self._values[str(key)] = value

    def get(self, key: str, default: Any = None) -> Any:
        return self._values.get(str(key), default)

    def __contains__(self, key: object) -> bool:
        return str(key) in self._values if isinstance(key, str) else False

    def __getitem__(self, key: str) -> Any:
        return self._values[key]

    @property
    def raw(self) -> dict[str, Any]:
        return self._values


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
        write_value_dataset(self.file, path, payload)

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


@dataclass(frozen=True)
class PipelineOutput:
    """Output namespace for the work H5 and sidecar artifacts."""

    manager: OutputManager | None
    h5: PipelineH5Output

    @property
    def available(self) -> bool:
        return self.manager is not None

    def dir_for(self, output_type):
        return self._manager().dir_for(output_type)

    def path_for(self, output_type, filename: str | None = None):
        return self._manager().path_for(output_type, filename)

    def open_h5(self, filename: str | None = None, mode: str = "w"):
        return self._manager().open_h5(filename, mode)

    def write_sidecar(self, output, output_type, filename: str | None = None):
        return self._manager().write_sidecar(output, output_type, filename)

    def write_json(self, output, filename: str | None = None):
        return self._manager().write_json(output, filename)

    def write_png(self, output, filename: str | None = None):
        return self._manager().write_png(output, filename)

    def _manager(self) -> OutputManager:
        if self.manager is None:
            raise ValueError("No output manager is available for this pipeline run.")
        return self.manager


@dataclass
class PipelineRuntime:
    """Pipeline engine bookkeeping for the current run."""

    work_h5: h5py.File
    preferred_input: str
    pipeline_name: str


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
        on_log: Callable[[str], None] | None = None,
    ) -> None:
        hd_config = dict(holodoppler_config or {})
        dv_config = dict(doppler_vision_config or {})
        self.runtime = PipelineRuntime(work_h5, preferred_input, pipeline_name)
        self.inputs = PipelineInputs(
            hd=PipelineInputSource(
                RawH5SourceReader(h5file=holodoppler_h5, label="HD"),
                hd_config,
            ),
            dv=PipelineInputSource(
                RawH5SourceReader(h5file=doppler_vision_h5, label="DV"),
                dv_config,
            ),
        )
        self.output = PipelineOutput(output_manager, PipelineH5Output(work_h5))
        self.state = PipelineState(variables)
        self._on_log = on_log
        self.attrs = MergedAttrs(
            work_h5,
            self._preferred_raw_source(),
            self._secondary_raw_source(),
            hd_config,
            dv_config,
        )

    def require_inputs(self, *inputs: str) -> None:
        requested = {name.lower() for name in inputs} or {"hd", "dv"}
        missing: list[str] = []
        if "hd" in requested and not self.inputs.hd.available:
            missing.append("HD")
        if "dv" in requested and not self.inputs.dv.available:
            missing.append("DV")
        if missing:
            raise ValueError(f"Missing required input(s): {', '.join(missing)}.")

    def log(self, message: str) -> None:
        if self._on_log is not None:
            self._on_log(message)

    @property
    def filename(self) -> str:
        primary = self._preferred_raw_source()
        if primary is not None and primary.filename is not None:
            return str(primary.filename)
        if self.runtime.work_h5.filename is not None:
            return str(self.runtime.work_h5.filename)
        return ""

    def _preferred_raw_source(self) -> h5py.File | None:
        hd_h5 = self.inputs.hd.h5.h5file
        dv_h5 = self.inputs.dv.h5.h5file
        if self.runtime.preferred_input == "dv":
            return dv_h5 or hd_h5
        return hd_h5 or dv_h5

    def _secondary_raw_source(self) -> h5py.File | None:
        hd_h5 = self.inputs.hd.h5.h5file
        dv_h5 = self.inputs.dv.h5.h5file
        preferred = self._preferred_raw_source()
        if preferred is hd_h5:
            return dv_h5
        if preferred is dv_h5:
            return hd_h5
        return None


def apply_pipeline_result(
    ctx: PipelineContext,
    result: ProcessResult | Mapping[str, Any] | None,
) -> None:
    """Persist a pipeline return value to the output H5."""

    if result is None:
        return
    if isinstance(result, ProcessResult):
        ctx.output.h5.set_attrs(result.attrs)
        ctx.output.h5.write_many(result.metrics)
        return
    if isinstance(result, Mapping):
        ctx.output.h5.write_many(result)
        return
    raise TypeError(
        "Pipeline must return None, a metrics dict, or ProcessResult. "
        f"Got: {type(result).__name__}"
    )


def finish_pipeline(ctx: PipelineContext, pipeline_name: str) -> None:
    """Record completion metadata after a pipeline succeeds."""

    ctx.runtime.pipeline_name = pipeline_name
    ctx.runtime.work_h5.attrs["last_pipeline"] = pipeline_name
    ctx.output.h5.flush()
