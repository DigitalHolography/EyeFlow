"""Pipeline context namespaces for inputs, runtime state, outputs, and run logs."""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import h5py

from app_settings import (
    DEFAULT_VELOCITY_ESTIMATION_METHOD,
    VelocityEstimationMethod,
    validate_velocity_estimation_method,
)
from input_output.h5_access import PipelineH5Output, PipelineInputSource, RawH5SourceReader
from input_output.inputs import MergedAttrs
from input_output.output_manager import OutputManager
from utils.logger import Logger
from velocity_calibration import (
    DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
    validate_band_ratio_frequency_scale_hz,
)

from .base import ExecutionVariant, ProcessResult


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
        pipeline_options: Mapping[str, Sequence[str]] | None = None,
        pipeline_order: Sequence[str] = (),
        pipeline_targets: Sequence[str] = (),
        velocity_estimation_method: str | None = DEFAULT_VELOCITY_ESTIMATION_METHOD,
        band_ratio_frequency_scale_hz: float = (DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ),
        output_manager: OutputManager | None = None,
        processing_root: str | None = None,
        output_provenance: Mapping[str, Any] | None = None,
        execution_variant: ExecutionVariant | None = None,
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
        self.output = PipelineOutput(
            output_manager,
            PipelineH5Output(
                work_h5,
                processing_root=processing_root,
                provenance=output_provenance,
            ),
        )
        self.state = PipelineState(variables)
        self.pipeline_options = {
            str(name): frozenset(str(option) for option in options)
            for name, options in (pipeline_options or {}).items()
        }
        self.pipeline_order = tuple(str(name) for name in pipeline_order)
        self.pipeline_targets = tuple(str(name) for name in pipeline_targets)
        self.execution_variant = execution_variant
        self._execution_variants: tuple[ExecutionVariant, ...] = ()
        self._execution_variant_failures: dict[str, str] = {}
        resolved_velocity_method = velocity_estimation_method
        if resolved_velocity_method is None and execution_variant is not None:
            resolved_velocity_method = execution_variant.provenance.get(
                "velocity_estimation_method"
            )
        self.velocity_estimation_method: VelocityEstimationMethod | None = (
            validate_velocity_estimation_method(str(resolved_velocity_method))
            if resolved_velocity_method is not None
            else None
        )
        self.band_ratio_frequency_scale_hz = validate_band_ratio_frequency_scale_hz(
            band_ratio_frequency_scale_hz
        )
        self.attrs = MergedAttrs(
            work_h5,
            self._preferred_raw_source(),
            self._secondary_raw_source(),
            hd_config,
            dv_config,
        )

    @property
    def execution_variants(self) -> tuple[ExecutionVariant, ...]:
        """Variants published by the pipeline currently using this context."""

        return self._execution_variants

    @property
    def execution_variant_failures(self) -> dict[str, str]:
        """Mutable failure registry shared with the runtime after publication."""

        return self._execution_variant_failures

    def publish_execution_variants(
        self,
        variants: Sequence[ExecutionVariant],
        *,
        failures: Mapping[str, object] | None = None,
    ) -> None:
        """Publish isolated downstream executions from the current pipeline."""

        items = tuple(variants)
        if not items:
            raise ValueError("At least one execution variant must be published.")
        names: set[str] = set()
        for variant in items:
            if not isinstance(variant, ExecutionVariant):
                raise TypeError(
                    "Execution variants must be ExecutionVariant instances."
                )
            name = str(variant.name).strip()
            if not name:
                raise ValueError("Execution variant names cannot be empty.")
            if name in names:
                raise ValueError(f"Duplicate execution variant name: '{name}'")
            names.add(name)
            if not isinstance(variant.state, dict):
                raise TypeError(f"Execution variant '{name}' state must be a dict.")
            if not isinstance(variant.provenance, Mapping):
                raise TypeError(
                    f"Execution variant '{name}' provenance must be a mapping."
                )
            _validate_variant_namespace(
                variant.output_namespace,
                variant_name=name,
                label="output",
            )
            _validate_variant_namespace(
                variant.artifact_namespace,
                variant_name=name,
                label="artifact",
            )
        self._execution_variants = items
        self._execution_variant_failures = {
            str(name): str(message) for name, message in (failures or {}).items()
        }

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
        Logger.log(message)

    def log_debug(self, message: str) -> None:
        Logger.log_debug(message)

    def log_info(self, message: str) -> None:
        Logger.log_info(message)

    def log_warning(self, message: str) -> None:
        Logger.log_warning(message)

    def log_error(self, message: str) -> None:
        Logger.log_error(message)

    def option_enabled(self, name: str, *, pipeline: str | None = None) -> bool:
        pipeline_name = pipeline or self.runtime.pipeline_name
        return str(name) in self.pipeline_options.get(pipeline_name, frozenset())

    def options_for(self, pipeline: str) -> frozenset[str]:
        return self.pipeline_options.get(str(pipeline), frozenset())

    def pipeline_scheduled(self, pipeline: str) -> bool:
        return str(pipeline) in self.pipeline_order

    def pipeline_targeted(self, pipeline: str) -> bool:
        """Return whether a pipeline was selected directly, not as a dependency."""

        if not self.pipeline_targets:
            return self.pipeline_scheduled(pipeline)
        return str(pipeline) in self.pipeline_targets

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


def _validate_variant_namespace(
    namespace: str,
    *,
    variant_name: str,
    label: str,
) -> None:
    raw = str(namespace)
    value = raw.replace("\\", "/").strip("/")
    if (
        not value
        or raw != value
        or ":" in value
        or any(part in {"", ".", ".."} for part in value.split("/"))
    ):
        raise ValueError(
            f"Execution variant '{variant_name}' has an invalid {label} namespace: "
            f"{namespace!r}."
        )


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
        f"Pipeline must return None, a metrics dict, or ProcessResult. Got: {type(result).__name__}"
    )


def finish_pipeline(ctx: PipelineContext, pipeline_name: str) -> None:
    """Record completion metadata after a pipeline succeeds."""

    ctx.runtime.pipeline_name = pipeline_name
    ctx.runtime.work_h5.attrs["last_pipeline"] = pipeline_name
    ctx.output.h5.flush()
