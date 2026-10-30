"""Pipeline context namespaces for inputs, runtime state, outputs, and run logs."""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Any

import h5py

from app_settings import (
    DEFAULT_VELOCITY_ESTIMATION_METHOD,
    DEFAULT_VELOCITY_WORKFLOW_SELECTION,
    selected_velocity_methods,
    VelocityEstimationMethod,
    validate_velocity_estimation_method,
)
from input_output.readers.h5 import MergedAttrs, PipelineInputSource, RawH5SourceReader
from input_output.writers.h5_output import PipelineH5Output
from input_output.output_manager import OutputManager
from pipeline_engine.scratch import RunArrayStorage
from input_output.writers.h5 import write_root_attributes
from utils.logger import Logger
from velocity_calibration import (
    DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
    validate_band_ratio_frequency_scale_hz,
)

from .base import ExecutionVariant, ProcessResult
from .background_artifacts import BackgroundArtifacts


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
    scratch: RunArrayStorage | None = None


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
        shared_inputs: PipelineInputs | None = None,
        preferred_input: str = "both",
        pipeline_name: str = "",
        variables: dict[str, Any] | None = None,
        pipeline_options: Mapping[str, Sequence[str]] | None = None,
        pipeline_order: Sequence[str] = (),
        pipeline_targets: Sequence[str] = (),
        velocity_estimation_method: str | None = DEFAULT_VELOCITY_ESTIMATION_METHOD,
        velocity_workflow_selection: str = DEFAULT_VELOCITY_WORKFLOW_SELECTION,
        band_ratio_frequency_scale_hz: float = (DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ),
        output_manager: OutputManager | None = None,
        processing_root: str | None = None,
        output_provenance: Mapping[str, Any] | None = None,
        execution_variant: ExecutionVariant | None = None,
        scratch: RunArrayStorage | None = None,
        background_artifacts: BackgroundArtifacts | None = None,
        computation_cache: dict[object, object] | None = None,
    ) -> None:
        hd_config = dict(holodoppler_config or {})
        dv_config = dict(doppler_vision_config or {})
        self.runtime = PipelineRuntime(work_h5, preferred_input, pipeline_name, scratch)
        self.inputs = shared_inputs or PipelineInputs(
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
        # Pure, completed computations may be reused across execution variants.
        # This per-input memo is separate from transactional pipeline products:
        # callers own compatibility keys and must store only immutable successes.
        self.computation_cache = computation_cache if computation_cache is not None else {}
        self.pipeline_options = {
            str(name): frozenset(str(option) for option in options)
            for name, options in (pipeline_options or {}).items()
        }
        self.pipeline_order = tuple(str(name) for name in pipeline_order)
        self.pipeline_targets = tuple(str(name) for name in pipeline_targets)
        self.execution_variant = execution_variant
        self._background_artifacts = background_artifacts
        self._pending_artifacts = []
        self._artifact_jobs = []
        self._figure_job = None
        self._figure_count = 0
        self._dag_execution = None
        self._graph_runtime_available = False
        self._execution_variants: tuple[ExecutionVariant, ...] = ()
        self._execution_variant_failures: dict[str, str] = {}
        self.velocity_workflows = selected_velocity_methods(velocity_workflow_selection)
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
    def dag(self):
        """The graph for this invocation; copies/workflows never borrow another's."""
        if self._dag_execution is None or self._dag_execution.context is not self:
            raise RuntimeError("No active pipeline DAG belongs to this execution context.")
        return self._dag_execution

    def __copy__(self):
        clone = object.__new__(type(self))
        clone.__dict__ = self.__dict__.copy()
        clone._dag_execution = None
        clone._pending_artifacts = []
        clone._artifact_jobs = []
        clone._figure_job = None
        clone._figure_count = 0
        return clone

    @property
    def can_defer_artifacts(self) -> bool:
        return self._background_artifacts is not None and self.output.available

    def defer_artifacts(self, name, exporter, snapshotter, **kwargs) -> None:
        """Snapshot now; render in a separate process after this pipeline commits.

        The exporter must be an importable function accepting a snapshot path,
        a staging OutputManager, and picklable keyword arguments. It returns
        the paths it created. The snapshotter runs here on the owning thread.
        """
        if not self.can_defer_artifacts:
            raise RuntimeError("This context has no background artifact lifecycle.")
        label = f"{self.runtime.pipeline_name}:{name}"
        if self._dag_execution is not None and self._dag_execution.active_node:
            label = f"{self.runtime.pipeline_name}/{self._dag_execution.active_node}:{name}"
        if self.execution_variant is not None and self.output.h5.processing_root:
            label += f" [{self.execution_variant.name}]"
        job = self._background_artifacts.prepare(
            label, self.output.manager, exporter, snapshotter, **kwargs,
        )
        self._pending_artifacts.append(job)
        self._artifact_jobs.append(job)

    def _capture_artifact_write(self, kind, payload, paths, kwargs) -> bool:
        """Shared writers use this automatically while a descriptor executes."""
        if not self.can_defer_artifacts:
            return False
        root = self.output.manager.layout.ef_dir.resolve()
        if any(not path.resolve().is_relative_to(root) for path in paths):
            return False
        from input_output.writers.deferred import export_figure_writes, snapshot_figure_write

        if self._figure_job is None:
            self.defer_artifacts("figures", export_figure_writes, lambda path: path.mkdir())
            self._figure_job = self._pending_artifacts[-1]
        snapshot = self._figure_job.snapshot / f"{self._figure_count:06d}.pickle"
        try:
            snapshot_figure_write(
                snapshot, kind, payload, [path.resolve().relative_to(root) for path in paths], kwargs,
            )
        except Exception as exc:
            # Serialization is export work: retain numerical results even if a
            # custom artist/formatter cannot be transferred to the renderer.
            snapshot.unlink(missing_ok=True)
            self._figure_job.kwargs.setdefault("errors", []).append(
                f"Cannot snapshot {paths[0].name}: {type(exc).__name__}: {exc}"
            )
        self._figure_count += 1
        return True

    @property
    def artifacts_pending(self) -> bool:
        return any(job.status in {"prepared", "running"} for job in self._artifact_jobs)

    def _start_artifacts(self) -> None:
        for job in self._pending_artifacts:
            self._background_artifacts.start(job)
        self._pending_artifacts.clear()

    def _cancel_artifacts(self) -> None:
        for job in self._pending_artifacts:
            self._background_artifacts.cancel(job)
        self._pending_artifacts.clear()

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
    write_root_attributes(ctx.runtime.work_h5, {"last_pipeline": pipeline_name})
    ctx.output.h5.flush()
