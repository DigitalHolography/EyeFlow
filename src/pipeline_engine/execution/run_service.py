"""Shared GUI/CLI orchestration for EyeFlow pipeline runs."""

from __future__ import annotations

import os
import json
from time import perf_counter
from collections import deque
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack, contextmanager
from dataclasses import dataclass, replace
from pathlib import Path

from app_settings import (
    DEFAULT_BAND_LIMITED_HARMONICS_SETTING,
    DEFAULT_COPY_HD_PASSTHROUGH,
    DEFAULT_VELOCITY_WORKFLOW_SELECTION,
    validate_velocity_workflow_selection,
    validate_band_limited_harmonics_setting,
    validate_copy_hd_passthrough,
    validate_input_backend,
    validate_input_unknown_transport,
)
from input_output import HoloRunLayout, resolve_selected_run_layouts
from input_output.output_manager import (
    OutputManager, OutputType, output_manager_for_layout, reject_duplicate_destinations,
)
from input_output.runs.inputs import (
    ExpandedRunInputs,
    batch_root as infer_batch_root,
    expand_run_inputs,
)
from utils.logger import Logger
from velocity_calibration import (
    DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
    validate_band_ratio_frequency_scale_hz,
)

from pipeline_engine.planning.base import PipelineDescriptor, PipelineSelection
from pipeline_engine.execution.background_artifacts import ArtifactWorkers, BackgroundArtifacts
from pipeline_engine.planning.dag import PipelineDAG, PipelineExecutionPlan
from pipeline_engine.io.manager import IOManager
from pipeline_engine.execution.runtime import PipelineRunError, run_pipelines_to_output

@dataclass(frozen=True)
class RunRequest:
    """One resolved HOLO input and its final output destination."""

    input_layout: HoloRunLayout
    output_manager: OutputManager
    hd_h5: Path
    dv_h5: Path


@dataclass(frozen=True)
class RunSpec:
    """Fully resolved, front-end-independent batch run specification."""

    plan: PipelineExecutionPlan
    requests: tuple[RunRequest, ...]
    pipeline_options: Mapping[str, tuple[str, ...]]
    band_ratio_frequency_scale_hz: float = DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ
    band_limited_harmonics: int = DEFAULT_BAND_LIMITED_HARMONICS_SETTING
    copy_hd_passthrough: bool = DEFAULT_COPY_HD_PASSTHROUGH
    velocity_estimation_method: str = DEFAULT_VELOCITY_WORKFLOW_SELECTION
    input_backend: str = "auto"
    input_unknown_transport: str = "coordinator"

    @property
    def total_pipeline_units(self) -> int:
        return len(self.plan.descriptors) * len(self.requests)

    @property
    def selection(self):
        return PipelineSelection(self.plan.names, self.plan.targets, self.pipeline_options)


@dataclass(frozen=True)
class RunFailure:
    input_path: Path
    message: str

    def __str__(self) -> str:
        return f"{self.input_path}: {self.message}"


@dataclass(frozen=True)
class RunResult:
    outputs: tuple[Path, ...]
    failures: tuple[RunFailure, ...]
    stopped: bool = False

    @property
    def succeeded(self) -> bool:
        return not self.failures and not self.stopped

    @property
    def last_output_path(self) -> Path | None:
        return self.outputs[-1] if self.outputs else None


def selectable_pipeline_registry(
    pipelines: Iterable[PipelineDescriptor],
) -> dict[str, PipelineDescriptor]:
    """Return descriptors that a user may select directly."""

    return {pipeline.name: pipeline for pipeline in pipelines if pipeline.visibility != "hidden"}


def resolve_run_spec(
    *,
    input_paths: Sequence[Path],
    target_names: Sequence[str],
    pipelines: Iterable[PipelineDescriptor],
    pipeline_options: Mapping[str, Iterable[str]] | None = None,
    band_ratio_frequency_scale_hz: float = (DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ),
    band_limited_harmonics: int = DEFAULT_BAND_LIMITED_HARMONICS_SETTING,
    copy_hd_passthrough: bool = DEFAULT_COPY_HD_PASSTHROUGH,
    velocity_estimation_method: str = DEFAULT_VELOCITY_WORKFLOW_SELECTION,
    input_backend: str = "auto",
    input_unknown_transport: str = "coordinator",
    output_root: Path | None = None,
    batch_root: Path | None = None,
) -> RunSpec:
    """Resolve targets, inputs, and deterministic output destinations."""
    input_backend = validate_input_backend(input_backend)
    input_unknown_transport = validate_input_unknown_transport(input_unknown_transport)
    resolved_velocity_selection = validate_velocity_workflow_selection(velocity_estimation_method)
    resolved_harmonics = validate_band_limited_harmonics_setting(band_limited_harmonics)
    resolved_copy_hd_passthrough = validate_copy_hd_passthrough(copy_hd_passthrough)

    resolved_band_ratio_scale_hz = validate_band_ratio_frequency_scale_hz(
        band_ratio_frequency_scale_hz
    )
    descriptors = tuple(pipelines)
    selectable = selectable_pipeline_registry(descriptors)
    hidden_targets = [name for name in target_names if name not in selectable]
    if hidden_targets:
        raise ValueError("Unknown or hidden pipeline target(s): " + ", ".join(hidden_targets))

    resolved_options = _resolve_pipeline_options(descriptors, pipeline_options)
    with Logger.perf_stage("run/dag_resolution"):
        plan = PipelineDAG(descriptors).resolve_targets(
            target_names,
            pipeline_options=resolved_options,
        )
    if not plan.targets:
        raise ValueError("Select at least one pipeline target.")
    unavailable = [pipeline for pipeline in plan.descriptors if not pipeline.available]
    if unavailable:
        details = []
        for pipeline in unavailable:
            reason = ", ".join(pipeline.missing_deps or pipeline.requires)
            details.append(f"{pipeline.name}" + (f" ({reason})" if reason else ""))
        raise ValueError("The DAG requires unavailable pipeline(s): " + ", ".join(details))
    resolved_options = {
        descriptor.name: resolved_options[descriptor.name]
        for descriptor in plan.descriptors
        if descriptor.name in resolved_options
    }

    with Logger.perf_stage("run/input_discovery_validation"):
        layouts = resolve_selected_run_layouts(input_paths)
    resolved_output_root = output_root.expanduser().resolve() if output_root is not None else None
    effective_batch_root = (
        batch_root.expanduser().resolve()
        if batch_root is not None
        else infer_batch_root([layout.holo_path for layout in layouts])
    )
    with Logger.perf_stage("run/companion_path_snapshot"):
        requests = tuple(
            RunRequest(
                input_layout=layout,
                hd_h5=layout.hd_h5,
                dv_h5=layout.dv_h5,
                output_manager=output_manager_for_layout(
                    layout,
                    output_root=resolved_output_root,
                    batch_root=effective_batch_root,
                ),
            )
            for layout in layouts
        )
    reject_duplicate_destinations(
        tuple((request.input_layout, request.output_manager) for request in requests)
    )
    return RunSpec(
        plan=plan,
        requests=requests,
        pipeline_options=resolved_options,
        band_ratio_frequency_scale_hz=resolved_band_ratio_scale_hz,
        band_limited_harmonics=resolved_harmonics,
        copy_hd_passthrough=resolved_copy_hd_passthrough,
        velocity_estimation_method=resolved_velocity_selection,
        input_backend=input_backend,
        input_unknown_transport=input_unknown_transport,
    )


def execute_run(
    spec: RunSpec,
    *,
    on_file_start: Callable[[Path, int, int], None] | None = None,
    on_pipeline_start: Callable[[str, int, int], None] | None = None,
    on_progress: Callable[[], None] | None = None,
    should_stop: Callable[[], bool] | None = None,
) -> RunResult:
    """Execute a batch, optionally stopping between input files."""

    started_at = perf_counter()
    Logger.log(f"[DAG] Targets -> {', '.join(spec.plan.targets)}")
    Logger.log(f"[DAG] Execution order -> {', '.join(spec.plan.names)}")
    from pipeline_engine.resources import ResourceBudget, current_resources

    Logger.log_perf(f"START pipeline=run variant=shared stage=batch at={started_at:.6f}")
    with Logger.perf_stage("run/budget_create"):
        budget = current_resources() or ResourceBudget()
    with ExitStack() as startup:
        with Logger.perf_stage("run/batch_resources"):
            startup.enter_context(budget.activate())
        with ArtifactWorkers() as workers, ThreadPoolExecutor(
            max_workers=1, thread_name_prefix="EyeFlowFileCompletion",
        ) as completions:
            result = _execute_batch(
                spec, workers, completions, on_file_start=on_file_start,
                on_pipeline_start=on_pipeline_start, on_progress=on_progress,
                should_stop=should_stop,
            )
    finished_at = perf_counter()
    elapsed_seconds = finished_at - started_at
    Logger.log_perf(
        f"END pipeline=run variant=shared stage=batch "
        f"at={finished_at:.6f} wall={elapsed_seconds:.3f}s"
    )
    Logger.log_info(f"Total app execution time: {elapsed_seconds:.3f} s")
    return result


def _execute_batch(
    spec, workers, completions, *, on_file_start, on_pipeline_start, on_progress, should_stop,
):
    outputs = []
    failures = []
    pending = deque()
    stopped = False

    def finish_oldest():
        failure = pending.popleft().result()
        if failure is not None:
            failures.append(failure)
            for request in spec.requests:
                if request.input_layout.input_path == failure.input_path:
                    path = request.output_manager.path_for(OutputType.H5)
                    if path in outputs and not path.is_file():
                        outputs.remove(path)

    for request_index, request in enumerate(spec.requests):
        # Permit the next input's analysis to overlap the preceding exports,
        # while bounding private snapshots to at most two unfinished inputs.
        if len(pending) >= 2:
            finish_oldest()
        if request_index and should_stop is not None and should_stop():
            stopped = True
            remaining = len(spec.requests) - request_index
            Logger.log(f"[STOP] Batch stopped; {remaining} input file(s) not started.")
            break
        input_layout = request.input_layout
        final_manager = request.output_manager
        final_path = final_manager.path_for(OutputType.H5)
        if on_file_start is not None:
            on_file_start(
                input_layout.input_path,
                request_index + 1,
                len(spec.requests),
            )
        Logger.log(f"[INPUT] SOURCE -> {input_layout.input_path}")
        Logger.log(f"[INPUT] DATA DIR -> {input_layout.root_dir}")
        Logger.log(f"[RESOLVED] HD -> {request.hd_h5}")
        Logger.log(f"[RESOLVED] DV -> {request.dv_h5}")
        Logger.log(f"[OUTPUT] {final_path}")

        artifacts = BackgroundArtifacts(workers)
        file_started_at = perf_counter()
        numerical_error = None
        finalizers: list[IOManager] = []
        try:
            final_dir = final_manager.layout.ef_dir
            if final_dir.exists() and not final_dir.is_dir():
                raise RuntimeError(
                    f"Refusing to replace non-directory EyeFlow output path: {final_dir}"
                )
            # Start clean so a failed run leaves only its own partial output,
            # rather than mixing new files with artifacts from an old run.
            with Logger.perf_stage("run/output_directory_reset"):
                final_manager.prepare(replace=True)
            run_pipelines_to_output(
                output_manager=final_manager,
                pipelines=spec.plan.descriptors,
                target_names=spec.plan.targets,
                pipeline_options=spec.pipeline_options,
                band_ratio_frequency_scale_hz=(spec.band_ratio_frequency_scale_hz),
                band_limited_harmonics=spec.band_limited_harmonics,
                copy_hd_passthrough=spec.copy_hd_passthrough,
                input_backend=spec.input_backend,
                input_unknown_transport=spec.input_unknown_transport,
                velocity_estimation_method=spec.velocity_estimation_method,
                holodoppler_h5=request.hd_h5,
                doppler_vision_h5=request.dv_h5,
                on_pipeline_start=on_pipeline_start,
                on_progress=on_progress,
                background_artifacts=artifacts,
                on_finalization_ready=finalizers.append,
            )
        except Exception as exc:  # noqa: BLE001
            if isinstance(exc, PipelineRunError):
                outputs.append(exc.output_h5_path)
            numerical_error = str(exc)
        except BaseException:
            # Runtime may already have transferred its payload store during
            # interruption. It must have an owner even without a completion job.
            _discard_file_output(artifacts, finalizers)
            raise
        else:
            outputs.append(final_path)
        if artifacts.pending:
            phase = "ended with errors" if numerical_error else "completed"
            Logger.log(
                f"[ANALYSIS COMPLETE] {input_layout.input_path.name}: "
                f"pipeline execution {phase}; artifact exports are still in progress."
            )
        try:
            pending.append(completions.submit(
                _finish_file, input_layout.input_path, final_path, artifacts, numerical_error,
                finalizers, file_started_at,
            ))
        except BaseException:
            _discard_file_output(artifacts, finalizers)
            raise

    while pending:
        finish_oldest()
    return RunResult(tuple(outputs), tuple(failures), stopped=stopped)


def _discard_file_output(artifacts, outputs):
    with ExitStack() as retained, artifacts:
        for output in outputs:
            retained.callback(output.close)


def _finish_file(input_path, output_path, artifacts, numerical_error, finalizers=(), started_at=None):
    """Drain retained output and publish after final sidecar statuses are known."""
    import h5py

    metrics = finalizers[0].model.storage.metrics if finalizers else None
    try:
        with ExitStack() as retained, artifacts:
            for output in finalizers:
                retained.callback(output.close)
            artifact_wait_started = perf_counter()
            Logger.log_perf(
                f"START pipeline=run variant=shared stage=background_artifact_join "
                f"at={artifact_wait_started:.6f}"
            )
            artifacts.collect(wait=True)
            Logger.log_perf(
                f"END pipeline=run variant=shared stage=background_artifact_join "
                f"at={perf_counter():.6f} wall={perf_counter() - artifact_wait_started:.3f}s"
            )
            for output in finalizers:
                output.finish()
            if not finalizers and artifacts.statuses and output_path.is_file():
                with h5py.File(output_path, "r+") as h5:
                    h5.attrs["artifact_status"] = json.dumps(artifacts.statuses, sort_keys=True)
                    h5.attrs["artifact_failures"] = json.dumps(artifacts.failures, sort_keys=True)
            _log_performance_summary(
                input_path, output_path, started_at, artifacts, metrics,
            )
            errors = [numerical_error] if numerical_error else []
            for name, error in artifacts.failures.items():
                message = f"{name}: {error}"
                if not numerical_error or message not in numerical_error:
                    errors.append(message)
            if errors:
                failure = RunFailure(input_path, "; ".join(errors))
                Logger.log_error(str(failure))
                return failure
    except Exception as exc:
        try:
            _log_performance_summary(input_path, output_path, started_at, artifacts, metrics)
        except Exception as summary_exc:
            Logger.log_error(f"Could not produce performance summary for {input_path.name}: {summary_exc}")
        failure = RunFailure(input_path, str(exc))
        Logger.log_error(str(failure))
        return failure
    Logger.log(f"Completed run for {input_path.name}: {output_path}")
    return None


def _log_performance_summary(input_path, output_path, started_at, artifacts, metrics):
    import h5py

    finished_at = perf_counter()
    wall_seconds = max(0.0, finished_at - started_at) if started_at is not None else 0.0
    stage_timings = {}
    artifact_timings = artifacts.timings
    if output_path.is_file():
        with h5py.File(output_path, "r") as output:
            if "stage_timings" in output.attrs:
                stage_timings = json.loads(output.attrs["stage_timings"])
            if "artifact_timings" in output.attrs:
                artifact_timings = json.loads(output.attrs["artifact_timings"])

    stage_rows, pipeline_rows, intervals = [], [], []
    for scope, stages in stage_timings.items():
        scope_starts, scope_ends = [], []
        for stage, timing in stages.items():
            duration = float(timing.get("wall_seconds", 0.0))
            stage_rows.append((duration, scope, stage, timing))
            if "start" in timing and "end" in timing:
                start, end = float(timing["start"]), float(timing["end"])
                if end >= start:
                    scope_starts.append(start)
                    scope_ends.append(end)
                    intervals.append((start, end))
        if scope_starts:
            pipeline_rows.append((max(scope_ends) - min(scope_starts), scope))

    artifact_rows = []
    for label, timing in artifact_timings.items():
        start = timing.get("prepared_at")
        end = timing.get("finished_at")
        duration = max(0.0, float(end) - float(start)) if start is not None and end is not None else 0.0
        artifact_rows.append((duration, label, timing))
        if start is not None and end is not None and float(end) >= float(start):
            intervals.append((float(start), float(end)))

    Logger.log_perf(
        f"SUMMARY run={input_path.name} wall={wall_seconds:.3f}s "
        f"pipelines={len(pipeline_rows)} stages={len(stage_rows)} artifacts={len(artifact_rows)}"
    )
    for duration, scope in sorted(pipeline_rows, reverse=True)[:10]:
        Logger.log_perf(f"SUMMARY pipeline={scope} wall={duration:.3f}s")
    for duration, scope, stage, timing in sorted(stage_rows, reverse=True)[:10]:
        Logger.log_perf(
            f"SUMMARY pipeline={scope} stage={stage} wall={duration:.3f}s "
            f"callback={float(timing.get('callback_wall_seconds', 0.0)):.3f}s "
            f"resource_wait={float(timing.get('resource_wait_seconds', 0.0)):.3f}s "
            f"artifact_wait={float(timing.get('artifact_wait_seconds', 0.0)):.3f}s"
        )
    for duration, label, timing in sorted(artifact_rows, reverse=True):
        Logger.log_perf(
            f"SUMMARY artifact={label} wall={duration:.3f}s "
            f"admission_queue={float(timing.get('admission_queue_seconds', 0.0)):.3f}s "
            f"resource_wait={float(timing.get('resource_admission_wait_seconds', 0.0)):.3f}s "
            f"render_queue={float(timing.get('renderer_queue_seconds', 0.0)):.3f}s "
            f"render={float(timing.get('rendering_seconds', 0.0)):.3f}s "
            f"publication={float(timing.get('publication_seconds', 0.0)):.3f}s"
        )

    clipped = sorted(
        (max(started_at, start), min(finished_at, end))
        for start, end in intervals
        if started_at is not None and end >= started_at and start <= finished_at
    )
    covered = 0.0
    if clipped:
        merged_start, merged_end = clipped[0]
        for start, end in clipped[1:]:
            if start <= merged_end:
                merged_end = max(merged_end, end)
            else:
                covered += merged_end - merged_start
                merged_start, merged_end = start, end
        covered += merged_end - merged_start
    unaccounted = max(0.0, wall_seconds - covered)
    threshold = max(1.0, wall_seconds * 0.1)
    Logger.log_perf(
        f"SUMMARY accounted_wall={min(covered, wall_seconds):.3f}s "
        f"unaccounted_wall={unaccounted:.3f}s threshold={threshold:.3f}s"
    )
    Logger.log_perf(
        "SUMMARY note=operation durations may overlap; only the interval-union "
        "accounted_wall value is overlap-adjusted"
    )
    if unaccounted > threshold:
        Logger.log_perf(
            f"WARNING significant unaccounted execution time for {input_path.name}: "
            f"{unaccounted:.3f}s ({unaccounted / wall_seconds:.0%} of wall time). "
            "Stage/artifact intervals exclude run setup, inter-stage coordination, and finalization."
        )

    if metrics is not None:
        values = metrics.snapshot()
        Logger.log_perf(
            "SUMMARY storage "
            f"source_reads={int(values.get('source_reads', 0))} "
            f"array_copies={int(values.get('array_copies', 0))} "
            f"array_copied_bytes={int(values.get('array_copied_bytes', 0))} "
            f"scratch_resident_bytes={int(values.get('scratch_resident_allocation_bytes', 0))} "
            f"scratch_disk_bytes={int(values.get('scratch_disk_allocation_bytes', 0))} "
            f"h5_array_bytes={int(values.get('h5_logical_bytes', 0))} "
            f"value_serializations={int(values.get('value_serializations', 0))} "
            f"manifest_serializations={int(values.get('serializations', 0))}"
        )
        for name in sorted(
            key for key, value in values.items()
            if key.endswith("_seconds") and value > 0
        ):
            Logger.log_perf(f"SUMMARY operation={name} duration={values[name]:.3f}s")


def _resolve_pipeline_options(
    pipelines: Sequence[PipelineDescriptor],
    selections: Mapping[str, Iterable[str]] | None,
) -> dict[str, tuple[str, ...]]:
    requested_by_pipeline = selections or {}
    resolved: dict[str, tuple[str, ...]] = {}
    for descriptor in pipelines:
        if not descriptor.options:
            continue
        known = {option.name for option in descriptor.options}
        requested = requested_by_pipeline.get(descriptor.name)
        if requested is None:
            selected = {option.name for option in descriptor.options if option.default_enabled}
        else:
            selected = {str(name).strip() for name in requested if str(name).strip()}
            unknown = sorted(selected - known)
            if unknown:
                raise ValueError(
                    f"Unknown option(s) for pipeline '{descriptor.name}': " + ", ".join(unknown)
                )
        options_by_name = {option.name: option for option in descriptor.options}
        pending = list(selected)
        while pending:
            option_name = pending.pop()
            for required_name in options_by_name[option_name].requires:
                if required_name not in selected:
                    selected.add(required_name)
                    pending.append(required_name)
        resolved[descriptor.name] = tuple(
            option.name for option in descriptor.options if option.name in selected
        )
    return resolved


__all__ = [
    "ExpandedRunInputs",
    "RunFailure",
    "RunRequest",
    "RunResult",
    "RunSpec",
    "execute_run",
    "expand_run_inputs",
    "resolve_run_spec",
    "selectable_pipeline_registry",
]
