"""Runtime execution helpers for resolved EyeFlow pipelines."""

import json
import multiprocessing
from collections.abc import Callable, Mapping, Sequence
from contextlib import ExitStack
from dataclasses import replace
from pathlib import Path

from app_settings import (
    DEFAULT_BAND_LIMITED_HARMONICS_SETTING,
    DEFAULT_COPY_HD_PASSTHROUGH,
    DEFAULT_VELOCITY_WORKFLOW_SELECTION,
    validate_velocity_workflow_selection,
    selected_velocity_methods,
    validate_band_limited_harmonics_setting,
    validate_copy_hd_passthrough,
)
from input_output.readers.h5 import PipelineInputSource, RawH5SourceReader
from input_output.readers.h5 import load_h5_sidecar_config
from input_output.output_manager import OutputManager, OutputType
from pipeline_engine.scratch import RunArrayStorage
from input_output.writers.h5 import (
    open_h5, write_root_attributes,
)
from input_output.writers.h5_initialization import initialize_output_h5
from utils.logger import Logger
from velocity_calibration import (
    DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
    validate_band_ratio_frequency_scale_hz,
)

from .base import ExecutionVariant, PipelineDescriptor
from .background_artifacts import BackgroundArtifacts
from .context import PipelineInputs
from .dag import PipelineDAG
from .invocation import _run_pipeline_descriptor
from .storage import ProductStorage
from .output_model import OutputModel
from .resources import ResourceBudget, current_resources
from .input_service import InputService
from .scheduler import CoordinatorScheduler
from .process_executor import ProcessExecutor
from .protocol import InvocationSettings
from .finalization import PendingOutput


class PipelineRunError(RuntimeError):
    """An incomplete run whose successful products remain available."""

    def __init__(self, message: str, output_h5_path: Path):
        super().__init__(message)
        self.output_h5_path = output_h5_path


def run_pipelines_to_output(
    *,
    output_manager: OutputManager,
    pipelines: Sequence[PipelineDescriptor],
    target_names: Sequence[str] = (),
    pipeline_options: Mapping[str, Sequence[str]] | None = None,
    band_ratio_frequency_scale_hz: float = (DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ),
    band_limited_harmonics: int = DEFAULT_BAND_LIMITED_HARMONICS_SETTING,
    copy_hd_passthrough: bool = DEFAULT_COPY_HD_PASSTHROUGH,
    velocity_estimation_method: str = DEFAULT_VELOCITY_WORKFLOW_SELECTION,
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
    on_pipeline_start: Callable[[str, int, int], None] | None = None,
    on_pipeline_success: Callable[[str], None] | None = None,
    on_progress: Callable[[], None] | None = None,
    background_artifacts: BackgroundArtifacts | None = None,
    on_finalization_ready: Callable[[PendingOutput], None] | None = None,
) -> Path:
    """Run resolved pipelines and write outputs through an OutputManager."""
    resolved_velocity_selection = validate_velocity_workflow_selection(velocity_estimation_method)
    resolved_harmonics = validate_band_limited_harmonics_setting(band_limited_harmonics)
    resolved_copy_hd_passthrough = validate_copy_hd_passthrough(copy_hd_passthrough)

    resolved_band_ratio_scale_hz = validate_band_ratio_frequency_scale_hz(
        band_ratio_frequency_scale_hz
    )
    output_manager.prepare()
    output_h5_path = output_manager.path_for(OutputType.H5)
    budget = current_resources() or ResourceBudget()
    with budget.activate():
        storage = ProductStorage(budget=budget)
        work_h5 = OutputModel(storage, output_h5_path)
        wait_for_exports = background_artifacts is None
        background_artifacts = background_artifacts or BackgroundArtifacts()
        pending = PendingOutput(work_h5, background_artifacts, owns_artifacts=wait_for_exports)
        try:
            with ExitStack() as stack:
                scratch = stack.enter_context(RunArrayStorage(output_manager.layout.ef_dir))
                return _execute_input(
                    work_h5=work_h5,
                    output_h5_path=output_h5_path,
                    output_manager=output_manager,
                    stack=stack,
                    scratch=scratch,
                    pipelines=pipelines,
                    target_names=target_names,
                    pipeline_options=pipeline_options or {},
                    band_ratio_frequency_scale_hz=resolved_band_ratio_scale_hz,
                    band_limited_harmonics=resolved_harmonics,
                    copy_hd_passthrough=resolved_copy_hd_passthrough,
                    velocity_workflow_selection=resolved_velocity_selection,
                    holodoppler_h5=holodoppler_h5,
                    doppler_vision_h5=doppler_vision_h5,
                    on_pipeline_start=on_pipeline_start,
                    on_pipeline_success=on_pipeline_success,
                    on_progress=on_progress,
                    background_artifacts=background_artifacts,
                    wait_for_exports=wait_for_exports,
                )
        finally:
            # Workers and input mappings have closed before payload ownership
            # is either transferred to the batch or completed here.
            if on_finalization_ready is None:
                pending.finish()
            else:
                try:
                    on_finalization_ready(pending)
                except BaseException:
                    pending.close()
                    raise



def _execute_input(
    *,
    work_h5,
    output_h5_path,
    output_manager: OutputManager,
    stack: ExitStack,
    scratch: RunArrayStorage,
    pipelines: Sequence[PipelineDescriptor],
    target_names: Sequence[str],
    pipeline_options: Mapping[str, Sequence[str]],
    band_ratio_frequency_scale_hz: float,
    band_limited_harmonics: int,
    copy_hd_passthrough: bool,
    velocity_workflow_selection: str,
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
    on_pipeline_start: Callable[[str, int, int], None] | None,
    on_pipeline_success: Callable[[str], None] | None,
    on_progress: Callable[[], None] | None,
    background_artifacts: BackgroundArtifacts,
    wait_for_exports: bool,
) -> Path:
    hd_h5, dv_h5 = _open_input_h5_sources(stack, holodoppler_h5, doppler_vision_h5)
    _initialize_work_h5(
        work_h5=work_h5,
        pipelines=pipelines,
        target_names=target_names,
        pipeline_options=pipeline_options,
        band_ratio_frequency_scale_hz=band_ratio_frequency_scale_hz,
        copy_hd_passthrough=copy_hd_passthrough,
        velocity_workflow_selection=velocity_workflow_selection,
        holodoppler_h5=holodoppler_h5,
        doppler_vision_h5=doppler_vision_h5,
        hd_source_h5=hd_h5,
        dv_source_h5=dv_h5,
    )
    input_service = InputService({"hd": hd_h5, "dv": dv_h5}, work_h5.storage,
                                 multiprocessing.get_context("spawn"))
    stack.callback(input_service.close)
    input_client = input_service.client()
    stack.callback(input_client.close)
    hd_config = load_h5_sidecar_config(hd_h5, source="hd")
    write_root_attributes(work_h5, {"band_limited_harmonics": band_limited_harmonics})
    dv_config = load_h5_sidecar_config(dv_h5, source="dv")
    files = input_service.files(input_client)
    hd_h5, dv_h5 = files["hd"], files["dv"]
    shared_inputs = PipelineInputs(
        hd=PipelineInputSource(RawH5SourceReader(h5file=hd_h5, label="HD"), hd_config),
        dv=PipelineInputSource(RawH5SourceReader(h5file=dv_h5, label="DV"), dv_config),
    )
    context_vars: dict[str, object] = {}
    computation_cache: dict[object, object] = {}

    dag = PipelineDAG(pipelines)
    execution_variants: list[ExecutionVariant] = []
    variant_dependents = {
        name
        for producer in pipelines if producer.produces_execution_variants
        for name in dag.dependents_of(
            producer.name, transitive=True, pipeline_options=pipeline_options
        )
    }
    variant_dependents.difference_update(
        pipeline.name for pipeline in pipelines if pipeline.run_once
    )
    variant_failures: dict[str, str] = {}
    shared_status: dict[str, str] = {}
    branch_status: dict[str, dict[str, str]] = {}
    pipeline_failures: dict[str, str] = {}
    pipeline_skips: dict[str, str] = {}

    pipeline_count = len(pipelines)
    executor = ProcessExecutor(
        context=multiprocessing.get_context("spawn"), input_service=input_service,
        inputs=files, model=work_h5, manager=output_manager, budget=current_resources(),
        state=context_vars,
        settings=InvocationSettings(
            configs={"hd": hd_config, "dv": dv_config},
            pipeline_options=dict(pipeline_options),
            pipeline_order=tuple(pipeline.name for pipeline in pipelines),
            pipeline_targets=tuple(target_names),
            velocity_workflow_selection=velocity_workflow_selection,
            band_ratio_frequency_scale_hz=band_ratio_frequency_scale_hz,
        ),
    )
    stack.callback(executor.close)
    scheduler = CoordinatorScheduler(
        pipelines=[replace(p, execution_mode="coordinator")
                   if p.name in variant_dependents or p.produces_execution_variants
                   or p.wait_for_artifacts
                   else p for p in pipelines],
        dag=dag, options=pipeline_options,
        executor=executor, statuses=shared_status, failures=pipeline_failures, skips=pipeline_skips,
        on_start=on_pipeline_start, on_success=on_pipeline_success, on_progress=on_progress,
    )
    for pipeline_index, pipeline_desc in scheduler.coordinator_jobs():
        # Only consumers of artifacts wait; ordinary calculations keep running.
        background_artifacts.collect(wait=pipeline_desc.wait_for_artifacts)
        _write_artifact_status(work_h5, background_artifacts)
        if pipeline_desc.wait_for_artifacts:
            work_h5.flush()
        if on_pipeline_start is not None:
            on_pipeline_start(pipeline_desc.name, pipeline_index, pipeline_count)
        arguments = {
            "work_h5": work_h5,
            "output_h5_path": output_h5_path,
            "output_manager": output_manager,
            "holodoppler_h5": hd_h5,
            "doppler_vision_h5": dv_h5,
            "holodoppler_config": hd_config,
            "doppler_vision_config": dv_config,
            "shared_inputs": shared_inputs,
            "scratch": scratch,
            "variables": context_vars,
            "computation_cache": computation_cache,
            "pipeline_options": pipeline_options,
            "pipeline_order": tuple(pipeline.name for pipeline in pipelines),
            "pipeline_targets": target_names,
            "velocity_estimation_method": None,
            "velocity_workflow_selection": velocity_workflow_selection,
            "execution_variant": None,
            "band_ratio_frequency_scale_hz": band_ratio_frequency_scale_hz,
            "on_pipeline_success": None,
            "on_progress": None,
            "background_artifacts": background_artifacts,
        }
        name = pipeline_desc.name
        dependencies = dag.dependencies_of(name, pipeline_options=pipeline_options)
        succeeded = False
        if name in variant_dependents and execution_variants:
            if pipeline_desc.produces_execution_variants:
                raise RuntimeError(
                    "Nested execution-variant producers are not supported."
                )
            variant_state_keys = {
                key for variant in execution_variants for key in variant.state
            }
            for variant in tuple(execution_variants):
                statuses = branch_status.setdefault(variant.name, {})
                label = f"{name} [{variant.name}]"
                blocked = [
                    dependency for dependency in dependencies
                    if statuses.get(dependency, shared_status.get(dependency)) != "succeeded"
                ]
                if pipeline_desc.wait_for_artifacts:
                    blocked.extend(background_artifacts.failures)
                if blocked:
                    statuses[name] = "skipped"
                    pipeline_skips[label] = "Unavailable dependencies: " + ", ".join(blocked)
                    variant_failures.setdefault(variant.name, pipeline_skips[label])
                    Logger.log_warning(f"Skipping {label}: {pipeline_skips[label]}")
                    continue
                # Variant state persists across descriptors; shared producers
                # are added without replacing variant-owned data.
                for key, value in context_vars.items():
                    if key not in variant_state_keys:
                        variant.state.setdefault(key, value)
                branch_arguments = {
                    **arguments,
                    "variables": variant.state,
                    "execution_variant": variant,
                    "processing_root": variant.output_namespace,
                    "output_provenance": variant.provenance,
                    "output_manager": output_manager.for_artifact_namespace(
                        variant.artifact_namespace
                    ),
                }
                try:
                    _run_pipeline_descriptor(pipeline_desc, **branch_arguments)
                except Exception as exc:
                    Logger.log_warning(f"Failed {label}; continuing safe pipelines: {exc}")
                    statuses[name] = "failed"
                    pipeline_failures[label] = str(exc)
                    variant_failures.setdefault(variant.name, str(exc))
                else:
                    statuses[name] = "succeeded"
                    succeeded = True
        else:
            # Prefer a complete workflow for unrelated shared pipelines.
            if execution_variants and not pipeline_desc.run_once:
                canonical = next(
                    (variant for variant in execution_variants
                     if variant.name not in variant_failures), execution_variants[0],
                )
                arguments["execution_variant"] = canonical
                context_vars.update(canonical.state)
            blocked = [
                dependency for dependency in dependencies
                if shared_status.get(dependency) != "succeeded"
            ]
            if pipeline_desc.wait_for_artifacts:
                blocked.extend(background_artifacts.failures)
            emitted_variants = ()
            if blocked:
                shared_status[name] = "skipped"
                pipeline_skips[name] = "Unavailable dependencies: " + ", ".join(blocked)
                Logger.log_warning(f"Skipping {name}: {pipeline_skips[name]}")
            else:
                try:
                    emitted_variants, emitted_failures = _run_pipeline_descriptor(
                        pipeline_desc, **arguments,
                    )
                except Exception as exc:
                    shared_status[name] = "failed"
                    pipeline_failures[name] = str(exc)
                    Logger.log_warning(f"Failed {name}; continuing safe pipelines: {exc}")
                else:
                    shared_status[name] = "succeeded"
                    succeeded = True
            if emitted_variants:
                if execution_variants:
                    raise RuntimeError(
                        "Multiple active execution-variant producers are not supported."
                    )
                execution_variants = list(emitted_variants)
                variant_failures.update(emitted_failures)
        background_artifacts.collect()
        _write_artifact_status(work_h5, background_artifacts)
        _write_execution_variant_status(
            work_h5,
            [variant for variant in execution_variants if variant.name not in variant_failures],
            variant_failures,
        )
        write_root_attributes(work_h5, {
            "pipeline_status": json.dumps(
                {"shared": shared_status, "variants": branch_status}, sort_keys=True,
            ),
            "pipeline_failures": json.dumps(pipeline_failures, sort_keys=True),
            "pipeline_skips": json.dumps(pipeline_skips, sort_keys=True),
        })
        work_h5.flush()
        if succeeded and on_pipeline_success is not None:
            on_pipeline_success(pipeline_desc.name)
        if on_progress is not None:
            on_progress()
    write_root_attributes(work_h5, {
        "pipeline_status": json.dumps({"shared": shared_status, "variants": branch_status}, sort_keys=True),
        "pipeline_failures": json.dumps(pipeline_failures, sort_keys=True),
        "pipeline_skips": json.dumps(pipeline_skips, sort_keys=True),
    })
    _write_execution_variant_status(work_h5,
        [variant for variant in execution_variants if variant.name not in variant_failures],
        variant_failures)
    # Reports, archives and run completion must see every completed sidecar.
    if wait_for_exports and background_artifacts.pending:
        Logger.log("Waiting for background figure exports before completing this run...")
    background_artifacts.collect(wait=wait_for_exports)
    _write_artifact_status(work_h5, background_artifacts)
    work_h5.flush()
    # Report failure only after every safe pipeline has had a chance to run.
    if (wait_for_exports and background_artifacts.failures) or "failed" in shared_status.values() or (
        execution_variants
        and all(variant.name in variant_failures for variant in execution_variants)
    ):
        errors = list(pipeline_failures.values()) + [
            f"{name}: {error}" for name, error in background_artifacts.failures.items()
        ]
        raise PipelineRunError("; ".join(errors), output_h5_path)
    return output_h5_path


def _open_input_h5_sources(
    stack: ExitStack,
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
):
    hd_h5 = (
        stack.enter_context(open_h5(holodoppler_h5, "r")) if holodoppler_h5 is not None else None
    )
    dv_h5 = (
        stack.enter_context(open_h5(doppler_vision_h5, "r"))
        if doppler_vision_h5 is not None
        else None
    )
    return hd_h5, dv_h5


def _initialize_work_h5(
    *,
    work_h5,
    pipelines: Sequence[PipelineDescriptor],
    target_names: Sequence[str],
    pipeline_options: Mapping[str, Sequence[str]],
    band_ratio_frequency_scale_hz: float,
    copy_hd_passthrough: bool,
    velocity_workflow_selection: str,
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
    hd_source_h5=None,
    dv_source_h5=None,
) -> None:
    initialize_output_h5(
        work_h5,
        holodoppler_source_file=(str(holodoppler_h5) if holodoppler_h5 is not None else None),
        doppler_vision_source_file=(
            str(doppler_vision_h5) if doppler_vision_h5 is not None else None
        ),
        holodoppler_source_h5=hd_source_h5,
        doppler_vision_source_h5=dv_source_h5,
        copy_hd_passthrough=copy_hd_passthrough,
    )
    write_root_attributes(work_h5, {
        "pipeline_targets": list(target_names),
        "pipeline_order": [pipeline.name for pipeline in pipelines],
        "velocity_estimation_methods": list(selected_velocity_methods(velocity_workflow_selection)),
        "velocity_quantity": "physical_velocity",
        "velocity_unit": "mm/s",
        "band_ratio_frequency_scale_hz": band_ratio_frequency_scale_hz,
        "copy_hd_passthrough": copy_hd_passthrough,
        "pipeline_options": json.dumps(
            {name: list(options) for name, options in pipeline_options.items()},
            sort_keys=True,
        ),
    })


def _write_artifact_status(work_h5, background_artifacts):
    if background_artifacts.statuses:
        work_h5.attrs["artifact_status"] = json.dumps(
            background_artifacts.statuses, sort_keys=True,
        )
        work_h5.attrs["artifact_failures"] = json.dumps(
            background_artifacts.failures, sort_keys=True,
        )


def _write_execution_variant_status(
    work_h5,
    variants: Sequence[ExecutionVariant],
    failures: Mapping[str, str],
) -> None:
    """Persist generic execution-variant status."""

    write_root_attributes(work_h5, {
        "execution_variant_failures": json.dumps(dict(failures), sort_keys=True),
        "execution_variants_completed": [variant.name for variant in variants],
    })
