"""Runtime execution helpers for resolved EyeFlow pipelines."""

import json
from collections.abc import Callable, Mapping, Sequence
from contextlib import ExitStack
from pathlib import Path
from time import perf_counter

from app_settings import VelocityEstimationMethod
from input_output.inputs import load_h5_sidecar_config
from input_output.output_manager import OutputManager, OutputType
from input_output.schema.eyeflow_output import VELOCITY_WORKFLOW_ROOTS
from input_output.writers.h5 import initialize_output_h5, open_h5
from utils.logger import Logger
from velocity_calibration import (
    DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
    validate_band_ratio_frequency_scale_hz,
)

from .base import PipelineDescriptor, ProcessResult
from .context import PipelineContext, apply_pipeline_result, finish_pipeline
from .dag import PipelineDAG
from .errors import format_pipeline_exception


def run_pipelines_to_output(
    *,
    output_manager: OutputManager,
    pipelines: Sequence[PipelineDescriptor],
    target_names: Sequence[str] = (),
    pipeline_options: Mapping[str, Sequence[str]] | None = None,
    band_ratio_frequency_scale_hz: float = (DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ),
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
    on_pipeline_start: Callable[[str, int, int], None] | None = None,
    on_pipeline_success: Callable[[str], None] | None = None,
    on_progress: Callable[[], None] | None = None,
) -> Path:
    """Run resolved pipelines and write outputs through an OutputManager."""

    resolved_band_ratio_scale_hz = validate_band_ratio_frequency_scale_hz(
        band_ratio_frequency_scale_hz
    )
    output_manager.prepare()
    output_h5_path = output_manager.path_for(OutputType.H5)
    with ExitStack() as stack:
        work_h5 = stack.enter_context(output_manager.open_h5())
        return _run_pipelines_with_work_h5(
            work_h5=work_h5,
            output_h5_path=output_h5_path,
            output_manager=output_manager,
            stack=stack,
            pipelines=pipelines,
            target_names=target_names,
            pipeline_options=pipeline_options or {},
            band_ratio_frequency_scale_hz=resolved_band_ratio_scale_hz,
            holodoppler_h5=holodoppler_h5,
            doppler_vision_h5=doppler_vision_h5,
            on_pipeline_start=on_pipeline_start,
            on_pipeline_success=on_pipeline_success,
            on_progress=on_progress,
        )


def _run_pipelines_with_work_h5(
    *,
    work_h5,
    output_h5_path,
    output_manager: OutputManager,
    stack: ExitStack,
    pipelines: Sequence[PipelineDescriptor],
    target_names: Sequence[str],
    pipeline_options: Mapping[str, Sequence[str]],
    band_ratio_frequency_scale_hz: float,
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
    on_pipeline_start: Callable[[str, int, int], None] | None,
    on_pipeline_success: Callable[[str], None] | None,
    on_progress: Callable[[], None] | None,
) -> Path:
    hd_h5, dv_h5 = _open_input_h5_sources(stack, holodoppler_h5, doppler_vision_h5)
    _initialize_work_h5(
        work_h5=work_h5,
        pipelines=pipelines,
        target_names=target_names,
        pipeline_options=pipeline_options,
        band_ratio_frequency_scale_hz=band_ratio_frequency_scale_hz,
        holodoppler_h5=holodoppler_h5,
        doppler_vision_h5=doppler_vision_h5,
    )
    hd_config = load_h5_sidecar_config(hd_h5, source="hd")
    dv_config = load_h5_sidecar_config(dv_h5, source="dv")
    context_vars: dict[str, object] = {}

    dag = PipelineDAG(pipelines)
    branched = (
        set(dag.dependents_of("velocity", transitive=True, pipeline_options=pipeline_options))
        if any(item.name == "velocity" for item in pipelines)
        else set()
    )

    pipeline_count = len(pipelines)
    for pipeline_index, pipeline_desc in enumerate(pipelines, start=1):
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
            "variables": context_vars,
            "pipeline_options": pipeline_options,
            "pipeline_order": tuple(pipeline.name for pipeline in pipelines),
            "pipeline_targets": target_names,
            "velocity_estimation_method": None,
            "band_ratio_frequency_scale_hz": band_ratio_frequency_scale_hz,
            "on_pipeline_success": None,
            "on_progress": None,
        }
        workflows = context_vars.get("velocity_workflows", {})
        if pipeline_desc.name in branched:
            for method, state in tuple(workflows.items()):
                # Method state persists across descriptors; shared producers
                # are added without replacing previously computed method data.
                for key, value in context_vars.items():
                    if key not in {"velocity", "velocity_source", "velocity_workflows"}:
                        state.setdefault(key, value)
                branch_arguments = {
                    **arguments,
                    "variables": state,
                    "velocity_estimation_method": method,
                    "processing_root": VELOCITY_WORKFLOW_ROOTS[method],
                    "velocity_provenance": state["velocity"].provenance,
                    "output_manager": output_manager.for_workflow(method),
                }
                try:
                    _run_pipeline_descriptor(pipeline_desc, **branch_arguments)
                except RuntimeError as exc:
                    Logger.log_warning(f"Skipping {method} workflow: {exc}")
                    context_vars.setdefault("velocity_workflow_failures", {})[method] = str(exc)
                    del workflows[method]
                    _discard_workflow(work_h5, output_manager, method)
            if not workflows:
                work_h5.attrs["velocity_workflow_failures"] = json.dumps(
                    context_vars.get("velocity_workflow_failures", {}),
                    sort_keys=True,
                )
                work_h5.attrs["velocity_workflows_completed"] = []
                raise RuntimeError("Neither velocity workflow completed downstream analysis.")
        else:
            # Shared consumers use the canonical surviving workflow's cycles
            # and geometry.
            if workflows:
                canonical = next(iter(workflows.values()))
                arguments["velocity_estimation_method"] = canonical["velocity"].provenance[
                    "velocity_estimation_method"
                ]
                context_vars.update(
                    {key: canonical[key] for key in ("velocity", "velocity_source")}
                )
            _run_pipeline_descriptor(pipeline_desc, **arguments)
        work_h5.attrs["velocity_workflow_failures"] = json.dumps(
            context_vars.get("velocity_workflow_failures", {}),
            sort_keys=True,
        )
        work_h5.attrs["velocity_workflows_completed"] = list(
            context_vars.get("velocity_workflows", {})
        )
        if on_pipeline_success is not None:
            on_pipeline_success(pipeline_desc.name)
        if on_progress is not None:
            on_progress()
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
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
) -> None:
    initialize_output_h5(
        work_h5,
        holodoppler_source_file=(str(holodoppler_h5) if holodoppler_h5 is not None else None),
        doppler_vision_source_file=(
            str(doppler_vision_h5) if doppler_vision_h5 is not None else None
        ),
    )
    work_h5.attrs["pipeline_targets"] = list(target_names)
    work_h5.attrs["pipeline_order"] = [pipeline.name for pipeline in pipelines]
    work_h5.attrs["velocity_estimation_methods"] = list(VELOCITY_WORKFLOW_ROOTS)
    work_h5.attrs["velocity_quantity"] = "physical_velocity"
    work_h5.attrs["velocity_unit"] = "mm/s"
    work_h5.attrs["band_ratio_frequency_scale_hz"] = band_ratio_frequency_scale_hz
    work_h5.attrs["pipeline_options"] = json.dumps(
        {name: list(options) for name, options in pipeline_options.items()},
        sort_keys=True,
    )


def _run_pipeline_descriptor(
    pipeline_desc: PipelineDescriptor,
    *,
    work_h5,
    output_h5_path,
    output_manager: OutputManager,
    holodoppler_h5,
    doppler_vision_h5,
    holodoppler_config,
    doppler_vision_config,
    variables: dict[str, object],
    pipeline_options: Mapping[str, Sequence[str]],
    pipeline_order: Sequence[str],
    pipeline_targets: Sequence[str],
    velocity_estimation_method: VelocityEstimationMethod | None,
    band_ratio_frequency_scale_hz: float,
    on_pipeline_success: Callable[[str], None] | None,
    on_progress: Callable[[], None] | None,
    processing_root: str | None = None,
    velocity_provenance: Mapping[str, object] | None = None,
) -> None:
    pipeline = pipeline_desc.instantiate()
    label = f"{pipeline.name} [{velocity_estimation_method}]" if processing_root else pipeline.name
    Logger.log(f"[START] {label}")
    ctx = PipelineContext(
        work_h5=work_h5,
        holodoppler_h5=holodoppler_h5,
        doppler_vision_h5=doppler_vision_h5,
        holodoppler_config=holodoppler_config,
        doppler_vision_config=doppler_vision_config,
        preferred_input=pipeline_desc.input_slot,
        output_manager=output_manager,
        pipeline_name=pipeline.name,
        processing_root=processing_root,
        velocity_provenance=velocity_provenance,
        variables=variables,
        pipeline_options=pipeline_options,
        pipeline_order=pipeline_order,
        pipeline_targets=pipeline_targets,
        velocity_estimation_method=velocity_estimation_method,
        band_ratio_frequency_scale_hz=band_ratio_frequency_scale_hz,
    )
    try:
        result = pipeline.run(ctx)
        if result is not None:
            write_started = perf_counter()
            Logger.log(f"Starting {pipeline.name} HDF5 output write...")
            apply_pipeline_result(ctx, result)
            Logger.log(
                f"Completed {pipeline.name} HDF5 output write in "
                f"{perf_counter() - write_started:.1f}s."
            )
    except Exception as exc:
        Logger.log_error(f"{pipeline.name}: {exc}")
        raise RuntimeError(format_pipeline_exception(exc, pipeline)) from exc
    if isinstance(result, ProcessResult):
        result.output_h5_path = str(output_h5_path)
    finish_pipeline(ctx, pipeline.name)
    Logger.log(f"[OK] {label}")
    if on_pipeline_success is not None:
        on_pipeline_success(pipeline.name)
    if on_progress is not None:
        on_progress()


def _discard_workflow(work_h5, output_manager, method):
    """Remove only the failed workflow's owned processing group and artifacts."""
    import shutil

    root = VELOCITY_WORKFLOW_ROOTS[method]
    if root in work_h5:
        del work_h5[root]
    manager = output_manager.for_workflow(method)
    workspace = output_manager.layout.ef_dir.resolve()
    for kind in OutputType:
        directory = manager.dir_for(kind).resolve()
        if not directory.is_relative_to(workspace):
            raise ValueError(f"Workflow directory is outside output root: {directory}")
        if directory.exists():
            shutil.rmtree(directory)
