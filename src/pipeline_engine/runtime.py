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

from .base import ExecutionVariant, PipelineDescriptor, ProcessResult
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
    execution_variants: list[ExecutionVariant] = []
    variant_dependents: set[str] = set()
    variant_failures: dict[str, str] = {}

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
            "execution_variant": None,
            "band_ratio_frequency_scale_hz": band_ratio_frequency_scale_hz,
            "on_pipeline_success": None,
            "on_progress": None,
        }
        if pipeline_desc.name in variant_dependents:
            if pipeline_desc.produces_execution_variants:
                raise RuntimeError(
                    "Nested execution-variant producers are not supported."
                )
            variant_state_keys = {
                key for variant in execution_variants for key in variant.state
            }
            for variant in tuple(execution_variants):
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
                except RuntimeError as exc:
                    Logger.log_warning(
                        f"Skipping {variant.name} execution variant: {exc}"
                    )
                    variant_failures[variant.name] = str(exc)
                    execution_variants.remove(variant)
                    _discard_execution_variant(work_h5, output_manager, variant)
            if not execution_variants:
                _write_execution_variant_status(
                    work_h5,
                    execution_variants,
                    variant_failures,
                )
                raise RuntimeError(
                    "No execution variant completed downstream analysis."
                )
        else:
            # Shared consumers use the first surviving variant as canonical.
            if execution_variants:
                canonical = execution_variants[0]
                arguments["execution_variant"] = canonical
                context_vars.update(canonical.state)
            emitted_variants, emitted_failures = _run_pipeline_descriptor(
                pipeline_desc,
                **arguments,
            )
            if emitted_variants:
                if execution_variants:
                    raise RuntimeError(
                        "Multiple active execution-variant producers are not supported."
                    )
                execution_variants = list(emitted_variants)
                variant_failures = emitted_failures
                variant_dependents = set(
                    dag.dependents_of(
                        pipeline_desc.name,
                        transitive=True,
                        pipeline_options=pipeline_options,
                    )
                )
        _write_execution_variant_status(
            work_h5,
            execution_variants,
            variant_failures,
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
    execution_variant: ExecutionVariant | None,
    band_ratio_frequency_scale_hz: float,
    on_pipeline_success: Callable[[str], None] | None,
    on_progress: Callable[[], None] | None,
    processing_root: str | None = None,
    output_provenance: Mapping[str, object] | None = None,
) -> tuple[tuple[ExecutionVariant, ...], dict[str, str]]:
    pipeline = pipeline_desc.instantiate()
    label = (
        f"{pipeline.name} [{execution_variant.name}]"
        if execution_variant is not None and processing_root
        else pipeline.name
    )
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
        output_provenance=output_provenance,
        execution_variant=execution_variant,
        variables=variables,
        pipeline_options=pipeline_options,
        pipeline_order=pipeline_order,
        pipeline_targets=pipeline_targets,
        velocity_estimation_method=velocity_estimation_method,
        band_ratio_frequency_scale_hz=band_ratio_frequency_scale_hz,
    )
    try:
        result = pipeline.run(ctx)
        if ctx.execution_variants and not pipeline_desc.produces_execution_variants:
            raise ValueError(
                f"Pipeline '{pipeline.name}' published execution variants without "
                "declaring produces_execution_variants=True."
            )
        if pipeline_desc.produces_execution_variants and not ctx.execution_variants:
            raise ValueError(
                f"Pipeline '{pipeline.name}' declared execution variants but "
                "published none."
            )
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
    return ctx.execution_variants, ctx.execution_variant_failures


def _discard_execution_variant(work_h5, output_manager, variant: ExecutionVariant):
    """Remove only one failed variant's output namespace and artifacts."""
    import shutil

    root = variant.output_namespace
    if root in work_h5:
        del work_h5[root]
    manager = output_manager.for_artifact_namespace(variant.artifact_namespace)
    workspace = output_manager.layout.ef_dir.resolve()
    for kind in OutputType:
        directory = manager.dir_for(kind).resolve()
        if not directory.is_relative_to(workspace):
            raise ValueError(f"Workflow directory is outside output root: {directory}")
        if directory.exists():
            shutil.rmtree(directory)


def _write_execution_variant_status(
    work_h5,
    variants: Sequence[ExecutionVariant],
    failures: Mapping[str, str],
) -> None:
    """Persist generic execution-variant status."""

    work_h5.attrs["execution_variant_failures"] = json.dumps(
        dict(failures),
        sort_keys=True,
    )
    work_h5.attrs["execution_variants_completed"] = [
        variant.name for variant in variants
    ]
