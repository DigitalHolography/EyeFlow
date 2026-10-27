"""Runtime execution helpers for resolved EyeFlow pipelines."""

import json
from collections.abc import Callable, Mapping, Sequence
from contextlib import ExitStack, contextmanager
from pathlib import Path
from time import perf_counter

from app_settings import (
    DEFAULT_BAND_LIMITED_HARMONICS_SETTING,
    VelocityEstimationMethod,
    validate_band_limited_harmonics_setting,
)
from input_output.h5_access import PipelineInputSource, RawH5SourceReader
from input_output.inputs import load_h5_sidecar_config
from input_output.output_manager import OutputManager, OutputType
from input_output.run_scratch import RunArrayStorage
from input_output.schema.eyeflow_output import VELOCITY_WORKFLOW_ROOTS
from input_output.writers.h5 import (
    h5_write_checkpoint, initialize_output_h5, open_h5, write_root_attributes,
)
from utils.logger import Logger
from velocity_calibration import (
    DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
    validate_band_ratio_frequency_scale_hz,
)

from .base import ExecutionVariant, PipelineDescriptor, ProcessResult
from .context import PipelineContext, PipelineInputs, apply_pipeline_result, finish_pipeline
from .dag import PipelineDAG
from .errors import format_pipeline_exception


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
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
    on_pipeline_start: Callable[[str, int, int], None] | None = None,
    on_pipeline_success: Callable[[str], None] | None = None,
    on_progress: Callable[[], None] | None = None,
) -> Path:
    """Run resolved pipelines and write outputs through an OutputManager."""
    resolved_harmonics = validate_band_limited_harmonics_setting(band_limited_harmonics)

    resolved_band_ratio_scale_hz = validate_band_ratio_frequency_scale_hz(
        band_ratio_frequency_scale_hz
    )
    output_manager.prepare()
    output_h5_path = output_manager.path_for(OutputType.H5)
    with ExitStack() as stack:
        work_h5 = stack.enter_context(output_manager.open_h5())
        scratch = stack.enter_context(RunArrayStorage(output_manager.layout.ef_dir))
        return _run_pipelines_with_work_h5(
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
    scratch: RunArrayStorage,
    pipelines: Sequence[PipelineDescriptor],
    target_names: Sequence[str],
    pipeline_options: Mapping[str, Sequence[str]],
    band_ratio_frequency_scale_hz: float,
    band_limited_harmonics: int,
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
        hd_source_h5=hd_h5,
        dv_source_h5=dv_h5,
    )
    hd_config = load_h5_sidecar_config(hd_h5, source="hd")
    write_root_attributes(work_h5, {"band_limited_harmonics": band_limited_harmonics})
    dv_config = load_h5_sidecar_config(dv_h5, source="dv")
    shared_inputs = PipelineInputs(
        hd=PipelineInputSource(RawH5SourceReader(h5file=hd_h5, label="HD"), hd_config),
        dv=PipelineInputSource(RawH5SourceReader(h5file=dv_h5, label="DV"), dv_config),
    )
    context_vars: dict[str, object] = {}

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
            "shared_inputs": shared_inputs,
            "scratch": scratch,
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
    # Report failure only after every safe pipeline has had a chance to run.
    if "failed" in shared_status.values() or (
        execution_variants
        and all(variant.name in variant_failures for variant in execution_variants)
    ):
        raise PipelineRunError("; ".join(pipeline_failures.values()), output_h5_path)
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
    )
    write_root_attributes(work_h5, {
        "pipeline_targets": list(target_names),
        "pipeline_order": [pipeline.name for pipeline in pipelines],
        "velocity_estimation_methods": list(VELOCITY_WORKFLOW_ROOTS),
        "velocity_quantity": "physical_velocity",
        "velocity_unit": "mm/s",
        "band_ratio_frequency_scale_hz": band_ratio_frequency_scale_hz,
        "pipeline_options": json.dumps(
            {name: list(options) for name, options in pipeline_options.items()},
            sort_keys=True,
        ),
    })


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
    shared_inputs: PipelineInputs | None = None,
    scratch: RunArrayStorage | None = None,
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
        shared_inputs=shared_inputs,
        scratch=scratch,
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
        with _pipeline_checkpoint(work_h5, variables, output_manager, scratch):
            return _execute_pipeline(pipeline_desc, pipeline, ctx, output_h5_path)
    except Exception as exc:
        Logger.log_error(f"{pipeline.name}: {exc}")
        raise RuntimeError(format_pipeline_exception(exc, pipeline)) from exc


def _execute_pipeline(pipeline_desc, pipeline, ctx, output_h5_path):
    """Calculate and persist a descriptor within its output checkpoint."""
    label = (
        f"{pipeline.name} [{ctx.execution_variant.name}]"
        if ctx.execution_variant is not None and ctx.output.h5.processing_root
        else pipeline.name
    )
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
    if isinstance(result, ProcessResult):
        result.output_h5_path = str(output_h5_path)
    finish_pipeline(ctx, pipeline.name)
    Logger.log(f"[OK] {label}")
    return ctx.execution_variants, ctx.execution_variant_failures


@contextmanager
def _pipeline_checkpoint(work_h5, variables, output_manager, scratch):
    """Restore state and remove only newly created artifacts on failure."""
    state = dict(variables)
    root = output_manager.layout.ef_dir
    artifacts = {path for path in root.rglob("*") if path.is_file()}
    try:
        with h5_write_checkpoint(work_h5):
            yield
    except Exception:
        variables.clear()
        variables.update(state)
        scratch_dir = scratch.scratch_directory if scratch is not None else None
        for path in root.rglob("*"):
            if (
                path.is_file() and path not in artifacts
                and (scratch_dir is None or scratch_dir not in path.parents)
            ):
                path.unlink()
        raise


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
