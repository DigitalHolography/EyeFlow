"""Runtime execution helpers for resolved EyeFlow pipelines."""

import json
from collections.abc import Callable, Mapping, Sequence
from contextlib import ExitStack
from pathlib import Path
from time import perf_counter

from input_output.inputs import load_h5_sidecar_config
from input_output.output_manager import OutputManager, OutputType
from input_output.writers.h5 import initialize_output_h5, open_h5
from utils.logger import Logger

from .base import PipelineDescriptor, ProcessResult
from .context import PipelineContext, apply_pipeline_result, finish_pipeline
from .errors import format_pipeline_exception


def run_pipelines_to_output(
    *,
    output_manager: OutputManager,
    pipelines: Sequence[PipelineDescriptor],
    target_names: Sequence[str] = (),
    pipeline_options: Mapping[str, Sequence[str]] | None = None,
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
    on_pipeline_start: Callable[[str, int, int], None] | None = None,
    on_pipeline_success: Callable[[str], None] | None = None,
    on_progress: Callable[[], None] | None = None,
) -> Path:
    """Run resolved pipelines and write outputs through an OutputManager."""

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
        holodoppler_h5=holodoppler_h5,
        doppler_vision_h5=doppler_vision_h5,
    )
    hd_config = load_h5_sidecar_config(hd_h5, source="hd")
    dv_config = load_h5_sidecar_config(dv_h5, source="dv")
    context_vars: dict[str, object] = {}

    pipeline_count = len(pipelines)
    for pipeline_index, pipeline_desc in enumerate(pipelines, start=1):
        if on_pipeline_start is not None:
            on_pipeline_start(
                pipeline_desc.name,
                pipeline_index,
                pipeline_count,
            )
        _run_pipeline_descriptor(
            pipeline_desc,
            work_h5=work_h5,
            output_h5_path=output_h5_path,
            output_manager=output_manager,
            holodoppler_h5=hd_h5,
            doppler_vision_h5=dv_h5,
            holodoppler_config=hd_config,
            doppler_vision_config=dv_config,
            variables=context_vars,
            pipeline_options=pipeline_options,
            pipeline_order=tuple(pipeline.name for pipeline in pipelines),
            on_pipeline_success=on_pipeline_success,
            on_progress=on_progress,
        )
    return output_h5_path


def _open_input_h5_sources(
    stack: ExitStack,
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
):
    hd_h5 = (
        stack.enter_context(open_h5(holodoppler_h5, "r"))
        if holodoppler_h5 is not None
        else None
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
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
) -> None:
    initialize_output_h5(
        work_h5,
        holodoppler_source_file=(
            str(holodoppler_h5) if holodoppler_h5 is not None else None
        ),
        doppler_vision_source_file=(
            str(doppler_vision_h5) if doppler_vision_h5 is not None else None
        ),
    )
    work_h5.attrs["pipeline_targets"] = list(target_names)
    work_h5.attrs["pipeline_order"] = [pipeline.name for pipeline in pipelines]
    work_h5.attrs["pipeline_options"] = json.dumps(
        {
            name: list(options)
            for name, options in pipeline_options.items()
        },
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
    on_pipeline_success: Callable[[str], None] | None,
    on_progress: Callable[[], None] | None,
) -> None:
    pipeline = pipeline_desc.instantiate()
    Logger.log(f"[START] {pipeline.name}")
    ctx = PipelineContext(
        work_h5=work_h5,
        holodoppler_h5=holodoppler_h5,
        doppler_vision_h5=doppler_vision_h5,
        holodoppler_config=holodoppler_config,
        doppler_vision_config=doppler_vision_config,
        preferred_input=pipeline_desc.input_slot,
        output_manager=output_manager,
        pipeline_name=pipeline.name,
        variables=variables,
        pipeline_options=pipeline_options,
        pipeline_order=pipeline_order,
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
    Logger.log(f"[OK] {pipeline.name}")
    if on_pipeline_success is not None:
        on_pipeline_success(pipeline.name)
    if on_progress is not None:
        on_progress()
