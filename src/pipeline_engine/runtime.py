from collections.abc import Callable, Sequence
from contextlib import ExitStack
from pathlib import Path

from input_output.inputs import load_h5_sidecar_config
from input_output.output_manager import OutputManager, OutputType
from input_output.writers.h5 import initialize_output_h5, open_h5

from .base import PipelineDescriptor, ProcessResult
from .context import PipelineContext
from .errors import format_pipeline_exception


def run_pipelines_to_output(
    *,
    output_manager: OutputManager,
    pipelines: Sequence[PipelineDescriptor],
    target_names: Sequence[str] = (),
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
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
            holodoppler_h5=holodoppler_h5,
            doppler_vision_h5=doppler_vision_h5,
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
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
    on_pipeline_success: Callable[[str], None] | None,
    on_progress: Callable[[], None] | None,
) -> Path:
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

    initialize_output_h5(
        work_h5,
        holodoppler_source_file=(
            str(holodoppler_h5) if holodoppler_h5 is not None else None
        ),
        doppler_vision_source_file=(
            str(doppler_vision_h5) if doppler_vision_h5 is not None else None
        ),
    )
    work_h5.attrs["trim_h5source"] = True
    work_h5.attrs["pipeline_targets"] = list(target_names)
    work_h5.attrs["pipeline_order"] = [pipeline.name for pipeline in pipelines]

    hd_config = load_h5_sidecar_config(hd_h5, source="hd")
    dv_config = load_h5_sidecar_config(dv_h5, source="dv")
    context_vars: dict[str, object] = {}

    for pipeline_desc in pipelines:
        pipeline = pipeline_desc.instantiate()
        ctx = PipelineContext(
            work_h5=work_h5,
            holodoppler_h5=hd_h5,
            doppler_vision_h5=dv_h5,
            holodoppler_config=hd_config,
            doppler_vision_config=dv_config,
            preferred_input=pipeline_desc.input_slot,
            output_manager=output_manager,
            pipeline_name=pipeline.name,
            variables=context_vars,
        )
        try:
            result = pipeline.run(ctx)
            ctx.apply_result(result)
        except Exception as exc:  # noqa: BLE001
            raise RuntimeError(format_pipeline_exception(exc, pipeline)) from exc
        if isinstance(result, ProcessResult):
            result.output_h5_path = str(output_h5_path)
        ctx.finish_pipeline(pipeline.name)
        if on_pipeline_success is not None:
            on_pipeline_success(pipeline.name)
        if on_progress is not None:
            on_progress()
    return output_h5_path
