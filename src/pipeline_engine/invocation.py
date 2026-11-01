"""Execute one descriptor with the same checkpoints in either process mode."""

from collections.abc import Callable, Mapping, Sequence
from contextlib import contextmanager, nullcontext
from time import perf_counter

from app_settings import DEFAULT_VELOCITY_WORKFLOW_SELECTION, VelocityEstimationMethod
from input_output.artifact_ownership import (
    ArtifactOwnership, ArtifactScope, artifact_transaction, current_artifact_transaction,
)
from input_output.output_manager import OutputManager
from input_output.writers.deferred import capture_artifact_writes
from input_output.writers.h5 import h5_write_checkpoint
from utils.logger import Logger

from .background_artifacts import BackgroundArtifacts
from .base import ExecutionVariant, PipelineDescriptor, ProcessResult
from .context import PipelineContext, PipelineInputs, apply_pipeline_result, finish_pipeline
from .errors import format_pipeline_exception
from .protocol import emit
from .resources import current_resources, pipeline_memory_bytes
from .scratch import RunArrayStorage


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
    background_artifacts: BackgroundArtifacts | None = None,
    computation_cache: dict[object, object] | None = None,
    velocity_workflow_selection: str = DEFAULT_VELOCITY_WORKFLOW_SELECTION,
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
        computation_cache=computation_cache,
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
        velocity_workflow_selection=velocity_workflow_selection,
        band_ratio_frequency_scale_hz=band_ratio_frequency_scale_hz,
        background_artifacts=background_artifacts,
    )
    ctx._graph_runtime_available = True
    try:
        with _pipeline_checkpoint(work_h5, variables, output_manager, background_artifacts) as checkpoint:
            # Consumers may create intermediate images they use immediately.
            # Their upstream artifacts have already been joined by the runtime.
            capture = None if pipeline_desc.wait_for_artifacts else ctx._capture_artifact_write
            with capture_artifact_writes(capture):
                result = _execute_pipeline(pipeline_desc, pipeline, ctx, output_h5_path, checkpoint)
    except Exception as exc:
        ctx._cancel_artifacts()
        Logger.log_error(f"{pipeline.name}: {exc}")
        raise RuntimeError(format_pipeline_exception(exc, pipeline)) from exc
    except BaseException:
        ctx._cancel_artifacts()
        raise
    # Publish the invocation snapshot only after its outer checkpoint accepts
    # every legacy result. Stage graphs publish their earlier handoffs separately.
    emit("ProductCommit")
    ctx._start_artifacts()
    if ctx.artifacts_pending:
        Logger.log(f"[OK] {label}: calculations complete; artifact exports continue in background.")
    else:
        Logger.log(f"[OK] {label}")
    return result


def _execute_pipeline(pipeline_desc, pipeline, ctx, output_h5_path, checkpoint):
    """Calculate and persist a descriptor within its output checkpoint."""
    budget = current_resources()
    memory = pipeline_memory_bytes(pipeline_desc, budget) if budget else 0
    with budget.compute(memory=memory) if budget else nullcontext():
        result = pipeline.run(ctx)
    from .graphs import PipelineGraph, GraphExecution

    if isinstance(result, PipelineGraph):
        GraphExecution(result, ctx, checkpoint, _pipeline_checkpoint).run()
        result = None
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
    return ctx.execution_variants, ctx.execution_variant_failures


@contextmanager
def _pipeline_checkpoint(work_h5, variables, output_manager, background_artifacts=None):
    """Restore one invocation and remove only artifacts explicitly owned by it."""
    state = dict(variables)
    parent_owner = current_artifact_transaction()
    registry = (background_artifacts.ownership if background_artifacts else
                parent_owner._registry if parent_owner is not None else ArtifactOwnership())
    transaction = ArtifactScope(registry, output_manager.layout.ef_dir)
    try:
        with artifact_transaction(transaction), h5_write_checkpoint(work_h5) as h5_checkpoint:
            class Checkpoint:
                def commit(self, paths=()):
                    h5_checkpoint.commit(paths)
                    state.clear()
                    state.update(variables)
                    transaction.commit()

            yield Checkpoint()
    except BaseException:
        variables.clear()
        variables.update(state)
        transaction.rollback()
        raise
