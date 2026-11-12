"""Runtime execution helpers for resolved EyeFlow pipelines."""

import json
import multiprocessing
from collections.abc import Callable, Mapping, Sequence
from contextlib import ExitStack
from dataclasses import dataclass, field, replace
from pathlib import Path
from time import perf_counter
from uuid import uuid4

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
from pipeline_engine.data.scratch import RunArrayStorage
from input_output.writers.h5 import (
    open_h5, write_root_attributes,
)
from input_output.writers.h5_initialization import initialize_output_h5, initialize_acquisition_metadata
from input_output.writers.h5_output import PipelineH5Output
from utils.logger import Logger
from velocity_calibration import (
    DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
    validate_band_ratio_frequency_scale_hz,
)

from pipeline_engine.planning.base import ExecutionVariant, PipelineDescriptor, PipelineSelection
from pipeline_engine.execution.background_artifacts import BackgroundArtifacts
from pipeline_engine.context import PipelineInputs
from pipeline_engine.planning.dag import PipelineDAG
from pipeline_engine.execution.invocation import _run_pipeline_descriptor
from pipeline_engine.data.storage import ProductStorage
from pipeline_engine.output.output_model import OutputModel
from pipeline_engine.resources import ResourceBudget, current_resources
from pipeline_engine.data.input_service import InputService
from pipeline_engine.execution.scheduler import CoordinatorScheduler
from pipeline_engine.execution.process_executor import ProcessExecutor
from pipeline_engine.execution.protocol import InvocationSettings
from pipeline_engine.stage_graph.model import StageGraph
from pipeline_engine.output.finalization import PendingOutput
from pipeline_engine.product_exchange import PipelineInvocationId, ProductId, ProductExchangeRegistry
from pipeline_engine.product_exchange.acceptance import ProductAcceptance


def _bind_product_requirements(registry, descriptors, identity, seen):
    """Declare qualified demands, including named variants not yet produced."""
    if identity in seen:
        return
    seen.add(identity)
    descriptor = descriptors[identity.pipeline]
    for requirement in descriptor.product_requires:
        product_id = requirement.resolve(identity)
        registry.declare(product_id)
        registry.add_dependent(product_id, identity)
        _bind_product_requirements(registry, descriptors, product_id.producer, seen)


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
                scratch = stack.enter_context(
                    RunArrayStorage(output_manager.layout.ef_dir, metrics=storage.metrics)
                )
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



@dataclass(frozen=True)
class _InputSources:
    """Run-owned replay service and the shared source views it supplies."""

    service: InputService
    files: dict
    configs: dict
    shared: PipelineInputs


@dataclass(frozen=True)
class _ExecutionPlan:
    """Preflight declarations reused by input planning and variant registration."""

    dag: PipelineDAG
    descriptors: dict[str, PipelineDescriptor]
    graphs: dict[str, StageGraph | None]
    variant_dependents: set[str]
    registry: ProductExchangeRegistry
    run_id: str


@dataclass
class _ExecutionStatus:
    """One input's outcomes; the scheduler shares these dictionaries by reference."""

    execution_variants: list[ExecutionVariant] = field(default_factory=list)
    variant_failures: dict[str, str] = field(default_factory=dict)
    shared_status: dict[str, str] = field(default_factory=dict)
    branch_status: dict[str, dict[str, str]] = field(default_factory=dict)
    pipeline_failures: dict[str, str] = field(default_factory=dict)
    pipeline_skips: dict[str, str] = field(default_factory=dict)


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
    """Prepare sources and execution, run the plan, then finalize input status.

    The caller's ExitStack remains responsible for worker, client and source
    cleanup. Final output publication remains with PendingOutput.
    """
    hd_h5, dv_h5 = _open_input_h5_sources(stack, holodoppler_h5, doppler_vision_h5)
    _initialize_work_h5(
        work_h5=work_h5, pipelines=pipelines, target_names=target_names,
        pipeline_options=pipeline_options,
        band_ratio_frequency_scale_hz=band_ratio_frequency_scale_hz,
        copy_hd_passthrough=copy_hd_passthrough,
        velocity_workflow_selection=velocity_workflow_selection,
        holodoppler_h5=holodoppler_h5, doppler_vision_h5=doppler_vision_h5,
        hd_source_h5=hd_h5, dv_source_h5=dv_h5,
    )
    sources = _prepare_input_replay(
        stack, work_h5, hd_h5, dv_h5, band_limited_harmonics,
    )
    plan = _build_execution_plan(pipelines, target_names, pipeline_options)
    sources.service.configure(
        plan.graphs, {"velocity_workflow_selection": velocity_workflow_selection},
    )
    status = _ExecutionStatus()
    variables: dict[str, object] = {}
    computation_cache: dict[object, object] = {}
    settings = InvocationSettings(
        configs=sources.configs, pipeline_options=dict(pipeline_options),
        pipeline_order=tuple(pipeline.name for pipeline in pipelines),
        pipeline_targets=tuple(target_names),
        velocity_workflow_selection=velocity_workflow_selection,
        band_ratio_frequency_scale_hz=band_ratio_frequency_scale_hz,
    )
    scheduler = _create_scheduler(
        stack=stack, work_h5=work_h5, output_manager=output_manager,
        pipelines=pipelines, pipeline_options=pipeline_options,
        sources=sources, plan=plan, status=status, variables=variables,
        settings=settings, on_pipeline_start=on_pipeline_start,
        on_pipeline_success=on_pipeline_success, on_progress=on_progress,
    )
    arguments = _coordinator_arguments(
        work_h5=work_h5, output_h5_path=output_h5_path, output_manager=output_manager,
        scratch=scratch, sources=sources, settings=settings,
        pipeline_options=pipeline_options, variables=variables,
        computation_cache=computation_cache, background_artifacts=background_artifacts,
    )
    _execute_coordinator_pipelines(
        scheduler, plan, status, sources, arguments, background_artifacts,
        on_pipeline_start=on_pipeline_start, on_pipeline_success=on_pipeline_success,
        on_progress=on_progress,
    )
    _finish_input_execution(
        work_h5, output_h5_path, status, background_artifacts,
        wait_for_exports=wait_for_exports,
    )
    return output_h5_path


def _prepare_input_replay(stack, work_h5, hd_h5, dv_h5, band_limited_harmonics) -> _InputSources:
    """Register replay cleanup and build the shared typed HD/DV source views."""
    started = perf_counter()
    Logger.log_perf(
        f"START pipeline=run variant=shared stage=input_replay_setup at={started:.6f}"
    )
    service = InputService(
        {"hd": hd_h5, "dv": dv_h5}, work_h5.storage, multiprocessing.get_context("spawn"),
    )
    stack.callback(service.close)
    client = service.client()
    stack.callback(client.close)
    hd_config = load_h5_sidecar_config(hd_h5, source="hd")
    write_root_attributes(work_h5, {"band_limited_harmonics": band_limited_harmonics})
    dv_config = load_h5_sidecar_config(dv_h5, source="dv")
    files = service.files(client)
    shared = PipelineInputs(
        hd=PipelineInputSource(RawH5SourceReader(h5file=files["hd"], label="HD"), hd_config),
        dv=PipelineInputSource(RawH5SourceReader(h5file=files["dv"], label="DV"), dv_config),
    )
    initialize_acquisition_metadata(work_h5, shared.hd)
    Logger.log_perf(
        f"END pipeline=run variant=shared stage=input_replay_setup "
        f"at={perf_counter():.6f} wall={perf_counter() - started:.3f}s"
    )
    return _InputSources(service, files, {"hd": hd_config, "dv": dv_config}, shared)


def _build_execution_plan(pipelines, target_names, pipeline_options) -> _ExecutionPlan:
    """Build each stage graph once and identify variant-scoped consumers."""
    dag = PipelineDAG(pipelines)
    selection = PipelineSelection(
        tuple(item.name for item in pipelines), tuple(target_names), pipeline_options,
    )
    variant_dependents = {
        name
        for producer in pipelines if producer.produces_execution_variants
        for name in dag.dependents_of(
            producer.name, transitive=True, pipeline_options=pipeline_options,
        )
    }
    variant_dependents.difference_update(p.name for p in pipelines if p.run_once)
    run_id = uuid4().hex
    registry = ProductExchangeRegistry()
    graphs = {}
    for descriptor in pipelines:
        identity = PipelineInvocationId(run_id, descriptor.name)
        registry.register_invocation(identity)
        graph = descriptor.build_graph(
            pipeline_options.get(descriptor.name, ()), selection=selection,
        )
        graphs[descriptor.name] = graph
        _declare_invocation_products(registry, identity, graph)
    return _ExecutionPlan(
        dag, {p.name: p for p in pipelines}, graphs, variant_dependents, registry, run_id,
    )


def _declare_invocation_products(registry, identity, graph) -> None:
    """Declare exposed graph products for an already registered invocation."""
    if graph is not None:
        for stage in graph.stages:
            for product in stage.products:
                if product.exposed:
                    registry.declare(ProductId(identity, product.id))


def _create_scheduler(
    *, stack, work_h5, output_manager, pipelines, pipeline_options,
    sources, plan, status, variables, settings,
    on_pipeline_start, on_pipeline_success, on_progress,
) -> CoordinatorScheduler:
    """Wire coordinator acceptance and worker transport to the existing scheduler."""
    bound = set()
    for descriptor in pipelines:
        _bind_product_requirements(
            plan.registry, plan.descriptors, PipelineInvocationId(plan.run_id, descriptor.name), bound,
        )
    acceptance = ProductAcceptance(
        model=work_h5, state=variables, destination=output_manager.layout.ef_dir,
        inputs=sources.files, registry=plan.registry,
    )
    executor = ProcessExecutor(
        context=multiprocessing.get_context("spawn"), input_service=sources.service,
        inputs=sources.files, model=work_h5, manager=output_manager, budget=current_resources(),
        state=variables, acceptance=acceptance, registry=plan.registry,
        run_id=plan.run_id, settings=settings,
    )
    stack.callback(executor.close)
    return CoordinatorScheduler(
        pipelines=[replace(p, execution_mode="coordinator")
                   if p.name in plan.variant_dependents or p.produces_execution_variants
                   or p.wait_for_artifacts else p for p in pipelines],
        dag=plan.dag, options=pipeline_options, executor=executor,
        statuses=status.shared_status, failures=status.pipeline_failures, skips=status.pipeline_skips,
        on_start=on_pipeline_start, on_success=on_pipeline_success, on_progress=on_progress,
        registry=plan.registry,
    )


def _coordinator_arguments(
    *, work_h5, output_h5_path, output_manager, scratch, sources, settings,
    pipeline_options, variables, computation_cache, background_artifacts,
) -> dict:
    """Assemble common invocation arguments; each dispatch receives its own dict."""
    return {
        **settings.context_arguments(),
        "work_h5": work_h5,
        "output_h5_path": output_h5_path,
        "output_manager": output_manager,
        "holodoppler_h5": sources.files["hd"],
        "doppler_vision_h5": sources.files["dv"],
        "holodoppler_config": sources.configs["hd"],
        "doppler_vision_config": sources.configs["dv"],
        "shared_inputs": sources.shared,
        "scratch": scratch,
        "variables": variables,
        "computation_cache": computation_cache,
        "pipeline_options": pipeline_options,
        "on_pipeline_success": None,
        "on_progress": None,
        "background_artifacts": background_artifacts,
    }


def _execute_coordinator_pipelines(
    scheduler, plan, status, sources, common_arguments, background_artifacts,
    *, on_pipeline_start, on_pipeline_success, on_progress,
) -> None:
    """Dispatch ready coordinator jobs, checkpoint outcomes, and notify progress."""
    work_h5 = common_arguments["work_h5"]
    executor = scheduler.executor
    for pipeline_index, descriptor in scheduler.coordinator_jobs():
        # Only consumers of artifacts wait; ordinary calculations keep running.
        background_artifacts.collect(wait=descriptor.wait_for_artifacts)
        _write_artifact_status(work_h5, background_artifacts)
        if descriptor.wait_for_artifacts:
            work_h5.flush()
        if on_pipeline_start is not None:
            on_pipeline_start(descriptor.name, pipeline_index, scheduler.total)
        sources.service.activate(descriptor.name)
        arguments = {
            **common_arguments,
            "product_inputs": executor.product_inputs(descriptor)
                if descriptor.name not in plan.variant_dependents else {},
            "on_product_commit": lambda stage, product, handle, descriptor=descriptor:
                executor.acceptance.accept_local(
                    executor.invocation_id(descriptor), stage, product, handle,
                ),
        }
        dependencies = plan.dag.dependencies_of(
            descriptor.name, pipeline_options=scheduler.options, include_products=False,
        )
        if descriptor.name in plan.variant_dependents and status.execution_variants:
            succeeded = _run_variant_pipeline(
                descriptor, dependencies, arguments, status, executor, background_artifacts,
            )
        else:
            succeeded = _run_shared_pipeline(
                descriptor, dependencies, arguments, plan, status, executor, background_artifacts,
            )
        background_artifacts.collect()
        _write_artifact_status(work_h5, background_artifacts)
        _write_execution_variant_status(
            work_h5,
            [v for v in status.execution_variants if v.name not in status.variant_failures],
            status.variant_failures,
        )
        _write_pipeline_status(work_h5, status)
        work_h5.flush()
        if succeeded and on_pipeline_success is not None:
            on_pipeline_success(descriptor.name)
        if on_progress is not None:
            on_progress()


def _run_variant_pipeline(
    descriptor, dependencies, arguments, status, executor, background_artifacts,
) -> bool:
    """Run a consumer independently in each published execution variant."""
    if descriptor.produces_execution_variants:
        raise RuntimeError("Nested execution-variant producers are not supported.")
    variables = arguments["variables"]
    variant_state_keys = {key for v in status.execution_variants for key in v.state}
    succeeded = False
    for variant in tuple(status.execution_variants):
        statuses = status.branch_status.setdefault(variant.name, {})
        label = f"{descriptor.name} [{variant.name}]"
        blocked = [dependency for dependency in dependencies
                   if statuses.get(dependency, status.shared_status.get(dependency)) != "succeeded"]
        if descriptor.wait_for_artifacts:
            blocked.extend(background_artifacts.failures)
        if blocked:
            statuses[descriptor.name] = "skipped"
            status.pipeline_skips[label] = "Unavailable dependencies: " + ", ".join(blocked)
            status.variant_failures.setdefault(variant.name, status.pipeline_skips[label])
            Logger.log_warning(f"Skipping {label}: {status.pipeline_skips[label]}")
            continue
        # Keep branch bindings isolated while adding newly available shared data.
        for key, value in variables.items():
            if key not in variant_state_keys:
                variant.state.setdefault(key, value)
        branch_arguments = {
            **arguments,
            "variables": variant.state,
            "execution_variant": variant,
            "processing_root": variant.output_namespace,
            "output_provenance": variant.provenance,
            "output_manager": arguments["output_manager"].for_artifact_namespace(
                variant.artifact_namespace,
            ),
            "product_inputs": executor.product_inputs(descriptor, variant.name),
            "on_product_commit": lambda stage, product, handle, variant_name=variant.name:
                executor.acceptance.accept_local(
                    executor.invocation_id(descriptor, variant_name), stage, product, handle,
                ),
        }
        try:
            outcome = _run_pipeline_descriptor(descriptor, **branch_arguments)
        except Exception as exc:
            Logger.log_warning(f"Failed {label}; continuing safe pipelines: {exc}")
            statuses[descriptor.name] = "failed"
            status.pipeline_failures[label] = str(exc)
            status.variant_failures.setdefault(variant.name, str(exc))
        else:
            statuses[descriptor.name] = outcome.status
            succeeded = outcome.status == "succeeded"
            if outcome.reason:
                status.pipeline_skips[label] = outcome.reason
        executor.acceptance.finalize(
            executor.invocation_id(descriptor, variant.name), status.pipeline_failures.get(label),
        )
    return succeeded


def _run_shared_pipeline(
    descriptor, dependencies, arguments, plan, status, executor, background_artifacts,
) -> bool:
    """Run one shared invocation and register any newly published variants."""
    name = descriptor.name
    if status.execution_variants and not descriptor.run_once:
        # Prefer a complete workflow for unrelated shared pipelines.
        canonical = next(
            (v for v in status.execution_variants if v.name not in status.variant_failures),
            status.execution_variants[0],
        )
        arguments["execution_variant"] = canonical
        arguments["variables"].update(canonical.state)
    blocked = [dependency for dependency in dependencies
               if status.shared_status.get(dependency) != "succeeded"]
    if descriptor.wait_for_artifacts:
        blocked.extend(background_artifacts.failures)
    succeeded = False
    if blocked:
        status.shared_status[name] = "skipped"
        status.pipeline_skips[name] = "Unavailable dependencies: " + ", ".join(blocked)
        Logger.log_warning(f"Skipping {name}: {status.pipeline_skips[name]}")
    else:
        try:
            outcome = _run_pipeline_descriptor(descriptor, **arguments)
        except Exception as exc:
            status.shared_status[name] = "failed"
            status.pipeline_failures[name] = str(exc)
            Logger.log_warning(f"Failed {name}; continuing safe pipelines: {exc}")
        else:
            status.shared_status[name] = outcome.status
            succeeded = outcome.status == "succeeded"
            if outcome.reason:
                status.pipeline_skips[name] = outcome.reason
            if outcome.variants:
                _register_execution_variants(
                    outcome.variants, outcome.variant_failures, arguments["work_h5"],
                    plan, status, executor,
                )
    executor.acceptance.finalize(
        executor.invocation_id(descriptor),
        status.pipeline_failures.get(name) or status.pipeline_skips.get(name),
    )
    return succeeded


def _register_execution_variants(variants, failures, work_h5, plan, status, executor) -> None:
    """Bind variant provenance and qualify downstream product identities."""
    if status.execution_variants:
        raise RuntimeError("Multiple active execution-variant producers are not supported.")
    status.execution_variants = list(variants)
    status.variant_failures.update(failures)
    for variant in status.execution_variants:
        work_h5.designate_metadata(variant.output_namespace, "runtime")
        PipelineH5Output(work_h5, processing_root=variant.output_namespace).set_attrs(variant.provenance)
    # Graph declarations are reused; only invocation qualification changes.
    for name in plan.descriptors:
        if name not in plan.variant_dependents:
            continue
        descriptor = plan.descriptors[name]
        plan.registry.remove_placeholder(executor.invocation_id(descriptor))
        for variant in status.execution_variants:
            identity = executor.invocation_id(descriptor, variant.name)
            plan.registry.register_invocation(identity)
            _declare_invocation_products(plan.registry, identity, plan.graphs[name])
    for name in plan.descriptors:
        if name in plan.variant_dependents:
            for identity in plan.registry.invocations_for(name):
                _bind_product_requirements(plan.registry, plan.descriptors, identity, set())


def _write_pipeline_status(work_h5, status) -> None:
    """Persist shared and branch outcomes without changing their ownership."""
    write_root_attributes(work_h5, {
        "pipeline_status": json.dumps(
            {"shared": status.shared_status, "variants": status.branch_status}, sort_keys=True,
        ),
        "pipeline_failures": json.dumps(status.pipeline_failures, sort_keys=True),
        "pipeline_skips": json.dumps(status.pipeline_skips, sort_keys=True),
    })


def _finish_input_execution(
    work_h5, output_h5_path, status, background_artifacts, *, wait_for_exports,
) -> None:
    """Join required exports and report failure only after safe work completes."""
    _write_pipeline_status(work_h5, status)
    _write_execution_variant_status(
        work_h5, [v for v in status.execution_variants if v.name not in status.variant_failures],
        status.variant_failures,
    )
    if wait_for_exports and background_artifacts.pending:
        Logger.log("Waiting for background figure exports before completing this run...")
    background_artifacts.collect(wait=wait_for_exports)
    _write_artifact_status(work_h5, background_artifacts)
    work_h5.flush()
    artifacts_failed = wait_for_exports and bool(background_artifacts.failures)
    shared_failed = "failed" in status.shared_status.values()
    all_variants_failed = bool(status.execution_variants) and all(
        v.name in status.variant_failures for v in status.execution_variants
    )
    if artifacts_failed or shared_failed or all_variants_failed:
        errors = list(status.pipeline_failures.values()) + [
            f"{name}: {error}" for name, error in background_artifacts.failures.items()
        ]
        raise PipelineRunError("; ".join(errors), output_h5_path)


def _open_input_h5_sources(
    stack: ExitStack,
    holodoppler_h5: Path | None,
    doppler_vision_h5: Path | None,
):
    started = perf_counter()
    Logger.log_perf(
        f"START pipeline=run variant=shared stage=input_hdf5_open at={started:.6f}"
    )
    hd_h5 = (
        stack.enter_context(open_h5(holodoppler_h5, "r")) if holodoppler_h5 is not None else None
    )
    dv_h5 = (
        stack.enter_context(open_h5(doppler_vision_h5, "r"))
        if doppler_vision_h5 is not None
        else None
    )
    Logger.log_perf(
        f"END pipeline=run variant=shared stage=input_hdf5_open "
        f"at={perf_counter():.6f} wall={perf_counter() - started:.3f}s"
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
    started = perf_counter()
    Logger.log_perf(
        f"START pipeline=run variant=shared stage=output_initialization at={started:.6f}"
    )
    with work_h5.checkpoint():
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
    Logger.log_perf(
        f"END pipeline=run variant=shared stage=output_initialization "
        f"at={perf_counter():.6f} wall={perf_counter() - started:.3f}s"
    )


def _write_artifact_status(work_h5, background_artifacts):
    if background_artifacts.statuses:
        work_h5.attrs["artifact_status"] = json.dumps(
            background_artifacts.statuses, sort_keys=True,
        )
        work_h5.attrs["artifact_failures"] = json.dumps(
            background_artifacts.failures, sort_keys=True,
        )
        work_h5.attrs["artifact_timings"] = json.dumps(background_artifacts.timings, sort_keys=True)


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
