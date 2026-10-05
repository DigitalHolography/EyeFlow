"""Shared GUI/CLI orchestration for EyeFlow pipeline runs."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from input_output import HoloRunLayout, resolve_selected_run_layouts
from input_output.output_manager import OutputManager, OutputType
from input_output.run_selection import (
    ExpandedRunInputs,
    batch_root as infer_batch_root,
    expand_run_inputs,
    output_manager_for_layout,
    reject_duplicate_destinations,
)
from utils.logger import Logger

from .base import PipelineDescriptor
from .dag import PipelineDAG, PipelineExecutionPlan
from .runtime import run_pipelines_to_output

@dataclass(frozen=True)
class RunRequest:
    """One resolved HOLO input and its final output destination."""

    input_layout: HoloRunLayout
    output_manager: OutputManager


@dataclass(frozen=True)
class RunSpec:
    """Fully resolved, front-end-independent batch run specification."""

    plan: PipelineExecutionPlan
    requests: tuple[RunRequest, ...]
    pipeline_options: Mapping[str, tuple[str, ...]]

    @property
    def total_pipeline_units(self) -> int:
        return len(self.plan.descriptors) * len(self.requests)


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

    return {
        pipeline.name: pipeline
        for pipeline in pipelines
        if pipeline.visibility != "hidden"
    }


def resolve_run_spec(
    *,
    input_paths: Sequence[Path],
    target_names: Sequence[str],
    pipelines: Iterable[PipelineDescriptor],
    pipeline_options: Mapping[str, Iterable[str]] | None = None,
    output_root: Path | None = None,
    batch_root: Path | None = None,
) -> RunSpec:
    """Resolve targets, inputs, and deterministic output destinations."""

    descriptors = tuple(pipelines)
    selectable = selectable_pipeline_registry(descriptors)
    hidden_targets = [name for name in target_names if name not in selectable]
    if hidden_targets:
        raise ValueError(
            "Unknown or hidden pipeline target(s): " + ", ".join(hidden_targets)
        )

    resolved_options = _resolve_pipeline_options(descriptors, pipeline_options)
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
        raise ValueError(
            "The DAG requires unavailable pipeline(s): " + ", ".join(details)
        )
    resolved_options = {
        descriptor.name: resolved_options[descriptor.name]
        for descriptor in plan.descriptors
        if descriptor.name in resolved_options
    }

    layouts = resolve_selected_run_layouts(input_paths)
    resolved_output_root = (
        output_root.expanduser().resolve() if output_root is not None else None
    )
    effective_batch_root = (
        batch_root.expanduser().resolve()
        if batch_root is not None
        else infer_batch_root([layout.holo_path for layout in layouts])
    )
    requests = tuple(
        RunRequest(
            input_layout=layout,
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

    Logger.log(f"[DAG] Targets -> {', '.join(spec.plan.targets)}")
    Logger.log(f"[DAG] Execution order -> {', '.join(spec.plan.names)}")
    outputs: list[Path] = []
    failures: list[RunFailure] = []
    stopped = False

    for request_index, request in enumerate(spec.requests):
        input_layout = request.input_layout
        final_manager = request.output_manager
        final_path = final_manager.path_for(OutputType.H5)
        if on_file_start is not None:
            on_file_start(
                input_layout.holo_path,
                request_index + 1,
                len(spec.requests),
            )
        Logger.log(f"[INPUT] HOLO -> {input_layout.holo_path}")
        Logger.log(f"[INPUT] DATA DIR -> {input_layout.root_dir}")
        Logger.log(f"[RESOLVED] HD -> {input_layout.hd_h5}")
        Logger.log(f"[RESOLVED] DV -> {input_layout.dv_h5}")
        Logger.log(f"[OUTPUT] {final_path}")

        try:
            final_dir = final_manager.layout.ef_dir
            if final_dir.exists() and not final_dir.is_dir():
                raise RuntimeError(
                    f"Refusing to replace non-directory EyeFlow output path: {final_dir}"
                )
            # Start clean so a failed run leaves only its own partial output,
            # rather than mixing new files with artifacts from an old run.
            final_manager.prepare(replace=True)
            run_pipelines_to_output(
                output_manager=final_manager,
                pipelines=spec.plan.descriptors,
                target_names=spec.plan.targets,
                pipeline_options=spec.pipeline_options,
                holodoppler_h5=input_layout.hd_h5,
                doppler_vision_h5=input_layout.dv_h5,
                on_pipeline_start=on_pipeline_start,
                on_progress=on_progress,
            )
        except Exception as exc:  # noqa: BLE001
            failure = RunFailure(input_layout.holo_path, str(exc))
            failures.append(failure)
            Logger.log_error(str(failure))
        else:
            outputs.append(final_path)
            Logger.log(f"Completed run for {input_layout.holo_path.name}: {final_path}")

        remaining_count = len(spec.requests) - request_index - 1
        if remaining_count and should_stop is not None and should_stop():
            stopped = True
            Logger.log(
                f"[STOP] Batch stopped after {input_layout.holo_path.name}; "
                f"{remaining_count} input file(s) not started."
            )
            break

    return RunResult(tuple(outputs), tuple(failures), stopped=stopped)


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
            selected = {
                option.name
                for option in descriptor.options
                if option.default_enabled
            }
        else:
            selected = {str(name).strip() for name in requested if str(name).strip()}
            unknown = sorted(selected - known)
            if unknown:
                raise ValueError(
                    f"Unknown option(s) for pipeline '{descriptor.name}': "
                    + ", ".join(unknown)
                )
        options_by_name = {
            option.name: option for option in descriptor.options
        }
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
