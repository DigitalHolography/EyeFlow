"""Shared GUI/CLI orchestration for EyeFlow pipeline runs."""

from __future__ import annotations

import os
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path

from input_output import INPUT_LIST_SUFFIX, HoloRunLayout, resolve_selected_run_layouts
from input_output.archives import extracted_zip_tree
from input_output.output_manager import OutputManager, OutputType
from utils.logger import Logger

from .base import PipelineDescriptor
from .dag import PipelineDAG, PipelineExecutionPlan
from .runtime import run_pipelines_to_output

HOLO_SUFFIX = ".holo"


@dataclass(frozen=True)
class ExpandedRunInputs:
    """Input paths produced from one CLI source selection."""

    paths: tuple[Path, ...]
    batch_root: Path | None = None


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

    @property
    def succeeded(self) -> bool:
        return not self.failures

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

    plan = PipelineDAG(descriptors).resolve_targets(target_names)
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
    resolved_options = _resolve_pipeline_options(plan, pipeline_options)

    layouts = resolve_selected_run_layouts(input_paths)
    resolved_output_root = (
        output_root.expanduser().resolve() if output_root is not None else None
    )
    effective_batch_root = (
        batch_root.expanduser().resolve()
        if batch_root is not None
        else _batch_root([layout.holo_path for layout in layouts])
    )
    requests = tuple(
        RunRequest(
            input_layout=layout,
            output_manager=_output_manager_for_layout(
                layout,
                output_root=resolved_output_root,
                batch_root=effective_batch_root,
            ),
        )
        for layout in layouts
    )
    _reject_duplicate_destinations(requests)
    return RunSpec(
        plan=plan,
        requests=requests,
        pipeline_options=resolved_options,
    )


def execute_run(
    spec: RunSpec,
    *,
    on_progress: Callable[[], None] | None = None,
) -> RunResult:
    """Execute a batch, writing each run directly to its final output folder."""

    Logger.log(f"[DAG] Targets -> {', '.join(spec.plan.targets)}")
    Logger.log(f"[DAG] Execution order -> {', '.join(spec.plan.names)}")
    outputs: list[Path] = []
    failures: list[RunFailure] = []

    for request in spec.requests:
        input_layout = request.input_layout
        final_manager = request.output_manager
        final_path = final_manager.path_for(OutputType.H5)
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
                on_progress=on_progress,
            )
        except Exception as exc:  # noqa: BLE001
            failure = RunFailure(input_layout.holo_path, str(exc))
            failures.append(failure)
            Logger.log_error(str(failure))
            continue

        outputs.append(final_path)
        Logger.log(f"Completed run for {input_layout.holo_path.name}: {final_path}")

    return RunResult(tuple(outputs), tuple(failures))


def _resolve_pipeline_options(
    plan: PipelineExecutionPlan,
    selections: Mapping[str, Iterable[str]] | None,
) -> dict[str, tuple[str, ...]]:
    requested_by_pipeline = selections or {}
    resolved: dict[str, tuple[str, ...]] = {}
    for descriptor in plan.descriptors:
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
        resolved[descriptor.name] = tuple(
            option.name for option in descriptor.options if option.name in selected
        )
    return resolved


@contextmanager
def expand_run_inputs(data_path: Path) -> Iterator[ExpandedRunInputs]:
    """Expand a CLI file, recursive folder, or ZIP source safely."""

    source = data_path.expanduser().resolve()
    if source.is_file() and source.suffix.lower() == ".zip":
        with extracted_zip_tree(source) as extracted_root:
            paths = _find_holo_inputs(extracted_root)
            if not paths:
                raise ValueError(f"No {HOLO_SUFFIX} files found in {source}")
            yield ExpandedRunInputs(tuple(paths), extracted_root)
        return

    paths = _find_holo_inputs(source)
    if not paths:
        raise ValueError(f"No {HOLO_SUFFIX} files found under {source}")
    batch_root = source if source.is_dir() else None
    yield ExpandedRunInputs(tuple(paths), batch_root)


def _find_holo_inputs(path: Path) -> list[Path]:
    if path.is_file():
        if path.suffix.lower() in {HOLO_SUFFIX, INPUT_LIST_SUFFIX}:
            return [path]
        raise ValueError(
            f"File is not a {HOLO_SUFFIX} or {INPUT_LIST_SUFFIX} file: {path}"
        )
    if path.is_dir():
        return sorted(
            candidate
            for candidate in path.rglob("*")
            if candidate.is_file() and candidate.suffix.lower() == HOLO_SUFFIX
        )
    raise FileNotFoundError(f"Input path does not exist: {path}")


def _output_manager_for_layout(
    layout: HoloRunLayout,
    *,
    output_root: Path | None,
    batch_root: Path,
) -> OutputManager:
    if output_root is None:
        return OutputManager(layout)
    relative_path = _relative_to_batch(layout.holo_path, batch_root)
    target_dir = output_root / relative_path.parent
    output_layout = HoloRunLayout.from_holo(
        layout.holo_path,
        output_root=target_dir,
    )
    return OutputManager(output_layout)


def _batch_root(holo_paths: Sequence[Path]) -> Path:
    if not holo_paths:
        return Path.cwd()
    if len(holo_paths) == 1:
        return holo_paths[0].parent
    try:
        return Path(os.path.commonpath([str(path.parent) for path in holo_paths]))
    except ValueError:
        return Path.cwd()


def _relative_to_batch(holo_path: Path, batch_root: Path) -> Path:
    try:
        return holo_path.relative_to(batch_root)
    except ValueError:
        anchor = Path(holo_path.anchor)
        drive_token = holo_path.drive.rstrip(":\\/") or "root"
        tail = holo_path.relative_to(anchor) if anchor != holo_path else Path()
        return Path(drive_token) / tail


def _reject_duplicate_destinations(requests: Sequence[RunRequest]) -> None:
    destinations: dict[str, Path] = {}
    for request in requests:
        destination = request.output_manager.layout.ef_dir.resolve(strict=False)
        key = os.path.normcase(str(destination))
        previous = destinations.get(key)
        if previous is not None:
            raise ValueError(
                "Multiple inputs resolve to the same EyeFlow output directory: "
                f"{previous} and {request.input_layout.holo_path} -> {destination}"
            )
        destinations[key] = request.input_layout.holo_path


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
