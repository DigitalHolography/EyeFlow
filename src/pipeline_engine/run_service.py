"""Shared GUI/CLI orchestration for EyeFlow pipeline runs."""

from __future__ import annotations

import os
import shutil
from collections.abc import Callable, Iterable, Iterator, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from uuid import uuid4

from input_output import INPUT_LIST_SUFFIX, HoloRunLayout, resolve_selected_run_layouts
from input_output.archives import extracted_zip_tree
from input_output.output_manager import OutputManager, OutputType

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
    return RunSpec(plan=plan, requests=requests)


def execute_run(
    spec: RunSpec,
    *,
    on_log: Callable[[str], None] | None = None,
    on_progress: Callable[[], None] | None = None,
) -> RunResult:
    """Execute a batch and transactionally commit each successful output."""

    _emit(on_log, f"[DAG] Targets -> {', '.join(spec.plan.targets)}")
    _emit(on_log, f"[DAG] Execution order -> {', '.join(spec.plan.names)}")
    outputs: list[Path] = []
    failures: list[RunFailure] = []

    for request in spec.requests:
        input_layout = request.input_layout
        final_manager = request.output_manager
        final_path = final_manager.path_for(OutputType.H5)
        _emit(on_log, f"[INPUT] HOLO -> {input_layout.holo_path}")
        _emit(on_log, f"[INPUT] DATA DIR -> {input_layout.root_dir}")
        _emit(on_log, f"[RESOLVED] HD -> {input_layout.hd_h5}")
        _emit(on_log, f"[RESOLVED] DV -> {input_layout.dv_h5}")
        _emit(on_log, f"[OUTPUT] {final_path}")

        stage_manager = _staging_output_manager(final_manager)

        def forward_runtime_log(message: str) -> None:
            # Pipelines write into staging, but staging is an implementation
            # detail and its paths become invalid immediately after commit.
            visible_message = message.replace(
                str(stage_manager.layout.root_dir),
                str(final_manager.layout.root_dir),
            )
            _emit(on_log, visible_message)

        try:
            run_pipelines_to_output(
                output_manager=stage_manager,
                pipelines=spec.plan.descriptors,
                target_names=spec.plan.targets,
                holodoppler_h5=input_layout.hd_h5,
                doppler_vision_h5=input_layout.dv_h5,
                on_log=forward_runtime_log,
                on_progress=on_progress,
            )
            _commit_staged_output(stage_manager, final_manager)
        except Exception as exc:  # noqa: BLE001
            _discard_staged_output(stage_manager)
            failure = RunFailure(input_layout.holo_path, str(exc))
            failures.append(failure)
            _emit(on_log, f"[FAIL] {failure}")
            continue

        outputs.append(final_path)
        _emit(on_log, f"Completed run for {input_layout.holo_path.name}: {final_path}")

    return RunResult(tuple(outputs), tuple(failures))


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


def _staging_output_manager(final_manager: OutputManager) -> OutputManager:
    layout = final_manager.layout
    stage_root = layout.root_dir.parent / (
        f".{layout.root_dir.name}.eyeflow-staging-{uuid4().hex}"
    )
    return OutputManager(layout.with_root_dir(stage_root))


def _commit_staged_output(
    stage_manager: OutputManager,
    final_manager: OutputManager,
) -> None:
    stage_dir = stage_manager.layout.ef_dir
    final_dir = final_manager.layout.ef_dir
    if not stage_dir.is_dir():
        raise RuntimeError(f"Pipeline run did not create its staged output: {stage_dir}")
    if final_dir.exists() and not final_dir.is_dir():
        raise RuntimeError(
            f"Refusing to replace non-directory EyeFlow output path: {final_dir}"
        )

    final_dir.parent.mkdir(parents=True, exist_ok=True)
    if final_dir.exists():
        shutil.rmtree(final_dir)
    try:
        stage_dir.replace(final_dir)
    finally:
        _remove_empty_stage_root(stage_manager.layout.root_dir)


def _discard_staged_output(stage_manager: OutputManager) -> None:
    stage_root = stage_manager.layout.root_dir
    if stage_root.exists():
        shutil.rmtree(stage_root, ignore_errors=True)


def _remove_empty_stage_root(stage_root: Path) -> None:
    try:
        stage_root.rmdir()
    except OSError:
        pass


def _emit(on_log: Callable[[str], None] | None, message: str) -> None:
    if on_log is not None:
        on_log(message)


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
