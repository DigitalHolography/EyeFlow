"""Pipeline-owned stage definitions and invocation-local product handoff."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from contextvars import ContextVar
from copy import deepcopy
from dataclasses import dataclass
import json
from threading import get_ident

from input_output.writers.h5 import h5_write_checkpoint
from input_output.writers.deferred import capture_artifact_writes
from utils.logger import Logger
from .protocol import emit
from .resources import current_resources, suspend_compute
from contextlib import nullcontext

from .base import DatasetValue, PipelineDescriptor
from .dag import PipelineDAG

_CURRENT = ContextVar("eyeflow_pipeline_dag", default=None)


def current_pipeline_dag():
    graph = _CURRENT.get()
    if graph is None:
        raise RuntimeError("There is no pipeline DAG in this execution scope.")
    graph._check_thread()
    return graph


def handoff(name, value):
    """Convenience for deeply nested helpers; resolves the current local DAG."""
    current_pipeline_dag().handoff(name, value)


@dataclass(frozen=True)
class PipelineProduct:
    name: str
    path: str | None = None
    attrs: Mapping | None = None
    pack: Callable | None = None  # (value, ctx) -> schema-path/payload mapping
    validate: Callable | None = None
    optional: bool = False


@dataclass(frozen=True)
class PipelineNode:
    name: str
    run: Callable
    requires: tuple[str, ...] = ()
    outputs: tuple[PipelineProduct, ...] = ()
    optional_requires: tuple[str, ...] = ()
    after: tuple[str, ...] = ()
    artifact: bool = False
    wait_for_artifacts: bool = False


@dataclass(frozen=True)
class PipelineGraph:
    nodes: tuple[PipelineNode, ...]
    cleanup: Callable | None = None

    def plan(self):
        """Use the existing planner; local names never enter the global catalog."""
        producers = {}
        nodes = {node.name: node for node in self.nodes}
        if len(nodes) != len(self.nodes):
            raise ValueError("Duplicate pipeline stage name.")
        for node in self.nodes:
            if not node.name or not callable(node.run):
                raise ValueError("Each stage needs a name and a callable.")
            for output in node.outputs:
                if not output.name or output.name in producers:
                    raise ValueError(f"Duplicate or empty stage product: {output.name!r}")
                if output.path is not None and output.pack is not None:
                    raise ValueError(f"Product {output.name} declares both path and pack.")
                producers[output.name] = node.name
        for node in self.nodes:
            unknown = set((*node.requires, *node.optional_requires)) - producers.keys()
            unknown_stages = set(node.after) - nodes.keys()
            if unknown or unknown_stages:
                missing = sorted(unknown) + [f"stage:{name}" for name in sorted(unknown_stages)]
                raise ValueError(f"Stage {node.name} has undeclared dependencies: {missing}")
            if node.name in node.after:
                raise ValueError(f"Stage {node.name} depends on itself.")
            if any(producers[key] == node.name for key in (*node.requires, *node.optional_requires)):
                raise ValueError(f"Stage {node.name} depends on its own output.")
        planner = PipelineDAG([
            PipelineDescriptor(
                name=f"node/{node.name}", description=node.name, available=True,
                dag_requires=tuple(f"product/{key}" for key in (*node.requires, *node.optional_requires))
                + tuple(f"node/{name}" for name in node.after),
                dag_produces=tuple(f"product/{output.name}" for output in node.outputs),
            ) for node in self.nodes
        ])
        return planner, nodes, producers


class StageSkipped(Exception):
    """An expected absence of source data, rather than a scientific failure."""


class GraphExecution:
    def __init__(self, definition, context, checkpoint, checkpoint_factory):
        self.definition = definition
        self.context = context
        self._outer = checkpoint
        self._checkpoint_factory = checkpoint_factory
        self._planner, self._nodes, self._producers = definition.plan()
        self._values = {}
        self._paths = {}
        self.products = {}
        self._active = []
        self._thread = get_ident()
        self.statuses = {name: "pending" for name in self._nodes}
        self.failures = {}
        self.skips = {}
        self._jobs = {}
        self._writing = False
        self._dispatching = False

    @property
    def active_node(self):
        return self._active[-1][0] if self._active else None

    def _check_thread(self):
        if get_ident() != self._thread:
            raise RuntimeError("DAG handoffs must run on the owning pipeline thread; return worker results to it.")
        if _CURRENT.get() is not self:
            raise RuntimeError("This DAG belongs to another or completed pipeline invocation.")

    def require(self, name):
        self._check_thread()
        node = self._nodes[self.active_node]
        if name not in (*node.requires, *node.optional_requires):
            raise ValueError(f"Stage {node.name} did not declare input {name!r}.")
        if name not in self._values:
            raise KeyError(f"Product {name!r} has not been handed off.")
        return deepcopy(self._values[name])

    def get(self, name, default=None):
        self._check_thread()
        if name not in self._values:
            # Validate the declared read even when its optional producer failed.
            node = self._nodes[self.active_node]
            if name not in (*node.requires, *node.optional_requires):
                raise ValueError(f"Stage {node.name} did not declare input {name!r}.")
            return default
        return self.require(name)

    def skip(self, reason):
        self._check_thread()
        raise StageSkipped(reason)

    def handoff(self, name, value):
        self._check_thread()
        if not self._active:
            raise RuntimeError("A handoff requires an executing stage.")
        node = self._nodes[self.active_node]
        declared = {output.name: output for output in node.outputs}
        if name not in declared or name in self._values:
            raise ValueError(f"Undeclared or duplicate product handoff: {name!r}")
        spec = declared[name]
        snapshot = deepcopy(value)
        if spec.validate is not None:
            spec.validate(snapshot)
        payload = spec.pack(snapshot, self.context) if spec.pack else (
            {spec.path: DatasetValue(snapshot, dict(spec.attrs or {}))} if spec.path else {}
        )
        if not isinstance(payload, Mapping):
            raise TypeError("A product packer must return a schema-path/payload mapping.")
        paths = [self.context.output.h5._path(path) for path in payload]
        seen = set(self._paths)
        for path in paths:
            if any(path == other or path.startswith(other + "/") or other.startswith(path + "/") for other in seen):
                raise ValueError("Handed-off HDF5 products must have distinct, non-overlapping dataset paths.")
            seen.add(path)
        # Packing/writing one product is atomic, even if its caller catches an error.
        self._writing = True
        try:
            if payload:
                with h5_write_checkpoint(self.context.runtime.work_h5):
                    self.context.output.h5.write_many(payload)
        finally:
            self._writing = False
        self._values[name] = snapshot
        self._paths.update({path: name for path in paths})
        self.products[name] = {"producer": node.name, "datasets": paths}
        self._record(paths)
        emit("ProductCommit", stage=node.name, product=name)
        Logger.log(f"[HANDOFF] {self.context.runtime.pipeline_name}/{node.name}: {name} committed.")
        self.context._start_artifacts()
        self._dispatch_artifacts()

    def _guard_write(self):
        if not self._writing:
            raise RuntimeError("Graph stages persist products through ctx.dag.handoff(), not direct HDF5 writes.")

    def run(self):
        ctx = self.context
        token = _CURRENT.set(self)
        ctx._dag_execution = self
        old_guard = getattr(ctx.output.h5, "_write_guard", None)
        ctx.output.h5._write_guard = self._guard_write
        try:
            for name in self._planner.execution_order:
                node = self._nodes[name.removeprefix("node/")]
                if self.statuses[node.name] != "pending":
                    continue
                missing = [key for key in node.requires if key not in self._values]
                missing.extend(name for name in node.after if self.statuses[name] != "succeeded")
                if missing:
                    self.statuses[node.name] = "skipped"
                    self.skips[node.name] = "Unavailable dependencies: " + ", ".join(missing)
                    self._record()
                    emit("StageSkipped", stage=node.name, error=self.skips[node.name])
                    continue
                self._run_node(node)
            if self.failures:
                raise RuntimeError("; ".join(f"{name}: {error}" for name, error in self.failures.items()))
        finally:
            ctx.output.h5._write_guard = old_guard
            ctx._dag_execution = None
            _CURRENT.reset(token)
            if self.definition.cleanup is not None:
                self.definition.cleanup()

    def _dispatch_artifacts(self):
        if self._dispatching:
            return
        self._dispatching = True
        try:
            for name in self._planner.execution_order:
                node = self._nodes[name.removeprefix("node/")]
                if not node.artifact or self.statuses[node.name] != "pending":
                    continue
                if any(key not in self._values for key in node.requires):
                    continue
                if any(key not in self._values and self.statuses[self._producers[key]] in {"pending", "running"}
                       for key in node.optional_requires):
                    continue
                if any(self.statuses[name] != "succeeded" for name in node.after):
                    continue
                self._run_node(node)
        finally:
            self._dispatching = False

    def _run_node(self, node):
        ctx = self.context
        if node.wait_for_artifacts and ctx._background_artifacts is not None:
            parents = self._planner.dependencies_of(f"node/{node.name}", transitive=True)
            jobs = [job for name in parents for job in self._jobs.get(name.removeprefix("node/"), ())
                    if job.future is not None]
            with suspend_compute():
                for job in jobs:
                    job.future.result()
            if any(job.status == "failed" for job in jobs):
                self.statuses[node.name] = "skipped"
                self.skips[node.name] = "Required artifact export failed."
                self._record()
                return
        previous_figure = (ctx._figure_job, ctx._figure_count)
        ctx._figure_job, ctx._figure_count = None, 0
        jobs_before = len(ctx._artifact_jobs)
        self.statuses[node.name] = "running"
        emit("StageStarted", stage=node.name)
        try:
            with self._checkpoint_factory(
                ctx.runtime.work_h5, ctx.state.raw, ctx.output.manager, ctx._background_artifacts,
            ) as checkpoint:
                self._active.append((node.name, checkpoint))
                try:
                    capture = None if node.wait_for_artifacts else ctx._capture_artifact_write
                    with capture_artifact_writes(capture):
                        budget = current_resources()
                        allowance = min(4, max(1, budget.cpu_limit - budget.cpu_used)) if budget else 1
                        memory = min(64 * 2**20, budget.memory_limit // 4) if budget else 0
                        with budget.compute(allowance, memory) if budget else nullcontext():
                            result = node.run(ctx)
                    if result is not None:
                        raise TypeError("Graph stages hand off products and return None.")
                    missing = [spec.name for spec in node.outputs if not spec.optional and spec.name not in self._values]
                    if missing:
                        raise ValueError("Missing handoffs: " + ", ".join(missing))
                    checkpoint.commit()
                finally:
                    self._active.pop()
        except StageSkipped as exc:
            ctx._cancel_artifacts()
            self.statuses[node.name] = "skipped"
            self.skips[node.name] = str(exc)
        except Exception as exc:
            ctx._cancel_artifacts()
            self.statuses[node.name] = "failed"
            self.failures[node.name] = str(exc)
            Logger.log_error(f"{ctx.runtime.pipeline_name}/{node.name}: {exc}")
        else:
            self.statuses[node.name] = "succeeded"
            ctx._start_artifacts()
        finally:
            self._jobs[node.name] = ctx._artifact_jobs[jobs_before:]
            ctx._figure_job, ctx._figure_count = previous_figure
        self._record()
        emit({"failed": "StageFailed", "skipped": "StageSkipped",
              "succeeded": "StageFinished"}[self.statuses[node.name]], stage=node.name,
             error=self.failures.get(node.name) or self.skips.get(node.name))
        self._dispatch_artifacts()

    def _record(self, paths=()):
        ctx = self.context
        scope = ctx.runtime.pipeline_name
        if ctx.execution_variant is not None and ctx.output.h5.processing_root:
            scope += f" [{ctx.execution_variant.name}]"
        for attribute, value in (("stage_status", self.statuses), ("stage_failures", self.failures),
                                 ("stage_skips", self.skips), ("stage_products", self.products)):
            records = json.loads(ctx.runtime.work_h5.attrs.get(attribute, "{}"))
            records[scope] = value
            ctx.runtime.work_h5.attrs[attribute] = json.dumps(records, sort_keys=True)
        for _, checkpoint in self._active:
            checkpoint.commit(paths)
        self._outer.commit(paths)
        ctx.output.h5.flush()
