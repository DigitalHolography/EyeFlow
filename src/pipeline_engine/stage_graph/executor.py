"""Execute stages and coordinate handoffs through explicit responsibility owners."""

from collections.abc import Mapping
from contextlib import nullcontext
from time import perf_counter

from data_manager import DataManager
from input_output.schema.output_payload import DatasetValue
from input_output.writers.deferred import capture_artifact_writes
from input_output.writers.h5 import h5_write_checkpoint
from utils.logger import Logger

from pipeline_engine.execution.protocol import emit, InvocationOutcome
from pipeline_engine.resources import current_resources
from pipeline_engine.stage_graph.artifacts import ArtifactDispatcher
from pipeline_engine.stage_graph.journal import StageJournal
from pipeline_engine.stage_graph.products import ProductStore
from pipeline_engine.stage_graph.scope import GraphScope


class StageSkipped(Exception):
    """Expected missing scientific input, rather than a calculation failure."""


class StageAccess:
    """Public stage operations, without output transaction or scheduling internals."""

    def __init__(self, executor):
        self._executor = executor
        self.context = executor.context

    @property
    def active_node(self):
        return self._executor.scope.stage

    def require(self, name):
        execution = self._executor
        execution.scope.validate()
        execution.products.validate_read(execution.scope.stage, name)
        started = perf_counter()
        try:
            return execution.data.read(execution.products.require(name))
        finally:
            execution._add_stage_timing("product_read_seconds", perf_counter() - started)

    def get(self, name, default=None):
        execution = self._executor
        execution.scope.validate()
        execution.products.validate_read(execution.scope.stage, name)
        return self.require(name) if execution.products.has(name) else default

    def require_mutable(self, name):
        execution = self._executor
        execution.scope.validate()
        execution.products.validate_read(execution.scope.stage, name)
        started = perf_counter()
        try:
            return execution.data.writable(execution.products.require(name))
        finally:
            execution._add_stage_timing("mutable_product_copy_seconds", perf_counter() - started)

    def handoff(self, name, value):
        return self._executor.handoff(name, value)

    def skip(self, reason):
        self._executor.scope.validate()
        raise StageSkipped(reason)


class StageGraphExecutor:
    def __init__(self, definition, context, checkpoint, checkpoint_factory):
        self.definition, self.context = definition, context
        self.outer, self.checkpoint_factory = checkpoint, checkpoint_factory
        self.plan = definition.plan(external_products=context.product_inputs)
        self.products = ProductStore(self.plan)
        self.data = DataManager(context.runtime.work_h5.storage, owner_id=context.product_owner_id)
        self.scope = GraphScope(context)
        label = context.runtime.pipeline_name
        if context.execution_variant is not None and context.output.h5.processing_root:
            label += f" [{context.execution_variant.name}]"
        self.journal = StageJournal(self.plan, label)
        self.artifacts = ArtifactDispatcher(self.plan, self.products, self.journal, context)
        self.api = StageAccess(self)
        self._checkpoints, self._paths = [], {}
        self._writing = False
        self.pipeline_name = context.runtime.pipeline_name
        self.variant_name = (
            context.execution_variant.name if context.execution_variant is not None else "shared"
        )
        for alias, handle in context.product_inputs.items():
            if handle is not None:
                self.products.record_accepted(alias, handle)
            else:
                self.products.record_failed(alias, "Upstream product unavailable.")

    def _record(self, paths=()):
        started = perf_counter()
        try:
            self.journal.persist(self.context.runtime.work_h5)
            self._advance_checkpoints(paths)
            self.context.output.h5.flush()
        finally:
            self.data.storage.metrics.add("stage_journal_write_seconds", perf_counter() - started)

    def _add_stage_timing(self, name, duration):
        stage = self.scope.stage
        if stage is None:
            return
        timing = self.journal.timings.get(stage)
        if timing is not None:
            timing[name] = timing.get(name, 0.0) + duration

    def _advance_checkpoints(self, paths=()):
        for checkpoint in self._checkpoints:
            checkpoint.commit(paths)
        self.outer.commit(paths)

    def handoff(self, name, value):
        self.scope.validate()
        stage = self.plan.stage_by_id[self.scope.stage]
        self.products.validate_handoff(stage.id, name)
        spec = next(item for item in stage.products if item.id == name)
        if spec.validate:
            spec.validate(value)
        started = perf_counter()
        Logger.log_perf(
            f"START pipeline={self.pipeline_name} variant={self.variant_name} "
            f"stage={stage.id} operation=handoff product={name} at={started:.6f}"
        )
        handle = self.data.accept(value, transportable=spec.exposed)
        snapshot = self.data.read(handle)
        snapshotted = perf_counter()
        payload = spec.pack(snapshot, self.context) if spec.pack else (
            {spec.path: DatasetValue(snapshot, dict(spec.attrs))} if spec.path else {}
        )
        if not isinstance(payload, Mapping):
            raise TypeError("A product packer must return a schema-path/payload mapping.")
        packed = perf_counter()
        paths = [self.context.output.h5._path(path) for path in payload]
        seen = set(self._paths)
        for path in paths:
            if any(path == other or path.startswith(other + "/") or other.startswith(path + "/") for other in seen):
                raise ValueError("Handed-off HDF5 products must have distinct, non-overlapping dataset paths.")
            seen.add(path)
        self._writing = True
        try:
            with h5_write_checkpoint(self.context.runtime.work_h5):
                self.context.output.h5.write_many(payload)
        finally:
            self._writing = False
        self._advance_checkpoints(paths)
        self.products.record_accepted(name, handle)
        self._paths.update({path: name for path in paths})
        self.journal.committed(name, stage.id, paths, timing={
            "snapshot_seconds": snapshotted - started,
            "packing_seconds": packed - snapshotted,
            "commit_seconds": perf_counter() - packed,
            "array_bytes": self.data.array_bytes(handle),
            "dataset_bytes": sum(self.context.runtime.work_h5[path].size *
                                 self.context.runtime.work_h5[path].dtype.itemsize for path in paths),
        })
        handoff_timing = self.journal.products[name]
        Logger.log_perf(
            f"END pipeline={self.pipeline_name} variant={self.variant_name} "
            f"stage={stage.id} operation=handoff product={name} at={perf_counter():.6f} "
            f"wall={perf_counter() - started:.3f}s "
            f"snapshot={handoff_timing['snapshot_seconds']:.3f}s "
            f"packing={handoff_timing['packing_seconds']:.3f}s "
            f"output_commit={handoff_timing['commit_seconds']:.3f}s"
        )
        self._record()
        self.data.storage.metrics.add("local_acceptance_seconds", perf_counter() - started)
        # Transport and coordinator acceptance are collaborators, never store policy.
        self.context._product_committed(stage.id, name, handle if spec.exposed else None)
        Logger.log(f"[HANDOFF] {self.context.runtime.pipeline_name}/{stage.id}: {name} committed.")
        self.context._start_artifacts()
        self.artifacts.notify()
        self.artifacts.dispatch(self._run_stage)
        return snapshot

    def _guard_write(self):
        if not self._writing:
            raise RuntimeError("Graph stages persist products through ctx.dag.handoff(), not direct HDF5 writes.")

    def run(self):
        ctx = self.context
        previous = ctx.output.h5._write_guard
        ctx.output.h5._write_guard = self._guard_write
        try:
            with self.scope.activate(self.api):
                self._record()
                emit("GraphStarted", journal=self.journal.attributes())
                for name in self.plan.order:
                    stage = self.plan.stage_by_id[name]
                    if self.journal.statuses[name] != "pending":
                        continue
                    missing = [key for key in stage.requires if not self.products.has(key)]
                    missing += [key for key in stage.after if self.journal.statuses[key] != "succeeded"]
                    if missing:
                        self.journal.mark(name, "skipped", "Unavailable dependencies: " + ", ".join(missing))
                        self._record()
                        emit("StageSkipped", stage=name, error=self.journal.skips[name])
                        Logger.log_perf(
                            f"pipeline={self.pipeline_name} variant={self.variant_name} "
                            f"stage={name} status=skipped at={self.journal.timings[name]['end']:.6f} "
                            f"reason=unavailable_dependencies"
                        )
                        self.artifacts.dispatch(self._run_stage)
                    else:
                        self._run_stage(stage)
                if self.journal.failures:
                    raise RuntimeError("; ".join(f"{name}: {error}" for name, error in self.journal.failures.items()))
                blocked = [stage.id for stage in self.definition.stages
                           if self.journal.statuses[stage.id] == "skipped" and not stage.artifact
                           and (stage.requires or stage.after or any(not p.optional for p in stage.products))]
                if blocked or (self.journal.statuses and set(self.journal.statuses.values()) == {"skipped"}):
                    supplied = any(self.products.has(p.id) for stage in self.definition.stages for p in stage.products)
                    reason = "; ".join(f"{name}: {value}" for name, value in self.journal.skips.items())
                    return InvocationOutcome(ctx.runtime.pipeline_name, status="partial" if supplied else "skipped", reason=reason)
                return InvocationOutcome(ctx.runtime.pipeline_name)
        finally:
            ctx.output.h5._write_guard = previous
            try:
                if self.definition.cleanup is not None:
                    self.definition.cleanup()
            finally:
                self.data.close_resources()

    def _run_stage(self, stage):
        ctx = self.context
        wait_started = perf_counter()
        if not self.artifacts.wait(stage):
            self.journal.mark(stage.id, "skipped", "Required artifact export failed.")
            self._record()
            emit("StageSkipped", stage=stage.id, error=self.journal.skips[stage.id])
            Logger.log_perf(
                f"pipeline={self.pipeline_name} variant={self.variant_name} "
                f"stage={stage.id} status=skipped artifact_wait={perf_counter() - wait_started:.3f}s"
            )
            return
        artifact_waited = perf_counter() - wait_started
        previous_figure = ctx._figure_job, ctx._figure_count
        ctx._figure_job, ctx._figure_count = None, 0
        jobs_before = len(ctx._artifact_jobs)
        self.journal.mark(stage.id, "running")
        timing = self.journal.timings[stage.id]
        ready_at = max([self.journal.created_at,
            *(self.products.resolved_at[name] for name in (*stage.requires, *stage.optional_requires)),
            *(self.journal.timings[name]["end"] for name in stage.after)])
        timing.update(ready_at=ready_at, dependency_wait_seconds=ready_at - self.journal.created_at,
                      dispatch_wait_seconds=max(0.0, timing["start"] - ready_at - artifact_waited),
                      artifact_wait_seconds=artifact_waited)
        Logger.log_perf(
            f"START pipeline={self.pipeline_name} variant={self.variant_name} "
            f"stage={stage.id} at={timing['start']:.6f}"
        )
        self._record()
        emit("StageStarted", stage=stage.id, journal=self.journal.attributes())
        try:
            with self.checkpoint_factory(ctx.runtime.work_h5, ctx.state.raw, ctx.output.manager, ctx._background_artifacts) as checkpoint:
                self._checkpoints.append(checkpoint)
                try:
                    capture = None if stage.wait_for_artifacts else ctx._capture_artifact_write
                    with self.scope.executing(stage.id), capture_artifact_writes(capture):
                        budget = current_resources()
                        acquired = perf_counter()
                        with budget.compute() if budget else nullcontext():
                            callback_started = perf_counter()
                            timing["resource_wait_seconds"] = callback_started - acquired
                            try:
                                result = stage.run(ctx)
                            finally:
                                timing["callback_wall_seconds"] = perf_counter() - callback_started
                    if result is not None:
                        raise TypeError("Graph stages hand off products and return None.")
                    missing = [p.id for p in stage.products if not p.optional and not self.products.has(p.id)]
                    if missing:
                        raise ValueError("Missing handoffs: " + ", ".join(missing))
                    checkpoint.commit()
                finally:
                    self._checkpoints.pop()
        except StageSkipped as exc:
            ctx._cancel_artifacts()
            self.journal.mark(stage.id, "skipped", exc)
        except Exception as exc:
            ctx._cancel_artifacts()
            self.journal.mark(stage.id, "failed", exc)
            Logger.log_error(f"{ctx.runtime.pipeline_name}/{stage.id}: {exc}")
        else:
            self.journal.mark(stage.id, "succeeded")
            ctx._start_artifacts()
        finally:
            self.artifacts.jobs[stage.id] = ctx._artifact_jobs[jobs_before:]
            ctx._figure_job, ctx._figure_count = previous_figure
        for product in stage.products:
            if not self.products.has(product.id):
                self.products.record_failed(product.id, "Producer finished without supplying product.")
        self._record()
        emit({"failed": "StageFailed", "skipped": "StageSkipped", "succeeded": "StageFinished"}[self.journal.statuses[stage.id]],
             stage=stage.id, error=self.journal.failures.get(stage.id) or self.journal.skips.get(stage.id))
        self.artifacts.dispatch(self._run_stage)
        Logger.log_perf(
            f"END pipeline={self.pipeline_name} variant={self.variant_name} "
            f"stage={stage.id} status={self.journal.statuses[stage.id]} "
            f"at={self.journal.timings[stage.id].get('end', perf_counter()):.6f} "
            f"wall={self.journal.timings[stage.id].get('wall_seconds', 0.0):.3f}s "
            f"callback={self.journal.timings[stage.id].get('callback_wall_seconds', 0.0):.3f}s "
            f"resource_wait={self.journal.timings[stage.id].get('resource_wait_seconds', 0.0):.3f}s "
            f"artifact_wait={self.journal.timings[stage.id].get('artifact_wait_seconds', 0.0):.3f}s "
            f"dispatch_wait={self.journal.timings[stage.id].get('dispatch_wait_seconds', 0.0):.3f}s "
            f"product_read={self.journal.timings[stage.id].get('product_read_seconds', 0.0):.3f}s "
            f"mutable_copy={self.journal.timings[stage.id].get('mutable_product_copy_seconds', 0.0):.3f}s"
        )
