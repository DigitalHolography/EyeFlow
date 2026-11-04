"""Choose dependency-ready work; process transport belongs to ProcessExecutor."""

from utils.logger import Logger
from .product_exchange import ProductStatus


class CoordinatorScheduler:
    def __init__(
        self, *, pipelines, dag, options, executor, statuses, failures, skips,
        on_start=None, on_success=None, on_progress=None, registry=None,
    ):
        self.pending = list(pipelines)
        self.dag, self.options, self.executor = dag, options, executor
        self.statuses, self.failures, self.skips = statuses, failures, skips
        self.on_start, self.on_success, self.on_progress = on_start, on_success, on_progress
        self.completed = set()
        self.total = len(pipelines)
        self.started = 0
        self.registry = registry

    def _dependencies(self, pipeline):
        return self.dag.dependencies_of(pipeline.name, pipeline_options=self.options, include_products=False)

    def _products_ready(self, pipeline):
        if not pipeline.product_requires:
            return True
        identities = self.registry.invocations_for(pipeline.name)
        return bool(identities) and all(
            self.registry.status(requirement.resolve(identity)) != ProductStatus.PENDING
            for identity in identities for requirement in pipeline.product_requires)

    def _finalize(self, pipeline, error=None):
        for identity in self.registry.invocations_for(pipeline.name):
            self.registry.finalize_producer(identity, error or "Producer did not supply product.")

    def _notify_start(self, pipeline):
        self.started += 1
        if self.on_start:
            self.on_start(pipeline.name, self.started, self.total)

    def _notify_completion(self):
        if self.on_progress:
            self.on_progress()

    def _collect(self, *, wait=False):
        for outcome in self.executor.collect(wait=wait):
            self.completed.add(outcome.name)
            self.statuses[outcome.name] = outcome.status
            if outcome.error:
                self.failures[outcome.name] = outcome.error
                Logger.log_warning(f"Failed {outcome.name}; keeping completed products: {outcome.error}")
            elif outcome.reason:
                self.skips[outcome.name] = outcome.reason
            elif outcome.status == "succeeded" and self.on_success:
                self.on_success(outcome.name)
            self._notify_completion()

    def _start_processes(self):
        for pipeline in tuple(self.pending):
            if pipeline.execution_mode != "process":
                continue
            dependencies = self._dependencies(pipeline)
            if not set(dependencies) <= self.completed:
                continue
            if not self._products_ready(pipeline):
                continue
            blocked = [name for name in dependencies if self.statuses.get(name) != "succeeded"]
            if blocked:
                self.pending.remove(pipeline)
                self.statuses[pipeline.name] = "skipped"
                self.skips[pipeline.name] = "Unavailable dependencies: " + ", ".join(blocked)
                self.completed.add(pipeline.name)
                self._finalize(pipeline, self.skips[pipeline.name])
                self._notify_start(pipeline)
                self._notify_completion()
                continue
            try:
                admitted = self.executor.start(pipeline)
            except Exception as exc:
                self.statuses[pipeline.name] = "failed"
                self.failures[pipeline.name] = str(exc)
                self.completed.add(pipeline.name)
                self._finalize(pipeline, str(exc))
                self.pending.remove(pipeline)
                self._notify_start(pipeline)
                self._notify_completion()
                continue
            if admitted:
                self.pending.remove(pipeline)
                self.statuses[pipeline.name] = "running"
                self._notify_start(pipeline)

    def _next_coordinator(self):
        if not self.executor.coordinator_ready:
            return None
        numerical_pending = self.executor.running or any(
            not pipeline.wait_for_artifacts for pipeline in self.pending
        )
        for pipeline in self.pending:
            if pipeline.execution_mode != "coordinator":
                continue
            if pipeline.wait_for_artifacts and numerical_pending:
                continue
            if set(self._dependencies(pipeline)) <= self.completed and self._products_ready(pipeline):
                return pipeline
        return None

    def coordinator_jobs(self):
        """Yield coordinator descriptors; their caller records invocation outcomes."""
        while self.pending or self.executor.running:
            self._collect()
            self._start_processes()
            selected = self._next_coordinator()
            if selected is not None:
                self.pending.remove(selected)
                # Coordinator callbacks remain with the invocation's caller.
                self.started += 1
                yield self.started, selected
                self.completed.add(selected.name)
                self._finalize(selected, self.failures.get(selected.name) or self.skips.get(selected.name))
            elif self.executor.running:
                self._collect(wait=True)
            elif self.pending:
                if any(
                    p.execution_mode == "process"
                    and set(self._dependencies(p)) <= self.completed
                    for p in self.pending
                ):
                    # Export/finalization reservations can temporarily fill the
                    # budget even when no numerical worker is active.
                    self._collect(wait=True)
                else:
                    raise RuntimeError("No schedulable pipeline; inspect active dependencies.")
