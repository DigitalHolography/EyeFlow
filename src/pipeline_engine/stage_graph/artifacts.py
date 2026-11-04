"""Early artifact readiness, nonrecursive dispatch and completion barriers."""

from ..resources import suspend_compute
from time import perf_counter


class ArtifactDispatcher:
    def __init__(self, plan, products, journal, context):
        self.plan, self.products, self.journal, self.context = plan, products, journal, context
        self.jobs = {}
        self._dispatching = False
        self._notified_at = None

    def notify(self):
        self._notified_at = perf_counter()

    def ready(self, stage):
        if any(not self.products.has(key) for key in stage.requires):
            return False
        if any(not self.products.has(key)
               and self.journal.statuses.get(self.plan.producer_by_product[key]) in {"pending", "running"}
               for key in stage.optional_requires):
            return False
        return all(self.journal.statuses[name] == "succeeded" for name in stage.after)

    def dispatch(self, run_stage):
        if self._dispatching:
            return
        self._dispatching = True
        try:
            for name in self.plan.order:
                stage = self.plan.stage_by_id[name]
                if stage.artifact and self.journal.statuses[name] == "pending" and self.ready(stage):
                    if self._notified_at is not None:
                        metrics = self.context.runtime.work_h5.storage.metrics
                        metrics.add("artifact_dispatch_seconds", perf_counter() - self._notified_at)
                        metrics.add("artifact_dispatches", 1)
                    run_stage(stage)
        finally:
            self._dispatching = False

    def wait(self, stage):
        if not stage.wait_for_artifacts:
            return True
        jobs = [job for parent in self.plan.ancestors(stage.id)
                for job in self.jobs.get(parent, ()) if job.future is not None]
        with suspend_compute():
            try:
                for job in jobs:
                    job.future.result()
            except Exception:
                return False
        return not any(job.status == "failed" for job in jobs)
