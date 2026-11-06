"""Stage diagnostics and their single reconciliation policy across invocations."""

import json
from enum import Enum
from contextlib import nullcontext
from os import getpid
from time import perf_counter


class StageState(str, Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    SKIPPED = "skipped"


class StageJournal:
    keys = ("stage_status", "stage_failures", "stage_skips", "stage_products", "stage_timings")

    def __init__(self, plan, scope):
        self.scope = scope
        self.statuses = {name: "pending" for name in plan.order}
        self.failures, self.skips, self.products = {}, {}, {}
        self.timings = {}
        self.created_at = perf_counter()

    def mark(self, stage, status, reason=None):
        status = StageState(status).value
        self.statuses[stage] = status
        now = perf_counter()
        timing = self.timings.setdefault(stage, {"process_id": getpid()})
        if status == "running":
            timing["start"] = now
        else:
            timing["end"] = now
            if "start" in timing:
                timing["wall_seconds"] = now - timing["start"]
        if reason is not None and status in {"failed", "skipped"}:
            (self.failures if status == "failed" else self.skips)[stage] = str(reason)

    def committed(self, product, stage, paths, *, timing=None):
        self.products[product] = {"producer": stage, "datasets": list(paths), **(timing or {})}

    def persist(self, model):
        self.apply(model, self.reconcile(model, self.attributes()))

    @classmethod
    def apply(cls, model, attributes):
        with model.writing_as("stage_journal") if hasattr(model, "writing_as") else nullcontext():
            for key, value in attributes.items():
                if key not in cls.keys:
                    raise ValueError("Only stage journal root attributes belong to StageJournal.")
                model.attrs[key] = value

    def attributes(self):
        return {key: json.dumps({self.scope: values}, sort_keys=True)
                for key, values in zip(self.keys, (self.statuses, self.failures, self.skips, self.products, self.timings))}

    @classmethod
    def reconcile(cls, model, attributes, *, allowed_scope=None):
        """Merge by qualified invocation *and* stage/product, never replace peers."""
        result = {key: value for key, value in attributes.items() if key not in cls.keys}
        for key in cls.keys:
            incoming = json.loads(attributes.get(key, "{}"))
            if allowed_scope is not None and set(incoming) - {allowed_scope}:
                raise ValueError("Worker diagnostics belong to another invocation scope.")
            if not incoming:
                continue
            existing = json.loads(model.attrs.get(key, "{}"))
            for scope, records in incoming.items():
                existing.setdefault(scope, {}).update(records)
            result[key] = json.dumps(existing, sort_keys=True)
        return result

    @classmethod
    def observe(cls, model, scope, stage, status, error=None):
        attrs = {"stage_status": json.dumps({scope: {stage: status}})}
        if error and status in {"failed", "skipped"}:
            key = "stage_failures" if status == "failed" else "stage_skips"
            attrs[key] = json.dumps({scope: {stage: error}})
        cls.apply(model, cls.reconcile(model, attrs))

    @classmethod
    def finalize(cls, model, scope, error):
        if not error:
            return
        current = json.loads(model.attrs.get("stage_status", "{}"))
        for stage, status in tuple(current.get(scope, {}).items()):
            if status in {"pending", "running"}:
                cls.observe(model, scope, stage, "failed" if status == "running" else "skipped", error)
