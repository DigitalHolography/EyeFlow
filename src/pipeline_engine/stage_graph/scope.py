"""Active graph/stage lifetime, context restoration and owning-thread checks."""

from contextlib import contextmanager
from contextvars import ContextVar
from threading import get_ident

_CURRENT = ContextVar("eyeflow_stage_graph_scope", default=None)


class GraphScope:
    def __init__(self, context):
        self.context = context
        self.thread = get_ident()
        self.stages = []
        self.api = None

    @property
    def stage(self):
        return self.stages[-1] if self.stages else None

    def validate(self):
        if get_ident() != self.thread:
            raise RuntimeError("DAG access must run on the owning pipeline thread.")
        if _CURRENT.get() is not self or self.context._dag_execution is not self.api:
            raise RuntimeError("This DAG belongs to another or completed pipeline invocation.")
        if self.stage is None:
            raise RuntimeError("A product operation requires an executing stage.")

    @contextmanager
    def activate(self, api):
        self.api = api
        previous = self.context._dag_execution
        token = _CURRENT.set(self)
        self.context._dag_execution = api
        try:
            yield
        finally:
            self.context._dag_execution = previous
            _CURRENT.reset(token)

    @contextmanager
    def executing(self, stage):
        self.stages.append(stage)
        try:
            yield
        finally:
            self.stages.pop()


def current_pipeline_dag():
    scope = _CURRENT.get()
    if scope is None:
        raise RuntimeError("There is no pipeline DAG in this execution scope.")
    scope.validate()
    return scope.api


def handoff(name, value):
    current_pipeline_dag().handoff(name, value)
