"""Small queue events; large invocation results travel through local files."""

from contextvars import ContextVar
from dataclasses import dataclass, field

_SINK = ContextVar("eyeflow_worker_events", default=None)


@dataclass(frozen=True)
class WorkerEvent:
    kind: str
    invocation: str
    stage: str | None = None
    detail: dict = field(default_factory=dict)


def emit(kind, *, stage=None, **detail):
    sink = _SINK.get()
    if sink is not None:
        sink(kind, stage, detail)
