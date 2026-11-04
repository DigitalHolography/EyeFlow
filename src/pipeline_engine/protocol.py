"""Small queue events; large invocation results travel through local files."""

from contextvars import ContextVar
from dataclasses import dataclass
from pathlib import Path
from typing import Literal
from .product_exchange.identity import PipelineInvocationId

_SINK = ContextVar("eyeflow_worker_events", default=None)


@dataclass(frozen=True)
class WorkerEvent:
    kind: Literal[
        "WorkerStarted", "GraphStarted", "StageStarted", "ProductCommit", "StageFinished",
        "StageFailed", "StageSkipped", "InvocationFinished", "LogMessage",
    ]
    invocation: str
    stage: str | None = None
    result_path: str | None = None
    error: str | None = None
    product: str | None = None
    message: str | None = None
    commit_id: str | None = None
    invocation_id: PipelineInvocationId | None = None
    handle: object = None
    emitted_at: float | None = None
    journal: dict | None = None


@dataclass(frozen=True)
class ProductAcceptedEvent(WorkerEvent):
    """A locally committed product; global availability still needs acceptance."""


@dataclass(frozen=True)
class TransportReceipt:
    """Receipt is deliberately not an authoritative acceptance acknowledgment."""

    commit_id: str | None
    phase: Literal["received"] = "received"


@dataclass(frozen=True)
class InvocationSettings:
    """Resolved context settings shared by coordinator and spawned invocations."""

    configs: dict
    pipeline_options: dict
    pipeline_order: tuple[str, ...]
    pipeline_targets: tuple[str, ...]
    velocity_workflow_selection: str
    band_ratio_frequency_scale_hz: float

    def context_arguments(self):
        return {
            "pipeline_options": self.pipeline_options,
            "pipeline_order": self.pipeline_order,
            "pipeline_targets": self.pipeline_targets,
            "velocity_workflow_selection": self.velocity_workflow_selection,
            "band_ratio_frequency_scale_hz": self.band_ratio_frequency_scale_hz,
            "velocity_estimation_method": None,
            "execution_variant": None,
        }


@dataclass(frozen=True)
class InvocationSpec:
    name: str
    module: str
    directory: Path
    output_manager: object
    input_metadata: dict
    input_client: object
    settings: InvocationSettings
    initial_result: str
    threads: int
    memory_account: object
    acknowledgements: object
    native_parallelism: bool = False
    invocation_id: PipelineInvocationId | None = None
    product_inputs: dict | None = None
    owner_id: str | None = None


@dataclass(frozen=True)
class InvocationResult:
    """Private accepted products and changed state bindings from one worker."""

    records: dict
    state: dict
    artifacts: tuple[Path, ...]
    artifact_root: Path
    messages: tuple[str, ...] = ()
    invocation_id: PipelineInvocationId | None = None
    commit_id: str | None = None
    products: tuple = ()
    metrics: dict | None = None


def emit(kind, *, stage=None, **detail):
    sink = _SINK.get()
    if sink is not None:
        sink(kind, stage, detail)
