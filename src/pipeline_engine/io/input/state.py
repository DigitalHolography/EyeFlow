"""Input control-thread state and scheduled operations; never scientific products."""

from dataclasses import dataclass, replace
from enum import Enum
from multiprocessing.shared_memory import SharedMemory
from typing import ClassVar

import numpy as np

SourceKey = tuple[str, str]
Selection = slice | tuple[slice, ...] | None


@dataclass(frozen=True)
class SmallRead:
    key: SourceKey
    kind: ClassVar[str] = "small"

    @property
    def identity(self):
        return self.kind, self.key


@dataclass(frozen=True)
class ReplayRead:
    key: SourceKey
    region: int
    kind: ClassVar[str] = "replay"

    @property
    def identity(self):
        return self.kind, self.key, self.region


@dataclass(frozen=True)
class DirectRead:
    owner: int
    key: SourceKey
    selection: tuple[slice, ...]
    kind: ClassVar[str] = "direct"

    @property
    def identity(self):
        return self.kind, self.owner


@dataclass(frozen=True)
class StreamRead:
    stream: int
    buffer: int
    kind: ClassVar[str] = "stream"

    @property
    def identity(self):
        return self.kind, self.stream, self.buffer


InputOperation = SmallRead | ReplayRead | DirectRead | StreamRead


@dataclass(frozen=True)
class ScheduledRead:
    operation: InputOperation
    priority: int
    sequence: int
    queued_at: float

    @property
    def order(self):
        return self.priority, self.sequence

    def promote(self, priority):
        return replace(self, priority=min(priority, self.priority))


@dataclass(kw_only=True)
class InFlightRead:
    token: int
    job: ScheduledRead
    issued_priority: int
    reserved: int = 0
    operation_type: ClassVar[type | None] = None

    def __post_init__(self):
        if self.operation_type is None:
            raise TypeError("Use an operation-specific in-flight record")
        if not isinstance(self.job.operation, self.operation_type):
            raise TypeError(f"{type(self).__name__} requires {self.operation_type.__name__}")


@dataclass(kw_only=True)
class SmallInFlight(InFlightRead):
    operation_type: ClassVar[type] = SmallRead


@dataclass(kw_only=True)
class DirectInFlight(InFlightRead):
    operation_type: ClassVar[type] = DirectRead
    shape: tuple[int, ...]
    started: float


@dataclass(kw_only=True)
class StreamInFlight(InFlightRead):
    operation_type: ClassVar[type] = StreamRead
    replay_key: SourceKey | None
    regions: range | tuple[int, ...]


@dataclass(kw_only=True)
class ReplayInFlight(InFlightRead):
    operation_type: ClassVar[type] = ReplayRead
    working_bytes: int
    staged: float
    cancel_followup: bool = False


class BufferState(str, Enum):
    FREE = "free"
    FILLING = "filling"
    READY = "ready"
    IN_USE = "in_use"
    FAILED = "failed"


@dataclass
class StreamBuffer:
    memory: SharedMemory
    state: BufferState = BufferState.FREE
    start: int | None = None
    stop: int | None = None
    speculative: bool = False
    consumed: bool = False

    def begin_fill(self, start, stop):
        if self.state is not BufferState.FREE:
            raise RuntimeError("Only a free stream buffer may be filled")
        self.start, self.stop, self.state = start, stop, BufferState.FILLING

    def ready(self, *, speculative):
        if self.state is not BufferState.FILLING:
            raise RuntimeError("Only a filling stream buffer may become ready")
        self.state, self.speculative, self.consumed = BufferState.READY, speculative, False

    def take(self):
        if self.state is not BufferState.READY:
            raise RuntimeError("Only a ready stream buffer may be leased")
        self.state, self.consumed = BufferState.IN_USE, True

    def release(self):
        if self.state is not BufferState.IN_USE:
            raise RuntimeError("Streaming buffer released without a lease")
        self.state = BufferState.FREE


@dataclass
class StreamState:
    owner: int | None
    key: SourceKey
    buffers: list[StreamBuffer]
    frames: int
    shape: tuple[int, ...]
    dtype: np.dtype
    reserved: int
    total: int
    replay: bool
    priority: int
    next_fill: int = 0
    next_take: int = 0
    waiting: float | None = None
    error: str | None = None
    closing: bool = False

    def adopt(self, owner):
        if self.owner is not None or self.closing:
            raise RuntimeError("Only an unowned active stream may be adopted")
        self.owner, self.priority = owner, 1


class TransferState(Enum):
    IDLE = "idle"
    LOADING = "loading"
    LEASED = "leased"


@dataclass
class DirectTransfer:
    memory: SharedMemory
    size: int
    state: TransferState = TransferState.IDLE
    token: int = 0
    closing: bool = False

    @property
    def busy(self):
        return self.state is not TransferState.IDLE

    @property
    def loading(self):
        return self.state is TransferState.LOADING

    def begin_read(self, token):
        if self.busy or self.closing:
            raise RuntimeError("Release the current direct transfer before reading")
        self.token, self.state = token, TransferState.LOADING

    def read_complete(self):
        if not self.loading:
            raise RuntimeError("Direct transfer completed without an active read")
        self.state = TransferState.LEASED

    def release(self, token):
        if token == self.token and self.state is TransferState.LEASED:
            self.state = TransferState.IDLE


class RegionState(str, Enum):
    NOT_REQUESTED = "not_requested"
    QUEUED = "queued"
    LOADING = "loading"
    READY = "ready"
    FAILED = "failed"
