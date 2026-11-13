"""Small spawn-safe output commands; scientific payloads stay in owned files."""

from dataclasses import dataclass, field
import os


def _rebuild_connection(handle, readable, writable):
    if os.name == 'nt':
        from multiprocessing.connection import PipeConnection as kind
    else:
        from multiprocessing.connection import Connection as kind
    return kind(handle, readable=readable, writable=writable)


@dataclass(frozen=True)
class WriterConnection:
    """Transfer ownership directly to the spawned process.

    Windows' ordinary pipe reducer leaves a duplicate in the parent for the
    child to steal. Interrupted startup can leak that duplicate. The spawning
    process already supports direct child duplication, as used by SemLock.
    """

    connection: object

    def __reduce__(self):
        from multiprocessing.context import get_spawning_popen
        spawning = get_spawning_popen()
        if spawning is None:
            raise RuntimeError('Writer connections may only be transferred during spawn')
        connection = self.connection
        return (_rebuild_connection, (spawning.duplicate_for_child(connection.fileno()),
                                     connection.readable, connection.writable))


@dataclass(frozen=True)
class OutputPolicy:
    max_pending_jobs: int = 32
    max_pending_bytes: int = 256 * 2**20
    admission_memory_bytes: int = 16 * 2**20
    hdf5_cache_bytes: int = 2**20
    shutdown_timeout: float = 5.
    operation_timeout: float = 300.

    def __post_init__(self):
        for name in ("max_pending_jobs", "max_pending_bytes", "admission_memory_bytes", "hdf5_cache_bytes"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer")
        if self.shutdown_timeout <= 0 or self.operation_timeout <= 0:
            raise ValueError("Output timeouts must be positive")


@dataclass(frozen=True)
class WriteCommand:
    sequence: int
    manifest: str
    digest: str
    final: bool = False
    accepted_at: float = 0.


@dataclass(frozen=True)
class WriteAcknowledgement:
    sequence: int
    error: str | None = None
    metrics: dict = field(default_factory=dict)
    paths: tuple = ()
    messages: tuple = ()


class OutputServiceError(OSError):
    """Accepted output could not be materialized completely."""
