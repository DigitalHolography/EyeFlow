"""Track exactly which pipeline invocation owns each newly reserved artifact."""

from contextlib import contextmanager
from contextvars import ContextVar
from pathlib import Path
from threading import RLock

_CURRENT = ContextVar("eyeflow_artifact_owner", default=None)


def current_artifact_transaction():
    transaction = _CURRENT.get()
    return getattr(transaction, "current", transaction)


@contextmanager
def artifact_transaction(transaction):
    token = _CURRENT.set(transaction)
    try:
        yield transaction
    finally:
        _CURRENT.reset(token)


def register_artifact_path(path, *, writing=False):
    """Writers claim writes; OutputManager also records new path reservations."""
    transaction = _CURRENT.get()
    if transaction is not None:
        transaction.reserve(path, writing=writing)


class ArtifactOwnership:
    """Run-local ownership shared by computations and artifact publication."""

    def __init__(self):
        self._owners = {}
        self._lock = RLock()

    def transaction(self, root):
        return ArtifactTransaction(self, Path(root).resolve())

    @property
    def paths(self):
        with self._lock:
            return [path for path in self._owners if path.is_file()]


class ArtifactTransaction:
    def __init__(self, registry, root):
        self._registry = registry
        self._root = root
        self._created = set()

    def reserve(self, path, *, writing=False):
        path = Path(path).resolve()
        if not path.is_relative_to(self._root):
            return
        with self._registry._lock:
            owner = self._registry._owners.get(path)
            if owner is not None and owner is not self:
                if writing:
                    raise FileExistsError(f"Artifact belongs to another pipeline invocation: {path}")
                return
            if owner is None:
                if path.exists():
                    if writing and path.is_file():
                        raise FileExistsError(f"Refusing to overwrite an existing unowned artifact: {path}")
                    return
                self._registry._owners[path] = self
                self._created.add(path)

    def rollback(self):
        """Remove only this invocation's files; never infer ownership from time."""
        with self._registry._lock:
            for path in self._created:
                if path.is_file():
                    path.unlink(missing_ok=True)
                self._registry._owners.pop(path, None)
            self._created.clear()


class ArtifactScope:
    """Rotate owners at a handoff; committed exports retain their own owner."""

    def __init__(self, registry, root):
        self.registry = registry
        self.root = root
        self.current = registry.transaction(root)

    def reserve(self, path, *, writing=False):
        self.current.reserve(path, writing=writing)

    def commit(self):
        self.current = self.registry.transaction(self.root)

    def rollback(self):
        self.current.rollback()
