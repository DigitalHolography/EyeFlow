"""Immutable executable stage definitions. No runtime collaborators live here."""

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field

from data_manager.handles import FrozenMapping


def _freeze(value):
    if isinstance(value, Mapping):
        return FrozenMapping({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_freeze(item) for item in value)
    if value is None or isinstance(value, (str, bytes, int, float, complex, bool)):
        return value
    raise TypeError("Graph metadata must contain immutable scalar/container values.")


@dataclass(frozen=True, slots=True)
class ProductSpec:
    id: str
    path: str | None = None
    attrs: Mapping = field(default_factory=dict)
    pack: Callable | None = None
    validate: Callable | None = None
    optional: bool = False
    exposed: bool = False

    def __post_init__(self):
        object.__setattr__(self, "attrs", _freeze(self.attrs or {}))


@dataclass(frozen=True, slots=True)
class StageSpec:
    id: str
    run: Callable
    requires: tuple[str, ...] = ()
    products: tuple[ProductSpec, ...] = ()
    optional_requires: tuple[str, ...] = ()
    after: tuple[str, ...] = ()
    artifact: bool = False
    wait_for_artifacts: bool = False

    def __post_init__(self):
        for key in ("requires", "products", "optional_requires", "after"):
            object.__setattr__(self, key, tuple(getattr(self, key)))


@dataclass(frozen=True, slots=True)
class StageGraph:
    stages: tuple[StageSpec, ...]
    cleanup: Callable | None = None

    def __post_init__(self):
        object.__setattr__(self, "stages", tuple(self.stages))

    def plan(self, *, external_products=()):
        from .planner import StageGraphPlanner

        return StageGraphPlanner().plan(self, external_products=external_products)
