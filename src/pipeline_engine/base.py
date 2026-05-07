from __future__ import annotations

import inspect
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from typing import Any

from dependency_utils import find_missing_dependencies

# Metadata registry populated by coded pipeline modules when they are imported.
PIPELINE_REGISTRY: dict[str, PipelineDescriptor] = {}


# Decorator to attach metadata to coded pipeline classes.
def registerPipeline(
    name: str,
    description: str = "",
    required_deps: list[str] | None = None,
    *,
    dag_requires: Iterable[str] | None = None,
    dag_produces: Iterable[str] | None = None,
):
    def decorator(cls):
        cls.name = name
        cls.description = description or getattr(cls, "description", "")
        cls.requires = list(required_deps or [])
        cls.dag_requires = _pipeline_keys(
            dag_requires if dag_requires is not None else getattr(cls, "dag_requires", ())
        )
        cls.dag_produces = _pipeline_keys(
            dag_produces if dag_produces is not None else getattr(cls, "dag_produces", ())
        )

        missing = find_missing_dependencies(cls.requires)
        cls.missing_deps = missing
        cls.available = len(missing) == 0

        PIPELINE_REGISTRY[name] = PipelineDescriptor(
            name=name,
            description=cls.description,
            available=cls.available,
            input_slot=getattr(cls, "input_slot", "both"),
            requires=list(cls.requires),
            missing_deps=list(cls.missing_deps),
            dag_requires=tuple(cls.dag_requires),
            dag_produces=tuple(cls.dag_produces),
            pipeline_factory=cls,
            source_path=_source_path(cls),
        )
        return cls

    return decorator


def pipeline(
    name: str,
    description: str = "",
    requires: list[str] | None = None,
    *,
    dag_requires: Iterable[str] | None = None,
    dag_produces: Iterable[str] | None = None,
    input_slot: str = "both",
):
    """Register a function pipeline.

    The decorated function receives a PipelineContext. It can write outputs via
    ctx.write()/ctx.write_many(), return a ProcessResult, return a metrics dict,
    or return None after writing directly to the context.
    """

    def decorator(func: Callable[[Any], ProcessResult | Mapping[str, Any] | None]):
        required = list(requires or [])
        missing = find_missing_dependencies(required)
        available = len(missing) == 0
        descriptor = PipelineDescriptor(
            name=name,
            description=description or (inspect.getdoc(func) or ""),
            available=available,
            input_slot=input_slot,
            requires=required,
            missing_deps=missing,
            dag_requires=_pipeline_keys(dag_requires or ()),
            dag_produces=_pipeline_keys(dag_produces or ()),
            pipeline_factory=lambda: FunctionPipeline(
                name=name,
                description=description or (inspect.getdoc(func) or ""),
                func=func,
                input_slot=input_slot,
                requires=required,
                missing_deps=missing,
                available=available,
                dag_requires=_pipeline_keys(dag_requires or ()),
                dag_produces=_pipeline_keys(dag_produces or ()),
            ),
            source_path=_source_path(func),
        )
        PIPELINE_REGISTRY[name] = descriptor
        return func

    return decorator


def _pipeline_keys(keys: Iterable[str]) -> tuple[str, ...]:
    return tuple(str(key).strip() for key in keys if str(key).strip())


@dataclass
class ProcessResult:
    metrics: dict[str, Any]
    attrs: dict[str, Any] | None = None  # attributes stored on the pipeline group
    output_h5_path: str | None = None


@dataclass
class DatasetValue:
    """Represents a dataset payload plus optional attributes for that dataset."""

    data: Any
    attrs: dict[str, Any] | None = None


def with_attrs(data: Any, attrs: dict[str, Any]) -> DatasetValue:
    """Convenience helper to attach attributes to a dataset value."""
    return DatasetValue(data=data, attrs=attrs)


class ProcessPipeline:
    name: str
    description: str
    available: bool
    missing_deps: list[str]
    requires: list[str]
    dag_requires: tuple[str, ...] = ()
    dag_produces: tuple[str, ...] = ()
    input_slot: str = "both"
    source_path: str | None = None

    def __init__(self) -> None:
        if not getattr(self, "name", None):
            module_name = (self.__class__.__module__ or "").rsplit(".", 1)[-1]
            self.name: str = module_name or self.__class__.__name__
        self.source_path = _source_path(self.__class__)

    def run(self, ctx: Any) -> ProcessResult | Mapping[str, Any] | None:
        raise NotImplementedError


class FunctionPipeline(ProcessPipeline):
    def __init__(
        self,
        *,
        name: str,
        description: str,
        func: Callable[[Any], ProcessResult | Mapping[str, Any] | None],
        input_slot: str,
        requires: list[str],
        missing_deps: list[str],
        available: bool,
        dag_requires: tuple[str, ...],
        dag_produces: tuple[str, ...],
    ) -> None:
        self.name = name
        self.description = description
        self.func = func
        self.input_slot = input_slot
        self.requires = requires
        self.missing_deps = missing_deps
        self.available = available
        self.dag_requires = dag_requires
        self.dag_produces = dag_produces
        self.source_path = _source_path(func)

    def run(self, ctx: Any) -> ProcessResult | Mapping[str, Any] | None:
        return self.func(ctx)


@dataclass
class PipelineDescriptor:
    name: str
    description: str
    available: bool
    input_slot: str = "both"
    # To avoid Python Mutable Default Arguments
    requires: list[str] = field(default_factory=list)
    missing_deps: list[str] = field(default_factory=list)
    dag_requires: tuple[str, ...] = ()
    dag_produces: tuple[str, ...] = ()
    pipeline_factory: Callable[[], ProcessPipeline] | None = None
    error_msg: str = ""
    source_path: str | None = None

    def instantiate(self) -> ProcessPipeline:
        """Factory method to create the actual pipeline instance."""
        if not self.available or self.pipeline_factory is None:
            return MissingPipeline(
                self.name,
                self.error_msg or self.description,
                self.missing_deps,
                self.requires,
            )
        return self.pipeline_factory()


class MissingPipeline(ProcessPipeline):
    """Placeholder for pipelines whose dependencies are missing."""

    available = False

    def __init__(
        self, name: str, description: str, missing_deps: list[str], requires: list[str]
    ) -> None:
        # super().__init__()
        self.name = name
        self.description = description or "Pipeline unavailable (missing dependencies)."
        self.missing_deps = missing_deps
        self.requires = requires

    def run(self, h5file):
        missing = ", ".join(
            self.missing_deps or self.requires or ["unknown dependency"]
        )
        raise ImportError(
            f"Pipeline '{self.name}' unavailable. Missing dependencies: {missing}"
        )


def _source_path(obj: object) -> str | None:
    try:
        return inspect.getsourcefile(obj)
    except (OSError, TypeError):
        return None
