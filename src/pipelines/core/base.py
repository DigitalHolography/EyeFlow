from __future__ import annotations

import csv
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

import h5py

from dependency_utils import find_missing_dependencies

# Metadata registry populated by coded pipeline modules when they are imported.
PIPELINE_REGISTRY: dict[str, type[ProcessPipeline]] = {}


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
        # metadata for the class
        cls.name = name
        cls.description = description or getattr(cls, "description", "")
        cls.requires = required_deps or []
        cls.dag_requires = _pipeline_keys(
            dag_requires if dag_requires is not None else getattr(cls, "dag_requires", ())
        )
        cls.dag_produces = _pipeline_keys(
            dag_produces if dag_produces is not None else getattr(cls, "dag_produces", ())
        )

        missing = find_missing_dependencies(cls.requires)
        cls.missing_deps = missing
        cls.available = len(missing) == 0

        # Add to registry
        PIPELINE_REGISTRY[name] = cls
        return cls

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


# +==========================================================================+ #
# |                            PIPELINES CLASSES                             | #
# +==========================================================================+ #


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
    pipeline_cls: type[ProcessPipeline] | None = None
    error_msg: str = ""

    def instantiate(self) -> ProcessPipeline:
        """Factory method to create the actual pipeline instance."""
        if not self.available or self.pipeline_cls is None:
            return MissingPipeline(
                self.name,
                self.error_msg or self.description,
                self.missing_deps,
                self.requires,
            )
        return self.pipeline_cls()


class ProcessPipeline:
    name: str
    description: str
    available: bool
    missing_deps: list[str]
    requires: list[str]
    dag_requires: tuple[str, ...] = ()
    dag_produces: tuple[str, ...] = ()
    input_slot: str = "both"

    def __init__(self) -> None:
        # Derive the pipeline name from the module filename (e.g., basic_stats.py -> basic_stats).
        if not getattr(self, "name", None):
            module_name = (self.__class__.__module__ or "").rsplit(".", 1)[-1]
            self.name: str = module_name or self.__class__.__name__

    def run(self, h5file: h5py.File) -> ProcessResult:
        raise NotImplementedError

    def export(self, result: ProcessResult, output_path: str) -> str:
        """Default CSV export for metrics."""
        with open(output_path, "w", newline="", encoding="utf-8") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(["metric", "value"])
            for key, value in result.metrics.items():
                writer.writerow([key, value])
        return output_path


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
