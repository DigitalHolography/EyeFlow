from __future__ import annotations

import importlib
import pkgutil
import sys
from pathlib import Path

from app_settings import runtime_pipelines_path
from pipeline_engine import (
    PIPELINE_REGISTRY,
    PipelineDescriptor,
)

_IGNORED_MODULES = {"utils"}
_BASE_PIPELINE_PATHS = tuple(Path(path).resolve() for path in __path__)


def _pipeline_search_paths() -> list[Path]:
    paths: list[Path] = []
    for path in (runtime_pipelines_path(), *_BASE_PIPELINE_PATHS):
        try:
            resolved = path.resolve()
        except OSError:
            resolved = path
        if not resolved.is_dir() or resolved in paths:
            continue
        paths.append(resolved)
    return paths


def _refresh_pipeline_search_path() -> None:
    __path__[:] = [str(path) for path in _pipeline_search_paths()]


def _discover_pipelines() -> tuple[list[PipelineDescriptor], list[PipelineDescriptor]]:
    available: list[PipelineDescriptor] = []
    missing: list[PipelineDescriptor] = []

    PIPELINE_REGISTRY.clear()
    _refresh_pipeline_search_path()
    importlib.invalidate_caches()

    discovered_modules: set[str] = set()
    for module_info in pkgutil.iter_modules(__path__):
        if module_info.name in _IGNORED_MODULES or module_info.name.startswith("_"):
            continue
        if module_info.name in discovered_modules:
            continue
        discovered_modules.add(module_info.name)

        module_name = f"{__name__}.{module_info.name}"
        try:
            if module_name in sys.modules:
                importlib.reload(sys.modules[module_name])
            else:
                importlib.import_module(module_name)
        except Exception as exc:  # noqa: BLE001
            missing.append(
                PipelineDescriptor(
                    name=module_info.name,
                    description=f"Import Error: {exc}",
                    available=False,
                    error_msg=str(exc),
                )
            )

    for descriptor in PIPELINE_REGISTRY.values():
        if descriptor.available:
            available.append(descriptor)
        else:
            missing.append(descriptor)

    available.sort(key=lambda item: item.name.lower())
    missing.sort(key=lambda item: item.name.lower())
    return available, missing


def load_pipeline_catalog() -> tuple[
    list[PipelineDescriptor], list[PipelineDescriptor]
]:
    """Return (available, missing) coded pipelines for UI/CLI surfaces."""
    return _discover_pipelines()


__all__ = [
    "PipelineDescriptor",
    "load_pipeline_catalog",
]
