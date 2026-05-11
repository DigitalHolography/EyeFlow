import importlib

from pipeline_engine import (
    PIPELINE_REGISTRY,
    MissingPipeline,
    PipelineDAG,
    PipelineDescriptor,
    PipelineExecutionPlan,
    ProcessPipeline,
    ProcessResult,
    pipeline,
    registerPipeline,
)

# Add new pipeline modules here. The catalog is explicit so helper files are
# not imported as runnable pipelines by accident.
PIPELINE_MODULES = (
    "dual_input_tutorial",
    "package_pipeline_tutorial",
    "test1",
    "test2",
    "waveform_shape_metrics",
)

_PIPELINE_IMPORT_ERRORS: list[PipelineDescriptor] = []


def _import_error_descriptor(name: str, exc: BaseException) -> PipelineDescriptor:
    return PipelineDescriptor(
        name=name,
        description=f"Import Error: {exc}",
        available=False,
        error_msg=str(exc),
    )


def _import_pipeline_modules() -> None:
    for module_name in PIPELINE_MODULES:
        try:
            importlib.import_module(f"{__name__}.{module_name}")
        except Exception as exc:  # noqa: BLE001
            _PIPELINE_IMPORT_ERRORS.append(_import_error_descriptor(module_name, exc))


_import_pipeline_modules()


def load_pipeline_catalog() -> tuple[
    list[PipelineDescriptor], list[PipelineDescriptor]
]:
    """Return (available, missing) coded pipelines for UI/CLI surfaces."""
    available: list[PipelineDescriptor] = []
    missing: list[PipelineDescriptor] = list(_PIPELINE_IMPORT_ERRORS)

    for descriptor in PIPELINE_REGISTRY.values():
        if descriptor.available:
            available.append(descriptor)
        else:
            missing.append(descriptor)

    available.sort(key=lambda item: item.name.lower())
    missing.sort(key=lambda item: item.name.lower())
    return available, missing


__all__ = [
    "MissingPipeline",
    "PIPELINE_MODULES",
    "PipelineDAG",
    "PipelineDescriptor",
    "PipelineExecutionPlan",
    "ProcessPipeline",
    "ProcessResult",
    "load_pipeline_catalog",
    "pipeline",
    "registerPipeline",
]
