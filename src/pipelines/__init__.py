import importlib

from pipeline_engine import (
    PIPELINE_REGISTRY,
    PipelineDescriptor,
)

# Add new pipeline modules here. The catalog is explicit so helper files are
# not imported as runnable pipelines by accident.
PIPELINE_MODULES = (
    "dual_input_tutorial",
    "matlab_pulse_poc",
    "package_pipeline_tutorial",
    "test1",
    "test2",
    "waveform_shape_metrics",
)

_PIPELINE_IMPORT_ERRORS: list[PipelineDescriptor] = []

for module_name in PIPELINE_MODULES:
    try:
        importlib.import_module(f"{__name__}.{module_name}")
    except Exception as exc:  # noqa: BLE001
        _PIPELINE_IMPORT_ERRORS.append(
            PipelineDescriptor(
                name=module_name,
                description=f"Import Error: {exc}",
                available=False,
                error_msg=str(exc),
            )
        )


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
    "PIPELINE_MODULES",
    "PipelineDescriptor",
    "load_pipeline_catalog",
]
