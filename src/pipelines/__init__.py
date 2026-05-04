from .core.base import (
    MissingPipeline,
    PipelineDescriptor,
    ProcessPipeline,
    ProcessResult,
    registerPipeline,
)

_PIPELINE_CLASSES: list[type[ProcessPipeline]] = []
_PIPELINE_IMPORT_ERRORS: list[PipelineDescriptor] = []


def _import_error_descriptor(name: str, exc: BaseException) -> PipelineDescriptor:
    return PipelineDescriptor(
        name=name,
        description=f"Import Error: {exc}",
        available=False,
        error_msg=str(exc),
    )


try:
    from .dual_input_tutorial import DualInputTutorial
except Exception as exc:  # noqa: BLE001
    _PIPELINE_IMPORT_ERRORS.append(_import_error_descriptor("dual_input_tutorial", exc))
else:
    _PIPELINE_CLASSES.append(DualInputTutorial)

try:
    from .waveform_shape_metrics import WaveformShapeMetrics
except Exception as exc:  # noqa: BLE001
    _PIPELINE_IMPORT_ERRORS.append(
        _import_error_descriptor("waveform_shape_metrics", exc)
    )
else:
    _PIPELINE_CLASSES.append(WaveformShapeMetrics)


def _descriptor_from_class(cls: type[ProcessPipeline]) -> PipelineDescriptor:
    return PipelineDescriptor(
        name=getattr(cls, "name", cls.__name__),
        description=getattr(cls, "description", ""),
        available=bool(getattr(cls, "available", True)),
        input_slot=getattr(cls, "input_slot", "both"),
        requires=list(getattr(cls, "requires", [])),
        missing_deps=list(getattr(cls, "missing_deps", [])),
        pipeline_cls=cls,
    )


def write_result_h5(*args, **kwargs):
    from input_output import write_result_h5 as _write_result_h5

    return _write_result_h5(*args, **kwargs)


def write_combined_results_h5(*args, **kwargs):
    from input_output import write_combined_results_h5 as _write_combined_results_h5

    return _write_combined_results_h5(*args, **kwargs)


def load_pipeline_catalog() -> tuple[
    list[PipelineDescriptor], list[PipelineDescriptor]
]:
    """Return (available, missing) coded pipelines for UI/CLI surfaces."""
    available: list[PipelineDescriptor] = []
    missing: list[PipelineDescriptor] = list(_PIPELINE_IMPORT_ERRORS)

    for cls in _PIPELINE_CLASSES:
        descriptor = _descriptor_from_class(cls)
        if descriptor.available:
            available.append(descriptor)
        else:
            missing.append(descriptor)

    available.sort(key=lambda p: p.name.lower())
    missing.sort(key=lambda p: p.name.lower())
    return available, missing


__all__ = [
    "ProcessPipeline",
    "ProcessResult",
    "registerPipeline",
    "write_result_h5",
    "write_combined_results_h5",
    "load_pipeline_catalog",
    "MissingPipeline",
    *[cls.__name__ for cls in _PIPELINE_CLASSES],
]
