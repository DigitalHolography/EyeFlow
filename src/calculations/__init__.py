"""Pure EyeFlow scientific calculations."""

def __getattr__(name):
    if name == "ArterialWaveformAnalysisStep":
        from .retinal_velocity import ArterialWaveformAnalysisStep
        return ArterialWaveformAnalysisStep
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

__all__ = [
    "ArterialWaveformAnalysisStep",
]
