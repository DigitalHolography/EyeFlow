"""Lazy public exports for independent scientific calculation packages."""

from importlib import import_module

_EXPORTS = {
    "CardiacCycleAnalysis": ".signal_analysis.cardiac_cycle",
    "SpectralCardiacCycleAnalysis": ".signal_analysis.cardiac_cycle",
    "SystoleDetectionResult": ".signal_analysis.cardiac_cycle",
    "analyze_cardiac_cycles": ".signal_analysis.cardiac_cycle",
    "cardiac_cycles_from_available_vessel": ".signal_analysis.cardiac_cycle",
    "find_systole_index": ".signal_analysis.cardiac_cycle",
    "missing_vessel_cardiac_cycles": ".signal_analysis.cardiac_cycle",
    "spectral_cardiac_cycle_analysis": ".signal_analysis.cardiac_cycle",
    "PerBeatAnalysisInput": ".signal_analysis.per_beat.runner",
    "PerBeatAnalysisResult": ".signal_analysis.per_beat.runner",
    "run_per_beat_analysis": ".signal_analysis.per_beat.runner",
    "PerBeatSegmentAnalysisResult": ".signal_analysis.per_beat.segments",
    "aggregate_per_beat_segment_analysis": ".signal_analysis.per_beat.segments",
    "per_beat_segment_analysis": ".signal_analysis.per_beat.segments",
    "PerBeatSignalAnalysisResult": ".signal_analysis.per_beat.signal",
    "per_beat_signal_analysis": ".signal_analysis.per_beat.signal",
    "ArterialWaveformAnalysis": ".signal_analysis.waveform",
    "PairedVesselCycles": ".signal_analysis.waveform",
    "PulseMetricData": ".signal_analysis.waveform",
    "VenousWaveformAnalysis": ".signal_analysis.waveform",
    "arterial_waveform_analysis": ".signal_analysis.waveform",
    "average_cycle": ".signal_analysis.waveform",
    "cycle_extrema": ".signal_analysis.waveform",
    "paired_vessel_cycles": ".signal_analysis.waveform",
    "pulse_metric": ".signal_analysis.waveform",
    "pulse_metric_from_signal": ".signal_analysis.waveform",
    "venous_waveform_analysis": ".signal_analysis.waveform"
}


def __getattr__(name):
    module = _EXPORTS.get(name)
    if module is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(module, __name__), name)
    globals()[name] = value
    return value


__all__ = [
    "ArterialWaveformAnalysis",
    "CardiacCycleAnalysis",
    "PairedVesselCycles",
    "PerBeatAnalysisInput",
    "PerBeatAnalysisResult",
    "PerBeatSegmentAnalysisResult",
    "PerBeatSignalAnalysisResult",
    "aggregate_per_beat_segment_analysis",
    "PulseMetricData",
    "SpectralCardiacCycleAnalysis",
    "SystoleDetectionResult",
    "VenousWaveformAnalysis",
    "arterial_waveform_analysis",
    "average_cycle",
    "cycle_extrema",
    "find_systole_index",
    "analyze_cardiac_cycles",
    "cardiac_cycles_from_available_vessel",
    "paired_vessel_cycles",
    "per_beat_segment_analysis",
    "per_beat_signal_analysis",
    "pulse_metric",
    "pulse_metric_from_signal",
    "missing_vessel_cardiac_cycles",
    "run_per_beat_analysis",
    "spectral_cardiac_cycle_analysis",
    "venous_waveform_analysis",
]
