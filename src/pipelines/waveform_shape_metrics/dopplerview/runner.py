from calculations.dopplerview_analysis import (
    ArterialWaveformAnalysisStep,
    VesselVelocityEstimatorStep,
)
from pipeline_engine.imports import HolodopplerTiming

from .constants import (
    LEGACY_FILTER_VELOCITY_SIGNALS,
    LEGACY_VELOCITY_SIGNAL_LOWPASS_HZ,
)
from .models import DopplerViewStepContext


def run_dopplerview_analysis(
    source_data,
) -> dict[str, object]:
    timing: HolodopplerTiming = source_data.timing
    step_context = DopplerViewStepContext(
        cache=source_data.dopplerview_cache(),
        holodoppler_config={
            "sampling_freq": timing.sampling_freq,
            "batch_stride": timing.batch_stride,
        },
        dopplerview_config={
            "VelocityEstimation": {
                "LocalBackgroundDist": source_data.local_background_dist,
            },
            "PulseAnalysis": {
                "FilterSignals": LEGACY_FILTER_VELOCITY_SIGNALS,
                "LowpassFreqHz": LEGACY_VELOCITY_SIGNAL_LOWPASS_HZ,
            },
        },
    )
    VesselVelocityEstimatorStep().run(step_context)
    ArterialWaveformAnalysisStep().run(step_context)
    return step_context.cache
