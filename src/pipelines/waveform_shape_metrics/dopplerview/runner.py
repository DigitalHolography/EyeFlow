from calculations.dopplerview_analysis import (
    ArterialWaveformAnalysisStep,
    run_chunked_velocity_estimator,
)
from pipeline_engine.imports import HolodopplerTiming

from .constants import (
    LEGACY_FILTER_VELOCITY_SIGNALS,
    LEGACY_VELOCITY_SIGNAL_LOWPASS_HZ,
)
from .models import DopplerViewStepContext


def run_dopplerview_analysis(
    source_data,
    scratch_h5,
    *,
    log=None,
) -> dict[str, object]:
    timing: HolodopplerTiming = source_data.timing
    cache = run_chunked_velocity_estimator(
        moment0=source_data.moment0,
        moment2=source_data.moment2,
        moment0_is_flat_fielded=source_data.moment0_is_flat_fielded,
        moment2_is_flat_fielded=source_data.moment2_is_flat_fielded,
        artery_mask=source_data.retinal_artery_mask,
        vein_mask=source_data.retinal_vein_mask,
        local_background_dist=source_data.local_background_dist,
        flat_field_gaussian_ratio=source_data.flat_field_gaussian_ratio,
        flat_field_border=source_data.flat_field_border,
        scratch_h5=scratch_h5,
        log=log,
    )
    step_context = DopplerViewStepContext(
        cache=cache,
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
    ArterialWaveformAnalysisStep().run(step_context)
    return step_context.cache
