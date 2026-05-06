from __future__ import annotations

from calculations.steps import ArterialWaveformAnalysisStep, VesselVelocityEstimatorStep
from input_output import DOPPLER_VIEW_SCHEMA
from input_output.input_access import HolodopplerTiming

from .models import DopplerViewStepContext
from .source_inputs import dopplerview_cache_from_h5, local_background_dist


def run_dopplerview_analysis(
    ctx,
    timing: HolodopplerTiming,
) -> dict[str, object]:
    step_context = DopplerViewStepContext(
        cache=dopplerview_cache_from_h5(ctx),
        holodoppler_config={
            "sampling_freq": timing.sampling_freq,
            "batch_stride": timing.batch_stride,
        },
        dopplerview_config={
            DOPPLER_VIEW_SCHEMA.config_value("local_background_dist").section: {
                DOPPLER_VIEW_SCHEMA.config_value(
                    "local_background_dist"
                ).json_key: local_background_dist(ctx),
            }
        },
    )
    VesselVelocityEstimatorStep().run(step_context)
    ArterialWaveformAnalysisStep().run(step_context)
    return step_context.cache
