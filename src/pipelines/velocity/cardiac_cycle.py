"""Workflow-specific cardiac-cycle detection from raw retinal RMS-frequency signals."""

from __future__ import annotations

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    CardiacCycleAnalysis,
    cardiac_cycles_from_available_vessel,
)
from input_output.schema import RetinalSourceData

from .signal_processing import DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ


def detect_source_cardiac_cycles(source: RetinalSourceData) -> tuple[CardiacCycleAnalysis, str]:
    """Detect this workflow's cycles from raw RMS frequency before estimation."""
    import numpy as np

    from .estimation import (
        iter_velocity_estimator_chunks,
        resolve_velocity_estimator_inputs,
    )

    inputs = resolve_velocity_estimator_inputs(
        source.image_maps,
        velocity_estimation_method=source.velocity_estimation_method,
        band_ratio_frequency_scale_hz=source.band_ratio_frequency_scale_hz,
    )
    vessels = source.segmentation.vessels
    mask = vessels.artery | vessels.vein
    signals = [np.full(inputs.first_volume.shape[0], np.nan, dtype=np.float32) for _ in range(2)]
    for chunk in iter_velocity_estimator_chunks(
        inputs,
        vessel_mask=mask,
        neighborhood_mask=~mask,
    ):
        for signal, vessel_mask in zip(signals, (vessels.artery, vessels.vein)):
            if np.any(vessel_mask):
                values = chunk.rms_frequency[:, vessel_mask]
                if np.any(np.isfinite(values)):
                    signal[chunk.frame_slice] = np.nanmean(values, axis=1, dtype=np.float64)
    return cardiac_cycles_from_available_vessel(
        *signals,
        dt_seconds=float(source.holodoppler.timing.dt_seconds),
        lowpass_freq_hz=DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ,
    )


__all__ = ["detect_source_cardiac_cycles"]
