"""Pipeline access to canonical HD/DV retinal source data."""

from dataclasses import replace

import numpy as np

from app_settings import validate_velocity_estimation_method
from input_output.readers.canonical import (
    align_mask, align_spectrogram_frame, apply_mask_alignment,
    apply_optional_mask_alignment, apply_spatial_alignment,
    load_retinal_source_data as _read_retinal_source_data,
)
from retinal_source_data import RetinalSourceData
from velocity_calibration import DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ


def load_retinal_source_data(
    hd,
    dv,
    *,
    velocity_estimation_method: str = "doppler_moments",
    include_spectrograms: bool = True,
    band_ratio_frequency_scale_hz: float = DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
) -> RetinalSourceData:
    """Apply the analysis policy for an unavailable optic-disc measurement."""
    source = _read_retinal_source_data(
        hd,
        dv,
        velocity_estimation_method=velocity_estimation_method,
        include_spectrograms=include_spectrograms,
        band_ratio_frequency_scale_hz=band_ratio_frequency_scale_hz,
    )
    segmentation = source.segmentation
    if not segmentation.optic_disc.is_fallback:
        return source
    # Keep the native vein support for velocity background inpainting, but
    # disable venous measurement when DV supplied no usable optic disc.
    vessels = replace(
        segmentation.vessels,
        vein=np.zeros_like(segmentation.vessels.vein, dtype=bool),
    )
    return replace(source, segmentation=replace(segmentation, vessels=vessels))


def load_vessel_topology_inputs(ctx) -> RetinalSourceData:
    """Resolve the method-specific source while preserving shared topology fallback."""
    ctx.require_inputs("hd", "dv")
    cached = ctx.state.get("velocity_source")
    if isinstance(cached, RetinalSourceData) and (
        getattr(ctx, "velocity_estimation_method", None) in {None, cached.velocity_estimation_method}
    ):
        return cached
    hd = ctx.inputs.hd.as_holodoppler()
    dv = ctx.inputs.dv.as_dopplerview()
    method = getattr(ctx, "velocity_estimation_method", None)
    scale = getattr(ctx, "band_ratio_frequency_scale_hz", DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ)
    if method is None:
        failures = {}
        for candidate in ("doppler_moments", "frequency_bands"):
            try:
                return load_retinal_source_data(
                    hd, dv, velocity_estimation_method=candidate,
                    include_spectrograms=False, band_ratio_frequency_scale_hz=scale,
                )
            except (KeyError, ValueError, TypeError) as exc:
                failures[candidate] = str(exc)
        raise ValueError(f"No usable retinal source for shared topology: {failures}")
    return load_retinal_source_data(
        hd, dv, velocity_estimation_method=validate_velocity_estimation_method(method),
        include_spectrograms=False, band_ratio_frequency_scale_hz=scale,
    )


__all__ = [
    "align_mask", "align_spectrogram_frame", "apply_mask_alignment",
    "apply_optional_mask_alignment", "apply_spatial_alignment",
    "load_retinal_source_data", "load_vessel_topology_inputs",
]
