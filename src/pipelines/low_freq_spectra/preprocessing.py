"""Reduce the exported HD spectra to spatial median band power in gap pixels."""

from __future__ import annotations

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from input_output.schema.holodoppler import HD_SPECTOGRAMS_PATH
from pipelines.vessel_inputs import align_spectrogram_frame
from utils.logger import Logger

from .masks import gap_velocity_mask


def spatial_median_band_power(spectrograms, mask: np.ndarray) -> np.ndarray:
    """Take the median of included pixels for every exported time/frequency pair.

    The mask uses the persisted EF lower-left (x, y) frame. Read only one
    time frame at a time from the lazy (time, freq, y, x) source, aligning
    its spatial coordinates at the canonical input boundary. Excluded pixels,
    including excluded NaN/Inf values, never enter the median.
    Finite negative spatial medians from HD's resampling are clamped to zero
    after reduction; individual pixel values are not clamped. Any nonfinite
    included sample yields NaN for the affected time/frequency pair. Callers must
    skip empty gap support.
    """
    if spectrograms.ndim != 4:
        raise ValueError("Spectrograms must have (time, freq, height, width) axes.")
    if np.dtype(spectrograms.dtype).kind not in "iuf":
        raise TypeError("Spectrograms must contain real numeric band-power values.")
    mask = np.asarray(mask)
    if mask.ndim != 2 or mask.dtype != np.bool_:
        raise TypeError("The application mask must be a 2-D boolean array.")
    if mask.shape != tuple(spectrograms.shape[-2:][::-1]):
        raise ValueError("The (x, y) application mask must match the aligned spectrogram shape.")

    result = np.full(spectrograms.shape[:2], np.nan, dtype=np.float64)
    if not np.any(mask):
        raise ValueError("An empty gap mask must be skipped before reducing spectra.")

    for time_index in range(spectrograms.shape[0]):
        samples = align_spectrogram_frame(spectrograms[time_index])[:, mask]
        valid = np.all(np.isfinite(samples), axis=-1)
        # Use float64 for the midpoint of an even number of included samples.
        result[time_index, valid] = np.median(
            np.asarray(samples[valid], dtype=np.float64), axis=-1,
        )
    np.maximum(result, 0.0, out=result)
    return result


def preprocess_gap_spectrograms(spectrograms, analysis_zones: np.ndarray) -> dict[str, object]:
    """Pack gap median power when valid spatial support is available."""
    prepared = prepare_gap_spectrograms(spectrograms, analysis_zones)
    if prepared is None:
        return {}
    mask, attrs = prepared
    schema = EyeFlowOutputPaths.active()
    return {schema.spectrograms_masked: (spatial_median_band_power(spectrograms, mask), attrs)}


def prepare_gap_spectrograms(spectrograms, analysis_zones: np.ndarray):
    """Prepare aligned support and provenance; None signals skipped gap analysis.

    Keep segmentation and shared topology available even if native or reduced
    gap support is empty. No spectral samples are read in either empty case.
    """
    schema = EyeFlowOutputPaths.active()
    if not np.any(analysis_zones):
        Logger.log_warning("Skipping gap analysis: the native gap mask is empty.")
        return None
    if spectrograms is None:
        Logger.log_warning(
            f"Omitting /{schema.spectrograms_masked}: HD /{HD_SPECTOGRAMS_PATH} "
            "is absent or is a legacy group rather than a frequency-bin video."
        )
        return None

    native_mask = gap_velocity_mask(analysis_zones, tuple(spectrograms.shape[-2:]))
    pixel_count = int(np.count_nonzero(native_mask))
    if pixel_count == 0:
        Logger.log_warning("Skipping gap analysis: the reduced application gap mask is empty.")
        return None
    # Match /Segmentation/Gaps/value's published coordinates at the reduced resolution.
    mask = np.flip(native_mask, axis=0).T
    source_attrs = getattr(spectrograms, "attrs", {})
    attrs = {
        "unit": source_attrs.get("unit", "a.u."),
        "unit_source": "HD_dataset" if "unit" in source_attrs else "uncalibrated_HD_power",
        "dimDesc": ["time", "freq"],
        "quantity": "spatial_median_band_power",
        "definition": "max(0, spatial median of exported band-power values over gap pixels)",
        "spatial_aggregation_method": "median",
        "source_path": f"/{HD_SPECTOGRAMS_PATH}",
        "mask_source_path": f"/{schema.segmentation.gaps}",
        "mask_application_coordinates": "EF_lower_left_xy",
        "mask_native_shape": analysis_zones.shape,
        "mask_application_shape": mask.shape,
        "mask_application_dimDesc": ["x", "y"],
        "spectrogram_spatial_transform": "flip_D2_then_swap_D2_D3",
        "mask_reduction": (
            "identity"
            if analysis_zones.shape == native_mask.shape
            else "nearest_neighbor_interpolation"
        ),
        "mask_interpolation_order": 0,
        "mask_interpolation_grid": "pixel_centers",
        "mask_target_shape_source": f"/{HD_SPECTOGRAMS_PATH} last two dimensions (height, width)",
        "gap_pixel_count": pixel_count,
        "frequency_bin_order": "HD_export_order",
        "negative_value_policy": "clamp_spatial_median_to_zero",
        "nonfinite_value_policy": "NaN_if_any_included_sample_is_nonfinite",
        "empty_mask_policy": "skip_gap_analysis",
    }
    return mask, attrs
