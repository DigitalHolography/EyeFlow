"""Reduce the exported HD spectra to spatial median band power in gap pixels."""

from __future__ import annotations

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from input_output.schema.holodoppler import HD_SPECTOGRAMS_PATH
from utils.logger import Logger

from .masks import gap_velocity_mask
from .spatial_reduction import DEFAULT_READ_BUFFER_BYTES, SpatialReducer, spectral_batches


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
    # Published masks remain in EF coordinates. Apply the equivalent native
    # support here, without flipping/transposing each spectral frame.
    mask = native_mask
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
        "mask_application_coordinates": "HD_export_yx",
        "mask_native_shape": analysis_zones.shape,
        "mask_application_shape": mask.shape,
        "mask_application_dimDesc": ["y", "x"],
        "mask_coordinate_strategy": "prepare_native_support_once_keep_published_EF_mask_unchanged",
        "spectrogram_spatial_transform": "none_HD_export_order",
        "spectrogram_read_buffer_bytes": DEFAULT_READ_BUFFER_BYTES,
        "spatial_sample_dtype": str(np.dtype(spectrograms.dtype)),
        "median_method": "inplace_partition_promote_only_central_values_to_float64",
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
