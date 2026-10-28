"""Reduce HD spectra over the full frame and optional gap support in one pass."""

from __future__ import annotations

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from input_output.schema.holodoppler import HD_SPECTOGRAMS_PATH
from utils.logger import Logger

from .preprocessing import prepare_gap_spectrograms
from .spatial_reduction import DEFAULT_READ_BUFFER_BYTES, SpatialReducer, spectral_batches


def preprocess_spectrograms(
    spectrograms, analysis_zones: np.ndarray | None, *, read_buffer_bytes=DEFAULT_READ_BUFFER_BYTES,
    workers=None,
) -> dict[str, object]:
    """Pack full-frame spatial mean/median and optional gap median band power.

    Read each HD frame once in bounded sequential slabs. Any nonfinite included
    pixel makes its reduction NaN; finite negative reductions are clamped only after aggregation.
    Empty gap support skips gap products while full-frame products remain.
    """
    if spectrograms is None:
        Logger.log_warning(
            f"Skipping spectrogram analysis: HD /{HD_SPECTOGRAMS_PATH} "
            "is absent or is a legacy group rather than a frequency-bin video."
        )
        return {}
    if spectrograms.ndim != 4 or any(size == 0 for size in spectrograms.shape):
        raise ValueError("Spectrograms must have nonempty (time, freq, y, x) axes.")
    if np.dtype(spectrograms.dtype).kind not in "iuf":
        raise TypeError("Spectrograms must contain real numeric band-power values.")
    prepared = (
        prepare_gap_spectrograms(spectrograms, analysis_zones)
        if analysis_zones is not None else None
    )
    average = np.full(spectrograms.shape[:2], np.nan, dtype=np.float64)
    median = np.full_like(average, np.nan)
    gap_median = np.full_like(average, np.nan) if prepared is not None else None
    gap_indices = np.flatnonzero(prepared[0]) if prepared is not None else None
    with SpatialReducer(workers=workers) as reducer:
        for start, block in spectral_batches(spectrograms, read_buffer_bytes=read_buffer_bytes):
            stop = start + block.shape[0]
            average[start:stop], median[start:stop], gap_values = reducer(block, gap_indices)
            if gap_median is not None:
                gap_median[start:stop] = gap_values
    np.maximum(average, 0.0, out=average)
    np.maximum(median, 0.0, out=median)
    schema = EyeFlowOutputPaths.active()
    source_attrs = getattr(spectrograms, "attrs", {})
    attrs = {
        "unit": source_attrs.get("unit", "a.u."),
        "unit_source": "HD_dataset" if "unit" in source_attrs else "uncalibrated_HD_power",
        "dimDesc": ["time", "freq"],
        "source_path": f"/{HD_SPECTOGRAMS_PATH}",
        "spatial_support": "all_exported_frame_pixels",
        "spatial_pixel_count": int(np.prod(spectrograms.shape[-2:])),
        "spectrogram_spatial_transform": "none_HD_export_order",
        "spatial_reduction_coordinates": "HD_export_yx",
        "spectrogram_read_buffer_bytes": read_buffer_bytes,
        "spatial_reduction_max_workers": reducer.workers,
        "spatial_sample_dtype": str(np.dtype(spectrograms.dtype)),
        "mean_accumulation_dtype": "float64",
        "median_method": "inplace_partition_promote_only_central_values_to_float64",
        "frequency_bin_order": "HD_export_order",
        "negative_value_policy": "clamp_spatial_reduction_to_zero_after_aggregation",
        "nonfinite_value_policy": "NaN_if_any_included_sample_is_nonfinite",
    }
    metrics = {
        schema.spectrograms_spatial_average: (average, {
            **attrs, "quantity": "spatial_mean_band_power",
            "definition": "max(0, arithmetic spatial mean over all exported pixels)",
        }),
        schema.spectrograms_spatial_median: (median, {
            **attrs, "quantity": "spatial_median_band_power",
            "definition": "max(0, spatial median over all exported pixels)",
        }),
    }
    if prepared is not None:
        np.maximum(gap_median, 0.0, out=gap_median)
        metrics[schema.spectrograms_masked] = (gap_median, {
            **prepared[1], "spectrogram_read_buffer_bytes": read_buffer_bytes,
            "spatial_reduction_max_workers": reducer.workers,
        })
    return metrics
