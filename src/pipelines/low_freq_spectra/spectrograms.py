"""Reduce HD spectra over the full frame and optional gap support in one pass."""

from __future__ import annotations

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from input_output.schema.holodoppler import HD_SPECTOGRAMS_PATH
from pipelines.vessel_inputs import align_spectrogram_frame
from utils.logger import Logger

from .preprocessing import prepare_gap_spectrograms


def preprocess_spectrograms(spectrograms, analysis_zones: np.ndarray | None) -> dict[str, object]:
    """Pack full-frame spatial mean/median and optional gap mean band power.

    Read each HD frame once. Any nonfinite included pixel makes its reduction
    NaN; finite negative reductions are clamped only after aggregation.
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
    gap_average = np.full_like(average, np.nan) if prepared is not None else None
    for frame_index in range(spectrograms.shape[0]):
        frame = np.asarray(align_spectrogram_frame(spectrograms[frame_index]), dtype=np.float64)
        full = frame.reshape(frame.shape[0], -1)
        valid = np.all(np.isfinite(full), axis=1)
        average[frame_index, valid] = np.mean(full[valid], axis=1)
        median[frame_index, valid] = np.median(full[valid], axis=1)
        if prepared is not None:
            samples = frame[:, prepared[0]]
            gap_valid = np.all(np.isfinite(samples), axis=1)
            gap_average[frame_index, gap_valid] = np.mean(samples[gap_valid], axis=1)
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
        "spectrogram_spatial_transform": "flip_D2_then_swap_D2_D3",
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
        np.maximum(gap_average, 0.0, out=gap_average)
        metrics[schema.spectrograms_masked] = (gap_average, prepared[1])
    return metrics
