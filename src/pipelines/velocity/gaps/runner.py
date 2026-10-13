"""Adapt moments-workflow timing to durable beat-resolved gap spectra."""

from __future__ import annotations

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from ..models import RetinalVelocity
from utils.logger import Logger

from .per_beat import DEFAULT_PHASE_SAMPLE_COUNT, gap_spectrograms_per_beat
from .signal_processing import (
    DEFAULT_GAP_SVD_REMOVED_MODES,
    band_limit_gap_spectrograms,
    svd_filter_gap_spectrograms,
)


def run_gap_analysis(ctx, moments: RetinalVelocity | None) -> dict[str, object]:
    """Process gaps once within velocity using its explicit moments result."""
    schema = EyeFlowOutputPaths.active()
    source = ctx.output.h5.get(schema.spectrograms_masked)
    if source is None:
        # The shared preprocessing stage already warned about empty/missing support.
        return {}
    if moments is None:
        Logger.log_warning(
            "Skipping per-beat gap analysis: Doppler moments timing is unavailable."
        )
        return {}
    boundaries = moments.cycle_boundary_indexes
    if moments.cardiac_cycle_source == "none" or boundaries.size < 2:
        Logger.log_warning(
            "Skipping per-beat gap analysis: no detected Doppler moments beat intervals."
        )
        return {}
    if source.shape[0] != moments.continuous("artery", raw=True).size:
        raise ValueError("Gap spectra and Doppler moments timing must share the same frame count.")

    per_beat = gap_spectrograms_per_beat(source[()], boundaries)
    attrs = {
        **dict(source.attrs),
        "dimDesc": ["phase", "freq", "beat"],
        "source_path": "/" + schema.spectrograms_masked,
        "definition": (
            "spatial mean gap band power linearly resampled between moments beat boundaries"
        ),
        "cycle_boundary_source_path": "/" + schema.cardiac_cycle.systolic_peak_frame_indices,
        "cycle_boundary_workflow": "doppler_moments",
        "cycle_boundary_vessel_source": moments.cardiac_cycle_source,
        "cycle_boundary_indices": boundaries,
        "cycle_boundary_index_base": 0,
        "beat_duration_seconds": moments.cycle_durations_seconds,
        "phase_samples": (
            np.arange(DEFAULT_PHASE_SAMPLE_COUNT, dtype=np.float32) / DEFAULT_PHASE_SAMPLE_COUNT
        ),
        "phase_endpoint": "excluded",
        "interpolation": "linear_without_bridging_NaN",
        "incomplete_edge_intervals": "omitted",
        "nonfinite_value_policy": "preserve_NaN_without_temporal_gap_filling",
    }
    attrs.pop("nameID", None)
    metrics = {schema.spectrograms_masked_per_beat: (per_beat, attrs)}
    band_limited = band_limit_gap_spectrograms(per_beat, moments.cycle_durations_seconds)
    filtered_attrs = {
        **attrs,
        "source_path": "/" + schema.spectrograms_masked_per_beat,
        "post_filter_negative_policy": "retain_undershoots_to_preserve_DC_power",
    }
    metrics[schema.spectrograms_masked_per_beat_band_limited] = (band_limited.filtered, {
        **filtered_attrs,
        "definition": "periodic per-beat temporal FFT low-pass retaining DC through 15 Hz",
        "filter_method": "normalized_rfft_harmonic_truncation_irfft",
        "filter_axis": "phase_independently_per_frequency_and_beat",
        "filter_cutoff_hz": band_limited.cutoff_hz,
        "filter_preserve_DC": True,
        "filter_boundary": "periodic_within_each_beat",
        "filter_highest_retained_harmonic_per_beat": band_limited.highest_harmonic,
        "filter_phase_sample_interval_seconds_per_beat": (
            np.asarray(moments.cycle_durations_seconds) / per_beat.shape[0]
        ),
        "filter_incomplete_series_policy": "leave_unchanged_without_gap_filling",
        "filter_series_count": int(np.count_nonzero(band_limited.filtered_series)),
        "filter_skipped_series_count": int(np.count_nonzero(~band_limited.filtered_series)),
    })
    if not np.all(band_limited.filtered_series):
        Logger.log_warning(
            "Gap FFT band limiting left "
            f"{np.count_nonzero(~band_limited.filtered_series)} incomplete beat/bin series "
            "unchanged to preserve invalid-data gaps."
        )
    filtered, svd = svd_filter_gap_spectrograms(per_beat)
    metrics[schema.spectrograms_masked_per_beat_svd_continuous] = (filtered, {
        **filtered_attrs,
        "dimDesc": ["time", "freq"],
        "time_axis_definition": "chronological_concatenation_of_phase_normalized_beats",
        "time_sampling": "uniform_within_each_beat_variable_between_beats",
        "relative_time_seconds_definition": (
            "sum(previous_beat_durations_seconds) + phase * current_beat_duration_seconds"
        ),
        "source_beat_count": per_beat.shape[2],
        "phase_samples_per_beat": per_beat.shape[0],
        "definition": "two leading centered SVD modes removed from concatenated beats; mean restored",
        "filter_method": "concatenate_beats_then_mean_preserving_SVD_mode_removal",
        "svd_matrix_dimDesc": ["beat_then_phase", "freq"],
        "svd_matrix_shape": [per_beat.shape[0] * per_beat.shape[2], per_beat.shape[1]],
        "svd_concatenation_order": "chronological_beats_with_phase_increasing_within_each_beat",
        "svd_input_time_coordinate": "phase_normalized_not_uniform_acquisition_time",
        "svd_modes_requested": DEFAULT_GAP_SVD_REMOVED_MODES,
        "svd_modes_removed": svd.modes_removed,
        "svd_mode_order": "descending_singular_value",
        "svd_preserve_temporal_mean": True,
        "svd_temporal_mean_per_frequency": svd.temporal_mean,
        "svd_singular_values": svd.singular_values,
        "svd_filter_applied": svd.applied,
        "svd_fit_row_count": int(np.count_nonzero(svd.valid_rows)),
        "svd_fit_frequency_count": int(np.count_nonzero(svd.valid_columns)),
        "svd_nonfinite_policy": "exclude_all_invalid_bins_and_incomplete_rows_leave_excluded_values_unchanged",
        "svd_insufficient_dimensions_policy": "leave_unchanged_and_warn_preserving_at_least_one_mode",
    })
    # Phase coordinates describe each source beat, rather than the new time axis.
    svd_attrs = metrics[schema.spectrograms_masked_per_beat_svd_continuous][1]
    svd_attrs["within_beat_phase_samples"] = svd_attrs.pop("phase_samples")
    svd_attrs["within_beat_phase_endpoint"] = svd_attrs.pop("phase_endpoint")
    if not svd.applied:
        Logger.log_warning(
            "Gap SVD filtering left per-beat spectra unchanged: insufficient finite "
            "matrix dimensions to remove two modes while retaining at least one mode."
        )
    return metrics
