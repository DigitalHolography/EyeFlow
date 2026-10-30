"""Descriptive gap-spectrum QC; no filtering, rejection, or SE/SEI modification."""

from __future__ import annotations

import numpy as np

from calculations.math import nanmedian
from input_output.schema import EyeFlowOutputPaths

from .metrics import _valley_neighborhood_indices, analyze_phase_spectrum, analyze_spectrum

_QC_SORT_CHUNK_VALUES = 500_000


def pack_gap_qc(spatial_spectrum, *, source_attrs, gap_outputs):
    """Pack QC from reduced acquisition-time spectra and this run's gap products.

    Use only reduced arrays, never reread the HD spatial video. All indicators
    are descriptive. Missing support produces NaN, not a passing QC decision.
    """
    schema = EyeFlowOutputPaths.active()
    paths = schema.spectrograms_masked_qc
    values = np.asarray(spatial_spectrum, dtype=np.float64).copy()
    values[~np.isfinite(values)] = np.nan
    output = {}

    def put(family, name, data, definition, *, dims=(), unit="1", source=None, **attrs):
        output[f"{family}/{name}"] = (data, {
            "definition": definition, "dimDesc": list(dims), "unit": unit,
            "source_path": "/" + (source or schema.spectrograms_masked),
            "qc_policy": "descriptive_only_no_automatic_rejection",
            "nonfinite_value_policy": "missing_or_undefined_indicators_remain_NaN",
            **attrs,
        })

    completeness = paths.data_completeness
    put(completeness, "AcquisitionValidFraction", np.isfinite(values).mean(),
        "fraction of finite acquisition-time by frequency spatial-median samples")

    # Normalize each bin by its own positive temporal median. Unlike a
    # MAD-based z score this remains defined for a constant positive baseline.
    scale = nanmedian(values, axis=0, dtype=np.float64)
    relative_change = np.full_like(values, np.nan)
    eligible = np.isfinite(scale) & (scale > 0)
    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        relative_change[1:, eligible] = np.abs(np.diff(values[:, eligible], axis=0)) / scale[eligible]
    relative_change[~np.isfinite(relative_change)] = np.nan
    motion_attrs = {
        "normalization": "per_bin_positive_acquisition_time_median_power",
        "time_assignment": "change_from_previous_frame_first_frame_undefined",
        "interpretation": "candidate_broadband_transient_indicator_not_a_saccade_detector",
    }
    put(paths.potential_motion_contamination, "RelativeBroadbandChange",
        nanmedian(relative_change, axis=1, dtype=np.float64),
        "median across bins of absolute consecutive-frame power change divided by each bin's temporal median",
        dims=("time",), **motion_attrs)
    shape = np.asarray(source_attrs.get("mask_application_shape", []), dtype=int)
    count = source_attrs.get("gap_pixel_count", np.nan)
    total = int(np.prod(shape)) if shape.size == 2 and np.all(shape > 0) else 0
    mask_source = source_attrs.get("mask_source_path", "/" + schema.segmentation.gaps).lstrip("/")
    put(paths.spatial_support_quality, "RetainedFieldFraction", count / total if total else np.nan,
        "retained application-mask pixels divided by the entire reduced rectangular field",
        source=mask_source, coverage_denominator="entire_rectangular_spectrogram_field",
        anatomical_suitability="visual_review_required_no_automatic_anatomical_score")

    raw_entry = gap_outputs.get(schema.spectrograms_masked_per_beat)
    if raw_entry is None:
        return output

    per_beat, phase_attrs = raw_entry
    raw = np.asarray(per_beat, dtype=np.float64).copy()
    raw[~np.isfinite(raw)] = np.nan
    phase = phase_attrs["phase_samples"]
    phase_source = schema.spectrograms_masked_per_beat
    phase_meta = {"phase_samples": phase, "phase_endpoint": "excluded"}
    se_entry = gap_outputs.get(f"{schema.spectrograms_masked_metrics_root}/SE")
    put(completeness, "ValidSEPhaseFraction", np.isfinite(se_entry[0]).mean() if se_entry else np.nan,
        "fraction of registered phases having finite SE",
        source=f"{schema.spectrograms_masked_metrics_root}/SE", se_available=se_entry is not None)

    if raw.shape[1] < 2:
        return output
    peak_entry = gap_outputs.get(f"{schema.spectrograms_masked_metrics_root}/index_max1")
    metric_attrs = peak_entry[1] if peak_entry else {}
    valley_length = metric_attrs.get("valley_neighborhood_length")
    bins = _valley_neighborhood_indices(valley_length, raw.shape[1])
    relative_mad = valley_repeatability(raw, bins)[3]
    put(paths.valley_repeatability, "RelativeMAD", relative_mad,
        "unscaled valley-height MAD divided by absolute median height at the same phase",
        dims=("phase",), source=phase_source, minimum_beats_for_dispersion=2,
        **phase_meta, valley_neighborhood_indices=bins)

    _pack_beat_sensitivity(raw, gap_outputs, metric_attrs, valley_length, put, schema)
    return output


def valley_repeatability(per_beat, bins):
    """Return per-beat heights, phase median, MAD and relative MAD for QC/plots."""
    raw = np.asarray(per_beat, dtype=np.float64).copy()
    raw[~np.isfinite(raw)] = np.nan
    heights = nanmedian(raw[:, bins, :], axis=1, dtype=np.float64)
    median = nanmedian(heights, axis=1, dtype=np.float64)
    mad = nanmedian(np.abs(heights - median[:, None]), axis=1, dtype=np.float64)
    mad[np.isfinite(heights).sum(axis=1) < 2] = np.nan
    relative_mad = np.full_like(mad, np.nan)
    np.divide(mad, np.abs(median), out=relative_mad, where=np.isfinite(median) & (median != 0))
    return heights, median, mad, relative_mad


def _pack_beat_sensitivity(raw, gap_outputs, metric_attrs, valley_length, put, schema):
    """Recompute the complete raw-spectrum SEI estimator with each beat omitted."""
    family = schema.spectrograms_masked_qc.sei_beat_sensitivity
    root = schema.spectrograms_masked_metrics_root
    entry = gap_outputs.get(f"{root}/SEI")
    reference = float(entry[0]) if entry else np.nan
    beat_count = raw.shape[2]
    estimates = np.full(beat_count, np.nan)
    phase_counts = np.zeros(beat_count, dtype=int)
    can_evaluate = entry is not None and beat_count >= 2
    if can_evaluate:
        calibration = dict(
            sampling_frequency_hz=metric_attrs["sampling_frequency_hz"],
            fft_bin_count=metric_attrs.get("fft_bin_count"),
            binning_factor=metric_attrs.get("frequency_binning_factor", 1),
            frequency_start=metric_attrs.get("spectrogram_frequency_start", 0.0),
            frequency_end=metric_attrs.get("spectrogram_frequency_end", 0.5),
            valley_neighborhood_length=valley_length,
        )
        for beat, (pooled, phase_spectrum) in enumerate(_leave_one_beat_medians(raw)):
            fixed = analyze_spectrum(pooled, **calibration)
            curve = analyze_phase_spectrum(
                phase_spectrum, index_max1=fixed["index_max1"], index_max2=fixed["index_max2"],
                **calibration,
            )["SE"]
            phase_counts[beat] = np.isfinite(curve).sum()
            estimates[beat] = np.float32(nanmedian(curve, dtype=np.float64))
    delta = estimates - reference
    attrs = {
        "source": schema.spectrograms_masked_per_beat,
        "reference_source_path": "/" + root + "/SEI",
        "omitted_beat_index_base": 0,
        "remaining_beat_count": max(0, beat_count - 1),
        "peak_selection": "reselect_pooled_peaks_for_each_omission",
        "valid_phase_count_per_omission": phase_counts,
        "evaluation_status": "evaluated" if can_evaluate else "requires_SEI_and_at_least_two_beats",
        "interpretation": "influence_diagnostic_not_a_confidence_interval_or_motion_classification",
    }
    put(family, "DeltaSEI", delta, "SEI without the indexed beat minus full-recording SEI",
        dims=("omitted_beat",), **attrs)


def _leave_one_beat_medians(raw):
    """Reuse sorted spectra for the pooled and phase medians of each omission.

    The production reductions use float32 medians. Sorting once avoids a
    separate NaN-aware sort of the same samples for every omitted beat. The
    caller has already normalized nonfinite input samples to NaN.
    """
    values = np.asarray(raw, dtype=np.float32)
    phase_count, frequency_count, beat_count = values.shape
    if np.any(np.isfinite(raw) & ~np.isfinite(values)):
        # Preserve nanmedian's treatment of finite float64 inputs that
        # overflow to infinity when it applies its float32 policy.
        for beat in range(beat_count):
            remaining = raw[:, :, np.arange(beat_count) != beat]
            yield nanmedian(remaining, axis=(0, 2)), nanmedian(remaining, axis=2)
        return
    pooled_by_beat = np.empty((beat_count, frequency_count), dtype=np.float32)
    phase_by_beat = np.empty((beat_count, phase_count, frequency_count), dtype=np.float32)
    chunk_width = max(1, _QC_SORT_CHUNK_VALUES // (phase_count * beat_count))

    for first in range(0, frequency_count, chunk_width):
        stop = min(first + chunk_width, frequency_count)
        chunk = values[:, first:stop, :]
        finite = np.isfinite(chunk)
        width = stop - first

        # The beat-axis order gives the location of the omitted sample, so
        # either middle rank can be read without sorting the other beats again.
        phase_order = np.argsort(np.where(finite, chunk, np.inf), axis=2, kind="stable")
        sorted_phase = np.take_along_axis(chunk, phase_order, axis=2)
        phase_ranks = np.argsort(phase_order, axis=2)
        phase_counts = finite.sum(axis=2)

        flattened = np.where(finite, chunk, np.inf).transpose(1, 0, 2).reshape(width, -1)
        pooled_order = np.argsort(flattened, axis=1, kind="stable")
        sorted_pooled = np.take_along_axis(flattened, pooled_order, axis=1)
        beat_ids = np.broadcast_to(
            np.arange(beat_count), (width, phase_count, beat_count)
        ).reshape(width, -1)
        sorted_beat_ids = np.take_along_axis(beat_ids, pooled_order, axis=1)
        pooled_valid = np.isfinite(sorted_pooled)
        pooled_counts = pooled_valid.sum(axis=1)
        frequency_rows = np.arange(width)

        for beat in range(beat_count):
            omitted_finite = finite[:, :, beat]
            remaining_phase_counts = phase_counts - omitted_finite
            omitted_rank = phase_ranks[:, :, beat]
            lower_rank = (remaining_phase_counts - 1) // 2
            upper_rank = remaining_phase_counts // 2
            lower_index = lower_rank + (omitted_finite & (omitted_rank <= lower_rank))
            upper_index = upper_rank + (omitted_finite & (omitted_rank <= upper_rank))
            lower = np.take_along_axis(sorted_phase, lower_index[:, :, None], axis=2)[:, :, 0]
            upper = np.take_along_axis(sorted_phase, upper_index[:, :, None], axis=2)[:, :, 0]
            phase_spectrum = np.median(np.stack((lower, upper), axis=-1), axis=-1)
            phase_spectrum[remaining_phase_counts == 0] = np.nan
            phase_by_beat[beat, :, first:stop] = phase_spectrum

            # A cumulative count selects the pooled middle ranks after the
            # omitted beat's phase samples have been removed.
            keep = pooled_valid & (sorted_beat_ids != beat)
            cumulative = np.cumsum(keep, axis=1)
            remaining_pooled_counts = pooled_counts - omitted_finite.sum(axis=0)
            lower_rank = (remaining_pooled_counts - 1) // 2 + 1
            upper_rank = remaining_pooled_counts // 2 + 1
            lower_index = np.count_nonzero(cumulative < lower_rank[:, None], axis=1)
            upper_index = np.count_nonzero(cumulative < upper_rank[:, None], axis=1)
            upper_index = np.minimum(upper_index, sorted_pooled.shape[1] - 1)
            lower = sorted_pooled[frequency_rows, lower_index]
            upper = sorted_pooled[frequency_rows, upper_index]
            pooled = np.median(np.stack((lower, upper), axis=-1), axis=-1)
            pooled[remaining_pooled_counts == 0] = np.nan
            pooled_by_beat[beat, first:stop] = pooled

    for beat in range(beat_count):
        yield pooled_by_beat[beat], phase_by_beat[beat]
