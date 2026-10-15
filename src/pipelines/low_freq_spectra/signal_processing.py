"""Temporal Fourier and concatenated-beat SVD filtering of gap band power."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from calculations.math.bandlimiting import (
    DEFAULT_BAND_LIMITED_HARMONICS,
    band_limit,
    validate_harmonics,
)
from calculations.math.low_rank import SvdFilterResult, remove_leading_svd_modes

DEFAULT_GAP_SVD_REMOVED_MODES = 2


@dataclass(frozen=True)
class BandLimitedResult:
    filtered: np.ndarray
    filtered_series: np.ndarray
    highest_harmonic: np.ndarray
    harmonics: int


def band_limit_gap_spectrograms(
    per_beat: np.ndarray,
    beat_durations_seconds: np.ndarray,
    *,
    harmonics: int = DEFAULT_BAND_LIMITED_HARMONICS,
) -> BandLimitedResult:
    """Keep DC and the first ``harmonics`` within each periodic beat.

    One phase sample spans duration / phase_count seconds for that beat.
    Filter only fully finite phase series; retain incomplete series unchanged
    rather than filling gaps or propagating NaNs through every phase sample.
    Use the shared signed Fourier reconstruction without taking abs or
    clipping undershoots, preserving the per-beat DC power.
    """
    values = _per_beat_values(per_beat)
    durations = np.asarray(beat_durations_seconds, dtype=np.float64)
    if durations.shape != (values.shape[2],) or not np.all(
        np.isfinite(durations) & (durations > 0)
    ):
        raise ValueError("Each beat must have one positive finite duration in seconds.")
    harmonics = validate_harmonics(harmonics)
    eligible = np.all(np.isfinite(values), axis=0)
    highest = np.full(values.shape[2], min(harmonics, values.shape[0] // 2), dtype=np.int64)
    for beat in range(values.shape[2]):
        columns = eligible[:, beat]
        if not columns.any():
            continue
        values[:, columns, beat] = band_limit(
            values[:, columns, beat],
            harmonics=harmonics,
            reconstruction="signed",
            axis=0,
        )
    return BandLimitedResult(values, eligible, highest, harmonics)


def svd_filter_gap_spectrograms(
    per_beat: np.ndarray,
    *,
    mode_count: int | None = None,
) -> tuple[np.ndarray, SvdFilterResult]:
    """Append beats chronologically and return filtered (beat*phase, freq) power.

    Concatenation is in normalized phase time, not uniform acquisition time.
    This does not affect SVD, which uses the matrix without a sampling rate.
    """
    values = _per_beat_values(per_beat)
    phases, frequencies, beats = values.shape
    matrix = values.transpose(2, 0, 1).reshape(beats * phases, frequencies)
    result = remove_leading_svd_modes(
        matrix,
        mode_count=DEFAULT_GAP_SVD_REMOVED_MODES if mode_count is None else mode_count,
    )
    return result.filtered, result


def _per_beat_values(per_beat):
    values = np.asarray(per_beat, dtype=np.float32).copy()
    if values.ndim != 3 or any(size == 0 for size in values.shape):
        raise ValueError("Gap spectra must have nonempty (phase, freq, beat) axes.")
    values[~np.isfinite(values)] = np.nan
    return values
