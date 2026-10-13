"""Temporal Fourier and concatenated-beat SVD filtering of gap band power."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from calculations.math import harmonic_pack
from calculations.math.low_rank import SvdFilterResult, remove_leading_svd_modes

DEFAULT_GAP_LOWPASS_HZ = 15.0
DEFAULT_GAP_SVD_REMOVED_MODES = 2


@dataclass(frozen=True)
class BandLimitedResult:
    filtered: np.ndarray
    filtered_series: np.ndarray
    highest_harmonic: np.ndarray
    cutoff_hz: float


def band_limit_gap_spectrograms(
    per_beat: np.ndarray,
    beat_durations_seconds: np.ndarray,
    *,
    cutoff_hz: float = DEFAULT_GAP_LOWPASS_HZ,
) -> BandLimitedResult:
    """Keep DC and temporal frequencies <= cutoff_hz within each periodic beat.

    One phase sample spans duration / phase_count seconds for that beat.
    Filter only fully finite phase series; retain incomplete series unchanged
    rather than filling gaps or propagating NaNs through every phase sample.
    Use the shared normalized real FFT reconstruction without taking abs or
    clipping undershoots, preserving the per-beat DC power.
    """
    values = _per_beat_values(per_beat)
    durations = np.asarray(beat_durations_seconds, dtype=np.float64)
    if durations.shape != (values.shape[2],) or not np.all(np.isfinite(durations) & (durations > 0)):
        raise ValueError("Each beat must have one positive finite duration in seconds.")
    if not np.isfinite(cutoff_hz) or cutoff_hz < 0:
        raise ValueError("cutoff_hz must be finite and non-negative.")
    eligible = np.all(np.isfinite(values), axis=0)
    highest = np.minimum(
        np.floor(cutoff_hz * durations + 1e-12).astype(np.int64), values.shape[0] // 2,
    )
    for beat in range(values.shape[2]):
        columns = eligible[:, beat]
        if not columns.any():
            continue
        pack = harmonic_pack(values[:, columns, beat], int(highest[beat]), axis=0)
        if pack["vb"] is not None:
            values[:, columns, beat] = pack["vb"]
    return BandLimitedResult(values, eligible, highest, float(cutoff_hz))


def svd_filter_gap_spectrograms(
    per_beat: np.ndarray, *, mode_count: int | None = None,
) -> tuple[np.ndarray, SvdFilterResult]:
    """Append beats chronologically and return filtered (beat*phase, freq) power.

    Concatenation is in normalized phase time, not uniform acquisition time.
    This does not affect SVD, which uses the matrix without a sampling rate.
    """
    values = _per_beat_values(per_beat)
    phases, frequencies, beats = values.shape
    matrix = values.transpose(2, 0, 1).reshape(beats * phases, frequencies)
    result = remove_leading_svd_modes(
        matrix, mode_count=DEFAULT_GAP_SVD_REMOVED_MODES if mode_count is None else mode_count,
    )
    return result.filtered, result


def _per_beat_values(per_beat):
    values = np.asarray(per_beat, dtype=np.float32).copy()
    if values.ndim != 3 or any(size == 0 for size in values.shape):
        raise ValueError("Gap spectra must have nonempty (phase, freq, beat) axes.")
    values[~np.isfinite(values)] = np.nan
    return values
