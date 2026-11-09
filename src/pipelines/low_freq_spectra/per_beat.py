"""Phase-align reduced and spatial Doppler power without filling invalid samples."""

from __future__ import annotations

import warnings

import numpy as np

from calculations.math import nanmedian
from calculations.math.cycle_boundaries import normalize_cycle_boundaries

DEFAULT_PHASE_SAMPLE_COUNT = 128


def beat_spectrum_std(per_beat, phase_mask) -> tuple[np.ndarray, np.ndarray]:
    """Return population beat STD after a per-beat phase median, plus valid beat counts.

    Input axes are (phase, freq, beat). Ignore NaNs in both reductions. A
    frequency with no valid beat spectra stays NaN; one valid beat has STD zero.
    """
    values = np.asarray(per_beat, dtype=np.float64)
    mask = np.asarray(phase_mask)
    if values.ndim != 3 or mask.shape != values.shape[:1] or mask.dtype != np.bool_:
        raise ValueError("Beat spectra require (phase, freq, beat) and a Boolean phase mask.")
    spectra = nanmedian(values[mask], axis=0, dtype=np.float64)
    valid_counts = np.count_nonzero(~np.isnan(spectra), axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        std = np.nanstd(spectra, axis=1, ddof=0, dtype=np.float64)
    return std.astype(np.float32), valid_counts


def gap_spectrograms_per_beat(
    values: np.ndarray,
    cycle_boundaries: np.ndarray,
    *,
    sample_count: int = DEFAULT_PHASE_SAMPLE_COUNT,
) -> np.ndarray:
    """Return (phase, freq, beat) power using zero-based moments boundaries.

    Each interval includes its next boundary as an interpolation anchor. The
    output phases are [0, 1), so the shared endpoint is not duplicated. Linear
    interpolation retains NaN gaps instead of bridging over invalid samples.
    Incomplete intervals outside the supplied boundary sequence are omitted.
    """
    values = np.asarray(values, dtype=np.float64).copy()
    if values.ndim != 2 or any(size == 0 for size in values.shape):
        raise ValueError("Gap spectra must have nonempty (time, freq) axes.")
    if not isinstance(sample_count, int) or sample_count < 2:
        raise ValueError("sample_count must be an integer of at least two.")
    values[~np.isfinite(values)] = np.nan
    boundaries = normalize_cycle_boundaries(cycle_boundaries, values.shape[0], index_base=0)
    phase = np.arange(sample_count, dtype=np.float64) / sample_count
    output = np.full((sample_count, values.shape[1], boundaries.size - 1), np.nan, np.float32)
    for beat, (start, stop) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        signal = values[int(start):int(stop) + 1]
        source_phase = np.linspace(0.0, 1.0, signal.shape[0])
        for frequency in range(values.shape[1]):
            output[:, frequency, beat] = np.interp(phase, source_phase, signal[:, frequency])
    return output


def align_power(values, boundaries, sample_count):
    """Linear phase interpolation, including next boundary, without NaN bridging.

    Vectorized spatial equivalent of gap_spectrograms_per_beat; retain float64
    power until the final output packing. Exact anchors do not touch neighbors.
    """
    values = np.array(values, dtype=np.float64, copy=True)
    values[~np.isfinite(values)] = np.nan
    phase = np.arange(sample_count) / sample_count
    aligned = np.empty((sample_count, len(boundaries) - 1, values.shape[1]))
    for beat, (start, stop) in enumerate(zip(boundaries[:-1], boundaries[1:])):
        source_phase = np.linspace(0, 1, int(stop - start) + 1)
        lower = np.searchsorted(source_phase, phase, side="right") - 1
        upper = np.minimum(lower + 1, source_phase.size - 1)
        weight = (phase - source_phase[lower]) / (source_phase[upper] - source_phase[lower])
        left = values[start + lower]
        right = values[start + upper]
        with np.errstate(invalid="ignore", over="ignore"):
            result = left + weight[:, None] * (right - left)
        exact = weight == 0
        result[exact] = left[exact]
        result[~np.isfinite(result)] = np.nan
        aligned[:, beat] = result
    return aligned
