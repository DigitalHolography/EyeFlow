"""Phase-align gap spectra without filling invalid temporal samples."""

from __future__ import annotations

import numpy as np

from calculations.math.cycle_boundaries import normalize_cycle_boundaries

DEFAULT_PHASE_SAMPLE_COUNT = 128


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
