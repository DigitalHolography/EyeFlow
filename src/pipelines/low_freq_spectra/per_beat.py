"""Phase-align gap spectra without filling invalid temporal samples."""

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


def early_late_phase_windows(phase_samples) -> dict[str, tuple[np.ndarray, list[float]]]:
    """Return closed early/late windows on normalized cardiac phase coordinates.

    A sample at one third belongs to both windows. Compare in the source dtype
    so a float32 coordinate representing one third keeps that boundary behavior.
    """
    phase = np.asarray(phase_samples)
    if phase.dtype.kind != "f":
        phase = phase.astype(np.float64)
    if (
        phase.ndim != 1 or phase.size == 0 or not np.all(np.isfinite(phase))
        or np.any(phase < 0) or np.any(phase >= 1) or np.any(np.diff(phase) <= 0)
    ):
        raise ValueError("Phase samples must be strictly increasing coordinates in [0, 1).")
    split = phase.dtype.type(1 / 3)
    return {"early": (phase <= split, [0.0, 1 / 3]), "late": (phase >= split, [1 / 3, 1.0])}


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
