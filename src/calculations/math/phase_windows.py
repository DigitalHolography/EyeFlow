"""Reusable closed phase windows, NaN-aware reductions, and late-minus-early deltas."""

from __future__ import annotations

import numpy as np

from .statistics import nanmedian

DEFAULT_EARLY_PHASE_INTERVAL = (0.0, 1 / 3)
DEFAULT_LATE_PHASE_INTERVAL = (1 / 3, 1.0)


def early_late_phase_windows(
    phase_samples, *, early_interval=DEFAULT_EARLY_PHASE_INTERVAL,
    late_interval=DEFAULT_LATE_PHASE_INTERVAL,
) -> dict[str, tuple[np.ndarray, list[float]]]:
    """Select closed early/late intervals expressed as fractions of a cardiac cycle.

    Preserve source-dtype comparisons so float32 coordinates at a shared boundary
    belong to both windows. Phase one remains an excluded source endpoint.
    """
    phase = np.asarray(phase_samples)
    if phase.dtype.kind != "f":
        phase = phase.astype(np.float64)
    if (
        phase.ndim != 1 or phase.size == 0 or not np.all(np.isfinite(phase))
        or np.any(phase < 0) or np.any(phase >= 1) or np.any(np.diff(phase) <= 0)
    ):
        raise ValueError("Phase samples must be strictly increasing coordinates in [0, 1).")
    windows = {}
    for name, interval in (("early", early_interval), ("late", late_interval)):
        bounds = np.asarray(interval, dtype=np.float64)
        if (
            bounds.shape != (2,) or not np.all(np.isfinite(bounds))
            or not 0 <= bounds[0] <= bounds[1] <= 1
        ):
            raise ValueError("Phase intervals must contain two ordered finite bounds in [0, 1].")
        lower, upper = (phase.dtype.type(bound) for bound in bounds)
        windows[name] = ((phase >= lower) & (phase <= upper), bounds.tolist())
    return windows


def phase_window_medians(values, windows, *, axis: int = 0) -> dict[str, np.ndarray]:
    """Reduce any named phase windows, ignoring NaNs and retaining missing cells."""
    values = np.asarray(values, dtype=np.float64)
    if values.ndim < 1:
        raise ValueError("Phase reductions require an array with a phase axis.")
    values = np.moveaxis(values, axis, 0)
    result = {}
    for name, (mask, _interval) in windows.items():
        mask = np.asarray(mask)
        if mask.dtype != np.bool_ or mask.shape != values.shape[:1]:
            raise ValueError("Phase masks must be Boolean and match the phase axis.")
        result[name] = nanmedian(values[mask], axis=0, dtype=np.float64)
    return result


def late_minus_early(early, late):
    """Return the signed late-minus-early difference, propagating missing values."""
    return np.asarray(late, dtype=np.float64) - np.asarray(early, dtype=np.float64)
