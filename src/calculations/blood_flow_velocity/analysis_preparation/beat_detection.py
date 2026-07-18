"""Robust systole detection from the derivative of the arterial waveform.

The primary pass intentionally preserves the historical detector.  A conservative
second pass only searches an approximately two-period gap for one locally
plausible missed peak; it does not lower the threshold over the whole trace.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from calculations.math import butter_lowpass_filtfilt

RECOVERY_GAP_RATIO_MIN = 1.75
RECOVERY_GAP_RATIO_MAX = 2.25
RECOVERY_SEARCH_RADIUS_PERIOD_FRACTION = 0.20
RECOVERY_ORIGINAL_THRESHOLD_FRACTION = 0.80
RECOVERY_MEDIAN_PRIMARY_HEIGHT_FRACTION = 0.50
RECOVERY_MEDIAN_PRIMARY_PROMINENCE_FRACTION = 0.50


@dataclass(frozen=True)
class SystoleDetectionResult:
    systole_indexes: np.ndarray
    artery_signal_filtered: np.ndarray
    derivative_signal: np.ndarray
    min_peak_distance: int
    min_peak_height: np.float32
    initial_systole_indexes: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int32)
    )
    recovered_systole_indexes: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int32)
    )
    nominal_period_samples: np.float32 = np.float32(np.nan)
    interval_period_ratio: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.float32)
    )
    interval_duration_valid: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.bool_)
    )


def find_systole_index(
    pulse_artery,
    *,
    dt: np.float32,
    lowpass_freq_hz: np.float32 = np.float32(15.0),
    min_duration_seconds: np.float32 = np.float32(0.5),
    validation_distance: int = 10,
    recover_missed_peaks: bool = True,
) -> SystoleDetectionResult:
    find_peaks, peak_prominences = _scipy_signal_dependencies()

    pulse = np.asarray(pulse_artery, dtype=np.float32).reshape(-1)
    filtered_pulse = butter_lowpass_filtfilt(
        pulse,
        dt_seconds=np.float32(dt),
        lowpass_freq_hz=np.float32(lowpass_freq_hz),
        order=4,
    )
    derivative = np.gradient(filtered_pulse).astype(np.float32)
    min_peak_height = np.float32(np.percentile(derivative, 95))
    min_peak_distance = _min_peak_distance(dt, min_duration_seconds)

    peaks, _ = find_peaks(
        derivative,
        height=min_peak_height,
        distance=min_peak_distance,
    )
    initial_indexes = _validate_peaks(peaks.astype(np.int32), validation_distance)
    if initial_indexes.size == 0:
        raise ValueError("No systole peaks detected. Check signal quality or parameters.")

    nominal_period = _nominal_period_samples(initial_indexes)
    recovered_indexes = np.empty(0, dtype=np.int32)
    indexes = initial_indexes
    if recover_missed_peaks:
        recovered_indexes = _recover_single_missed_peaks(
            derivative,
            initial_indexes,
            nominal_period,
            min_peak_height,
            min_peak_distance,
            find_peaks,
            peak_prominences,
        )
        if recovered_indexes.size:
            indexes = np.sort(
                np.concatenate((initial_indexes, recovered_indexes))
            ).astype(np.int32, copy=False)

    interval_period_ratio = _interval_period_ratio(indexes, nominal_period)
    return SystoleDetectionResult(
        systole_indexes=indexes,
        initial_systole_indexes=initial_indexes,
        recovered_systole_indexes=recovered_indexes,
        artery_signal_filtered=filtered_pulse,
        derivative_signal=derivative,
        min_peak_distance=min_peak_distance,
        min_peak_height=min_peak_height,
        nominal_period_samples=np.float32(nominal_period),
        interval_period_ratio=interval_period_ratio,
        interval_duration_valid=_valid_intervals(interval_period_ratio),
    )


def _min_peak_distance(dt_seconds: np.float32, min_duration_seconds: np.float32) -> int:
    if dt_seconds <= 0:
        raise ValueError("dt_seconds must be positive for systole detection.")
    return max(1, int(np.floor(float(min_duration_seconds) / float(dt_seconds))))


def _validate_peaks(peaks: np.ndarray, min_distance: int) -> np.ndarray:
    if peaks.size == 0:
        return peaks.astype(np.int32, copy=False)
    validated = [int(peaks[0])]
    for peak in peaks[1:]:
        if int(peak) - validated[-1] >= int(min_distance):
            validated.append(int(peak))
    return np.asarray(validated, dtype=np.int32)


def _nominal_period_samples(peaks: np.ndarray) -> float:
    periods = np.diff(peaks).astype(np.float64, copy=False)
    if periods.size == 0:
        return float("nan")

    nominal = float(np.median(periods))
    for _ in range(2):
        typical = periods[periods <= 1.5 * nominal]
        if typical.size < 2:
            break
        nominal = float(np.median(typical))
    return nominal


def _recover_single_missed_peaks(
    derivative: np.ndarray,
    primary_peaks: np.ndarray,
    nominal_period: float,
    min_peak_height: np.float32,
    min_peak_distance: int,
    find_peaks,
    peak_prominences,
) -> np.ndarray:
    if primary_peaks.size < 4 or not np.isfinite(nominal_period) or nominal_period <= 0:
        return np.empty(0, dtype=np.int32)

    primary_heights = derivative[primary_peaks]
    prominence_window = _prominence_window(nominal_period, derivative.size)
    primary_prominences = peak_prominences(
        derivative,
        primary_peaks,
        wlen=prominence_window,
    )[0]
    median_primary_height = float(np.median(primary_heights))
    median_primary_prominence = float(np.median(primary_prominences))

    recovered: list[int] = []
    for left, right in zip(primary_peaks[:-1], primary_peaks[1:], strict=True):
        gap_ratio = (int(right) - int(left)) / nominal_period
        if not RECOVERY_GAP_RATIO_MIN <= gap_ratio <= RECOVERY_GAP_RATIO_MAX:
            continue

        midpoint = 0.5 * (int(left) + int(right))
        radius = max(
            1,
            int(round(RECOVERY_SEARCH_RADIUS_PERIOD_FRACTION * nominal_period)),
        )
        search_start = max(int(left) + min_peak_distance, int(round(midpoint)) - radius)
        search_stop = min(int(right) - min_peak_distance, int(round(midpoint)) + radius)
        if search_stop - search_start < 2:
            continue

        local_peaks, _ = find_peaks(derivative[search_start : search_stop + 1])
        if local_peaks.size == 0:
            continue
        candidates = local_peaks.astype(np.int32, copy=False) + search_start
        candidate = int(candidates[np.argmax(derivative[candidates])])
        candidate_height = float(derivative[candidate])
        candidate_prominence = float(
            peak_prominences(
                derivative,
                np.asarray([candidate], dtype=np.int32),
                wlen=prominence_window,
            )[0][0]
        )

        if candidate_height < RECOVERY_ORIGINAL_THRESHOLD_FRACTION * float(
            min_peak_height
        ):
            continue
        if (
            candidate_height
            < RECOVERY_MEDIAN_PRIMARY_HEIGHT_FRACTION * median_primary_height
        ):
            continue
        if (
            candidate_prominence
            < RECOVERY_MEDIAN_PRIMARY_PROMINENCE_FRACTION
            * median_primary_prominence
        ):
            continue
        recovered.append(candidate)

    return np.asarray(recovered, dtype=np.int32)


def _prominence_window(nominal_period: float, signal_size: int) -> int | None:
    window = min(signal_size, max(3, int(round(2.0 * nominal_period))))
    if window % 2 == 0:
        window -= 1
    return window if window >= 3 else None


def _interval_period_ratio(peaks: np.ndarray, nominal_period: float) -> np.ndarray:
    if peaks.size < 2:
        return np.empty(0, dtype=np.float32)
    if not np.isfinite(nominal_period) or nominal_period <= 0:
        return np.full(peaks.size - 1, np.nan, dtype=np.float32)
    return (
        np.diff(peaks).astype(np.float32, copy=False) / np.float32(nominal_period)
    ).astype(np.float32, copy=False)


def _valid_intervals(period_ratios: np.ndarray) -> np.ndarray:
    ratios = np.asarray(period_ratios, dtype=np.float32)
    return (np.isfinite(ratios) & (ratios >= 0.55) & (ratios <= 1.60)).astype(
        np.bool_,
        copy=False,
    )


def _scipy_signal_dependencies():
    try:
        from scipy.signal import find_peaks, peak_prominences
    except ModuleNotFoundError as exc:
        raise ImportError("Systole detection requires scipy.") from exc
    return find_peaks, peak_prominences
