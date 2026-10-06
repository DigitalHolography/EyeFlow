"""Adaptive systolic-boundary detection for retinal velocity signals."""

from __future__ import annotations

from itertools import pairwise

import numpy as np

from calculations.math import butter_lowpass_filtfilt

from .models import SuspectedMissedBeatGap, SystoleDetectionResult

DEFAULT_MIN_PERIOD_SECONDS = np.float32(0.2)
DEFAULT_MAX_PERIOD_SECONDS = np.float32(2.0)
DEFAULT_LOWPASS_FREQUENCY_HZ = np.float32(15.0)
DEFAULT_ADAPTIVE_DISTANCE_RATIO = np.float32(0.5)
DEFAULT_PROMINENCE_MAD_MULTIPLIER = np.float32(0.75)
DEFAULT_PERIOD_AGREEMENT_RATIO = np.float32(0.25)
DEFAULT_PERIOD_MARGIN_RATIO = np.float32(0.1)
DEFAULT_PERIOD_VARIABILITY_SIGMAS = np.float32(3.0)
ROBUST_SIGMA_FROM_MAD = np.float32(1.4826)
MISSED_BEAT_MULTIPLES = (2, 3)
MISSED_BEAT_MULTIPLE_TOLERANCE = np.float32(0.25)


class SystoleDetectionError(ValueError):
    """Raised when a signal cannot provide two usable systolic boundaries."""


def find_systole_index(
    signal,
    *,
    dt: np.float32,
    lowpass_freq_hz: np.float32 = DEFAULT_LOWPASS_FREQUENCY_HZ,
    min_duration_seconds: np.float32 = DEFAULT_MIN_PERIOD_SECONDS,
    max_duration_seconds: np.float32 = DEFAULT_MAX_PERIOD_SECONDS,
    adaptive_distance_ratio: np.float32 = DEFAULT_ADAPTIVE_DISTANCE_RATIO,
    prominence_mad_multiplier: np.float32 = DEFAULT_PROMINENCE_MAD_MULTIPLIER,
) -> SystoleDetectionResult:
    """Detect systolic upstrokes without assuming a fixed heart-rate ceiling.

    A permissive physiological spacing produces provisional peaks. Their median
    interval is cross-checked against the signal autocorrelation. The final
    refractory distance is the estimated period minus a robust allowance for
    observed period variability. Missing boundaries are retained as long cycles
    and reported on the result rather than synthesized or discarded.
    """

    find_peaks, correlate = _scipy_signal_dependencies()

    pulse = np.asarray(signal, dtype=np.float32).reshape(-1)
    _validate_detection_inputs(
        pulse,
        dt,
        min_duration_seconds,
        max_duration_seconds,
        adaptive_distance_ratio,
        prominence_mad_multiplier,
    )
    filtered_pulse = butter_lowpass_filtfilt(
        pulse,
        dt_seconds=np.float32(dt),
        lowpass_freq_hz=np.float32(lowpass_freq_hz),
        order=4,
    )
    derivative = np.gradient(filtered_pulse, np.float32(dt)).astype(np.float32)
    min_peak_prominence = _robust_peak_prominence(
        derivative,
        prominence_mad_multiplier,
    )
    candidate_distance = _min_peak_distance(dt, min_duration_seconds)

    provisional_peaks, _ = find_peaks(
        derivative,
        prominence=float(min_peak_prominence),
        distance=candidate_distance,
    )
    median_period = _median_peak_period(provisional_peaks)
    autocorrelation_period = _autocorrelation_period(
        filtered_pulse,
        candidate_distance,
        _max_peak_distance(dt, max_duration_seconds, pulse.size),
        find_peaks,
        correlate,
    )
    estimated_period = _resolve_period_estimate(
        median_period,
        autocorrelation_period,
        agreement_ratio=DEFAULT_PERIOD_AGREEMENT_RATIO,
    )
    min_peak_distance = _adaptive_peak_distance(
        candidate_distance,
        estimated_period,
        adaptive_distance_ratio,
        provisional_peaks,
    )

    peaks, _ = find_peaks(
        derivative,
        prominence=float(min_peak_prominence),
        distance=min_peak_distance,
    )
    indexes = peaks.astype(np.int32, copy=False)
    if indexes.size < 2:
        raise SystoleDetectionError(
            "Fewer than two systole peaks were detected. Check signal quality or parameters."
        )
    final_median_period = _median_peak_period(indexes)
    estimated_period = _resolve_period_estimate(
        final_median_period,
        autocorrelation_period,
        agreement_ratio=DEFAULT_PERIOD_AGREEMENT_RATIO,
    )
    missed_beat_gaps = _suspected_missed_beat_gaps(indexes, estimated_period)
    return SystoleDetectionResult(
        systole_indexes=indexes,
        signal_filtered=filtered_pulse,
        derivative_signal=derivative,
        min_peak_distance=min_peak_distance,
        min_peak_height=np.float32(np.nan),
        min_peak_prominence=min_peak_prominence,
        estimated_period_samples=np.float32(estimated_period),
        suspected_missed_beat_gaps=missed_beat_gaps,
    )


def _min_peak_distance(dt_seconds: np.float32, min_duration_seconds: np.float32) -> int:
    if dt_seconds <= 0:
        raise ValueError("dt_seconds must be positive for systole detection.")
    return max(1, int(np.floor(float(min_duration_seconds) / float(dt_seconds))))


def _max_peak_distance(
    dt_seconds: np.float32,
    max_duration_seconds: np.float32,
    signal_length: int,
) -> int:
    return min(
        max(1, int(np.ceil(float(max_duration_seconds) / float(dt_seconds)))),
        max(1, int(signal_length) // 2),
    )


def _validate_detection_inputs(
    pulse: np.ndarray,
    dt_seconds: np.float32,
    min_duration_seconds: np.float32,
    max_duration_seconds: np.float32,
    adaptive_distance_ratio: np.float32,
    prominence_mad_multiplier: np.float32,
) -> None:
    values = np.asarray(
        (
            dt_seconds,
            min_duration_seconds,
            max_duration_seconds,
            adaptive_distance_ratio,
            prominence_mad_multiplier,
        ),
        dtype=np.float64,
    )
    if pulse.size < 2:
        raise ValueError("At least two signal samples are required for systole detection.")
    if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
        raise ValueError("Systole detection timing and threshold parameters must be positive.")
    if min_duration_seconds >= max_duration_seconds:
        raise ValueError("min_duration_seconds must be less than max_duration_seconds.")
    if adaptive_distance_ratio > 1.0:
        raise ValueError("adaptive_distance_ratio must not exceed one.")


def _robust_peak_prominence(
    derivative: np.ndarray,
    multiplier: np.float32,
) -> np.float32:
    finite = derivative[np.isfinite(derivative)]
    if finite.size == 0:
        return np.float32(np.inf)
    center = np.median(finite)
    mad = np.median(np.abs(finite - center))
    robust_sigma = ROBUST_SIGMA_FROM_MAD * np.float32(mad)
    scale = max(float(np.max(np.abs(finite))), 1.0)
    numerical_floor = np.finfo(np.float32).eps * scale
    return np.float32(max(float(multiplier * robust_sigma), numerical_floor))


def _median_peak_period(peaks: np.ndarray) -> float:
    indexes = np.asarray(peaks, dtype=np.int32).reshape(-1)
    if indexes.size < 2:
        return np.nan
    return float(np.median(np.diff(indexes)))


def _autocorrelation_period(
    signal: np.ndarray,
    min_lag: int,
    max_lag: int,
    find_peaks,
    correlate,
) -> float:
    values = np.asarray(signal, dtype=np.float64).reshape(-1)
    if values.size < 3 or not np.all(np.isfinite(values)):
        return np.nan
    centered = values - np.mean(values)
    total_energy = float(np.dot(centered, centered))
    if total_energy <= np.finfo(np.float64).eps:
        return np.nan

    correlation = correlate(centered, centered, mode="full", method="fft")[
        values.size - 1 :
    ]
    squared = centered * centered
    cumulative = np.concatenate(([0.0], np.cumsum(squared)))
    lags = np.arange(correlation.size, dtype=np.int32)
    left_energy = cumulative[values.size - lags]
    right_energy = cumulative[values.size] - cumulative[lags]
    denominator = np.sqrt(left_energy * right_energy)
    normalized = np.divide(
        correlation,
        denominator,
        out=np.full_like(correlation, np.nan),
        where=denominator > np.finfo(np.float64).eps,
    )

    start = max(1, int(min_lag))
    stop = min(int(max_lag) + 1, normalized.size)
    if stop - start < 3:
        return np.nan
    window = normalized[start:stop]
    local_peaks, _ = find_peaks(window, prominence=0.05)
    if local_peaks.size == 0:
        return np.nan
    strengths = window[local_peaks]
    maximum = float(np.nanmax(strengths))
    if not np.isfinite(maximum) or maximum <= 0.1:
        return np.nan
    strong = local_peaks[strengths >= 0.5 * maximum]
    if strong.size == 0:
        return np.nan
    return float(start + int(strong[0]))


def _resolve_period_estimate(
    median_period: float,
    autocorrelation_period: float,
    *,
    agreement_ratio: np.float32 = DEFAULT_PERIOD_AGREEMENT_RATIO,
) -> float:
    median_valid = np.isfinite(median_period) and median_period > 0
    autocorrelation_valid = np.isfinite(autocorrelation_period) and autocorrelation_period > 0
    if median_valid and not autocorrelation_valid:
        return float(median_period)
    if autocorrelation_valid and not median_valid:
        return float(autocorrelation_period)
    if not median_valid and not autocorrelation_valid:
        return np.nan
    relative_difference = abs(median_period - autocorrelation_period) / float(
        autocorrelation_period
    )
    if relative_difference <= float(agreement_ratio):
        return float(median_period)
    return float(autocorrelation_period)


def _adaptive_peak_distance(
    minimum_distance: int,
    estimated_period: float,
    minimum_period_ratio: np.float32,
    provisional_peaks: np.ndarray,
) -> int:
    if not np.isfinite(estimated_period) or estimated_period <= 0:
        return int(minimum_distance)
    period = float(estimated_period)
    ratio_floor = int(np.floor(period * float(minimum_period_ratio)))
    variability = _robust_interval_variability(provisional_peaks)
    if not np.isfinite(variability):
        return max(int(minimum_distance), ratio_floor, 1)
    variability_margin = max(
        period * float(DEFAULT_PERIOD_MARGIN_RATIO),
        float(DEFAULT_PERIOD_VARIABILITY_SIGMAS) * variability,
    )
    statistical_distance = int(np.floor(period - variability_margin))
    return max(int(minimum_distance), ratio_floor, statistical_distance, 1)


def _robust_interval_variability(peaks: np.ndarray) -> float:
    indexes = np.asarray(peaks, dtype=np.int32).reshape(-1)
    if indexes.size < 3:
        return np.nan
    intervals = np.diff(indexes).astype(np.float64)
    center = np.median(intervals)
    mad = np.median(np.abs(intervals - center))
    return float(ROBUST_SIGMA_FROM_MAD) * float(mad)


def _suspected_missed_beat_gaps(
    indexes: np.ndarray,
    estimated_period: float,
) -> tuple[SuspectedMissedBeatGap, ...]:
    if not np.isfinite(estimated_period) or estimated_period <= 0:
        return ()
    gaps = []
    for start, stop in pairwise(indexes):
        interval = int(stop) - int(start)
        ratio = interval / float(estimated_period)
        multiple = int(np.rint(ratio))
        if multiple not in MISSED_BEAT_MULTIPLES:
            continue
        if abs(ratio - multiple) > float(MISSED_BEAT_MULTIPLE_TOLERANCE):
            continue
        gaps.append(
            SuspectedMissedBeatGap(
                start_index=int(start),
                stop_index=int(stop),
                interval_samples=interval,
                estimated_period_samples=float(estimated_period),
                estimated_multiple=multiple,
            )
        )
    return tuple(gaps)


def _scipy_signal_dependencies():
    try:
        from scipy.signal import correlate, find_peaks
    except ModuleNotFoundError as exc:
        raise ImportError("Systole detection requires scipy.") from exc
    return find_peaks, correlate
