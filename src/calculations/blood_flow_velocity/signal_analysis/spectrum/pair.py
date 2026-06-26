"""Two-signal spectrum, coherence, transfer, and delay calculations."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, signal

from calculations.blood_flow_velocity.signal_analysis.waveform import (
    average_cycle,
    nan_to_mean,
    rescale,
    standardize,
)

from .signal import spectrum_signal_analysis


@dataclass(frozen=True)
class CorrelationData:
    first: np.ndarray
    second: np.ndarray
    lags_seconds: np.ndarray
    cross_corr: np.ndarray
    time_lag: float
    max_corr: float
    coherence_freq: np.ndarray
    coherence: np.ndarray
    heart_rate_hz: float
    gamma_0: float


@dataclass(frozen=True)
class TransferData:
    frequencies: np.ndarray
    transfer: np.ndarray


@dataclass(frozen=True)
class DelayFitData:
    time: np.ndarray
    first: np.ndarray
    second_scaled: np.ndarray
    model: np.ndarray
    tau_samples: float
    tau_milliseconds: float
    shift_milliseconds: float


@dataclass(frozen=True)
class PairedSpectrumAnalysisResult:
    correlation: CorrelationData
    transfer: TransferData
    delay: DelayFitData | None


def paired_spectrum_analysis(
    first: np.ndarray,
    second: np.ndarray,
    dt_seconds: float,
    beat_indexes: np.ndarray | None = None,
) -> PairedSpectrumAnalysisResult:
    delay = None
    if beat_indexes is not None:
        peaks = np.asarray(beat_indexes, dtype=np.int32).reshape(-1)
        first_cycle = average_cycle(first, peaks, 128)
        second_cycle = average_cycle(second, peaks, 128)
        if first_cycle is not None and second_cycle is not None:
            delay = delay_fit_analysis(first_cycle, second_cycle, peaks, dt_seconds)
    return PairedSpectrumAnalysisResult(
        correlation=correlation_data(first, second, dt_seconds),
        transfer=transfer_function(first, second, dt_seconds),
        delay=delay,
    )


def correlation_data(first: np.ndarray, second: np.ndarray, dt_seconds: float) -> CorrelationData:
    a = standardize(first)
    v = standardize(second)
    cross = signal.correlate(a, v, mode="full")
    denom = max(float(a.size * np.nanstd(a) * np.nanstd(v)), np.finfo(np.float32).eps)
    cross = cross / denom
    lags = signal.correlation_lags(a.size, v.size, mode="full")
    max_idx = int(np.nanargmax(cross))
    fs = 1.0 / dt_seconds
    nperseg = min(64, a.size)
    freqs, coherence = signal.coherence(
        a,
        v,
        fs=fs,
        window=np.hamming(nperseg),
        nperseg=nperseg,
        noverlap=None,
        nfft=max(256, nperseg),
    )
    signal_fft = spectrum_signal_analysis(a, dt_seconds)
    heart_rate = (
        float(signal_fft.frequencies[signal_fft.peak_indexes[0]])
        if signal_fft.peak_indexes.size
        else np.nan
    )
    coherence_at_heart_rate = gamma_0(coherence, freqs, heart_rate)
    return CorrelationData(
        a.astype(np.float32),
        v.astype(np.float32),
        (lags * dt_seconds).astype(np.float32),
        cross.astype(np.float32),
        float(lags[max_idx] * dt_seconds),
        float(cross[max_idx]),
        freqs.astype(np.float32),
        coherence.astype(np.float32),
        heart_rate,
        coherence_at_heart_rate,
    )


def transfer_function(
    input_signal: np.ndarray,
    output_signal: np.ndarray,
    dt_seconds: float,
) -> TransferData:
    n = max(np.asarray(input_signal).size, 1) * 10
    input_fft = np.fft.fft(nan_to_mean(input_signal), n)
    output_fft = np.fft.fft(nan_to_mean(output_signal), n)
    transfer = output_fft / np.where(np.abs(input_fft) == 0, np.nan + 0j, input_fft)
    freqs = np.fft.fftfreq(n, d=dt_seconds)
    order = np.argsort(freqs)
    return TransferData(freqs[order].astype(np.float32), transfer[order].astype(np.complex64))


def delay_fit_analysis(
    first_cycle: np.ndarray,
    second_cycle: np.ndarray,
    beat_indexes: np.ndarray,
    dt_seconds: float,
) -> DelayFitData:
    first = rescale(first_cycle)
    second = rescale(second_cycle)
    amin = int(np.nanargmin(second_cycle))
    samples = np.arange(first.size, dtype=np.float32)

    def model_tau(tau):
        kernel = np.exp(-samples / max(tau, 1e-6)) / max(tau, 1e-6)
        result = np.convolve(first, kernel, mode="full")[: first.size]
        return rescale(np.roll(result, amin))

    target = np.zeros(second.size, dtype=np.float32)
    target[amin:] = 1.0

    def objective(tau):
        return float(np.sum((rescale(second) - model_tau(float(tau[0]))) ** 2 * target))

    fitted = optimize.minimize(objective, x0=np.asarray([40.0]), method="Nelder-Mead")
    tau = float(fitted.x[0]) if fitted.success else 40.0
    if beat_indexes.size > 1:
        period = float(np.nanmean(np.diff(beat_indexes)) * dt_seconds)
    else:
        period = first.size * dt_seconds
    model = model_tau(tau)
    return DelayFitData(
        time=np.linspace(0, period, first.size, dtype=np.float32),
        first=first.astype(np.float32),
        second_scaled=(second * float(np.nanmax(model))).astype(np.float32),
        model=model.astype(np.float32),
        tau_samples=tau,
        tau_milliseconds=float(tau / first.size * period * 1000),
        shift_milliseconds=float(amin / first.size * period * 1000),
    )


def gamma_0(coherence: np.ndarray, frequencies: np.ndarray, heart_rate_hz: float) -> float:
    if not np.isfinite(heart_rate_hz):
        return np.nan
    valid = (frequencies < heart_rate_hz + 0.3) & (frequencies > heart_rate_hz - 0.3)
    if not np.any(valid):
        return np.nan
    return float(np.nanmean(coherence[valid]))
