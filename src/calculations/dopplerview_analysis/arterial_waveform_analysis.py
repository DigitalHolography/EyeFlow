"""Arterial waveform analysis, using BloodFlowVelocity/find_systole_index.m."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.analysis_preparation.beat_detection import find_systole_index
from calculations.math import butter_lowpass_filtfilt


def analyze_arterial_waveforms(
    artery_signal,
    vein_signal,
    *,
    dt_seconds: np.float32,
    lowpass_freq_hz: np.float32 = np.float32(15.0),
) -> dict[str, object]:
    """Return refreshed waveform analysis and auditable beat diagnostics."""

    artery = np.asarray(artery_signal, dtype=np.float32).reshape(-1)
    vein = np.asarray(vein_signal, dtype=np.float32).reshape(-1)
    detection = find_systole_index(
        artery,
        dt=np.float32(dt_seconds),
        lowpass_freq_hz=np.float32(lowpass_freq_hz),
    )
    peaks = detection.systole_indexes
    vein_filtered = butter_lowpass_filtfilt(
        vein,
        dt_seconds=np.float32(dt_seconds),
        lowpass_freq_hz=np.float32(lowpass_freq_hz),
        order=4,
    )
    return {
        "retinal_artery_velocity_signal_filtered_perbeat": slice_interp_beats(
            peaks,
            detection.artery_signal_filtered,
        ),
        "retinal_artery_velocity_signal_filtered": detection.artery_signal_filtered,
        "retinal_artery_velocity_signal_derivative": detection.derivative_signal,
        "retinal_vein_velocity_signal_filtered": vein_filtered,
        "retinal_vein_velocity_signal_derivative": np.gradient(vein_filtered).astype(
            np.float32
        ),
        "beat_indices": peaks,
        "time_per_beat": (
            np.diff(peaks).astype(np.float32) * np.float32(dt_seconds)
        ).astype(np.float32),
        "beat_detection_min_peak_distance": detection.min_peak_distance,
        "beat_detection_min_peak_height": detection.min_peak_height,
        "beat_detection_initial_indices": detection.initial_systole_indexes,
        "beat_detection_recovered_indices": detection.recovered_systole_indexes,
        "beat_detection_nominal_period_samples": detection.nominal_period_samples,
        "beat_detection_interval_period_ratio": detection.interval_period_ratio,
        "beat_detection_interval_duration_valid": detection.interval_duration_valid,
    }


def slice_interp_beats(peaks, signal, ninterp: int = 128) -> np.ndarray:
    peaks = np.asarray(peaks, dtype=np.int32).reshape(-1)
    signal = np.asarray(signal, dtype=np.float32).reshape(-1)
    nbeat = max(0, peaks.size - 1)
    signal_per_beat = np.zeros(shape=(nbeat, ninterp), dtype=np.float32)

    for beat_index in range(nbeat):
        beat_signal = signal[peaks[beat_index] : peaks[beat_index + 1]]
        beat_interp = np.interp(
            np.linspace(0, 1, ninterp, dtype=np.float32),
            np.linspace(0, 1, len(beat_signal), dtype=np.float32),
            beat_signal,
        ).astype(np.float32)
        signal_per_beat[beat_index, :] = beat_interp
    return signal_per_beat


class ArterialWaveformAnalysisStep:
    def _relevant_config(self, ctx):
        return {
            "sampling_freq": ctx.hd_config_value("sampling_freq"),
            "stride": ctx.hd_config_value("batch_stride"),
            "LowpassFreqHz": ctx.dv_config_value(
                "PulseAnalysis",
                "LowpassFreqHz",
                15.0,
            ),
        }

    def slice_interp_beats(self, peaks, sig, ninterp=128):
        return slice_interp_beats(peaks, sig, ninterp)

    def run(self, ctx):
        # ---- Requires ----
        sig = ctx.require("retinal_artery_velocity_signal")
        stride = np.float32(ctx.hd_config_value("batch_stride"))
        fs = np.float32(ctx.hd_config_value("sampling_freq"))
        dt = stride / fs

        outputs = analyze_arterial_waveforms(
            sig,
            ctx.require("retinal_vein_velocity_signal"),
            dt_seconds=dt,
            lowpass_freq_hz=np.float32(
                ctx.dv_config_value("PulseAnalysis", "LowpassFreqHz", 15.0)
            ),
        )
        for key, value in outputs.items():
            ctx.set(key, value)
