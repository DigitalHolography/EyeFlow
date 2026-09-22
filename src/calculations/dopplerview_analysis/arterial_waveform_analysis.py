"""Arterial waveform analysis, using BloodFlowVelocity/find_systole_index.m."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.heartbeat import (
    HeartbeatAnalysisResult,
    SystoleDetectionResult,
    run_heartbeat_analysis,
    spectral_heartbeat_analysis,
)
from calculations.math import butter_lowpass_filtfilt


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
        nbeat = max(0, len(peaks) - 1)

        sig_perbeat = np.zeros(shape=(nbeat, ninterp), dtype=np.float32)

        for i in range(nbeat):
            beat_sig = sig[peaks[i]:peaks[i+1]]
            beat_sig_interp = np.interp(
                np.linspace(0, 1, ninterp, dtype=np.float32),
                np.linspace(0, 1, len(beat_sig), dtype=np.float32),
                beat_sig,
            ).astype(np.float32)
            sig_perbeat[i,:] = beat_sig_interp

        return sig_perbeat

    def run(self, ctx):
        # ---- Requires ----
        sig = np.asarray(
            ctx.require("retinal_artery_velocity_signal"),
            dtype=np.float32,
        )
        vein_sig = np.asarray(
            ctx.require("retinal_vein_velocity_signal"),
            dtype=np.float32,
        )
        stride = np.float32(ctx.hd_config_value("batch_stride"))
        fs = np.float32(ctx.hd_config_value("sampling_freq"))
        dt = stride / fs
        lowpass_freq_hz = np.float32(
            ctx.dv_config_value("PulseAnalysis", "LowpassFreqHz", 15.0)
        )
        sig_filtered = butter_lowpass_filtfilt(
            sig,
            dt_seconds=dt,
            lowpass_freq_hz=lowpass_freq_hz,
            order=4,
        )
        vein_filtered = butter_lowpass_filtfilt(
            vein_sig,
            dt_seconds=dt,
            lowpass_freq_hz=lowpass_freq_hz,
            order=4,
        )
        heartbeat, beat_detection_source = _heartbeat_from_available_vessel(
            sig,
            vein_sig,
            dt_seconds=float(dt),
            lowpass_freq_hz=float(lowpass_freq_hz),
        )
        detection = heartbeat.systole
        peaks = detection.systole_indexes
        artery_derivative = np.gradient(sig_filtered, dt).astype(np.float32)

        sig_perbeat = self.slice_interp_beats(peaks, sig_filtered)

        ctx.set("retinal_artery_velocity_signal_filtered_perbeat", sig_perbeat)
        ctx.set("retinal_artery_velocity_signal_filtered", sig_filtered)
        ctx.set("retinal_artery_velocity_signal_derivative", artery_derivative)
        ctx.set("retinal_vein_velocity_signal_filtered", vein_filtered)
        ctx.set(
            "retinal_vein_velocity_signal_derivative",
            np.gradient(vein_filtered, dt).astype(np.float32),
        )
        ctx.set("beat_indices", peaks)
        ctx.set(
            "time_per_beat",
            (np.diff(peaks).astype(np.float32) * dt).astype(np.float32),
        ) # TODO parametrize look for params
        ctx.set("beat_detection_min_peak_distance", detection.min_peak_distance)
        ctx.set("beat_detection_min_peak_height", detection.min_peak_height)
        ctx.set("beat_detection_source", beat_detection_source)
        ctx.set("_heartbeat_analysis_result", heartbeat)


def _heartbeat_from_available_vessel(
    artery_signal: np.ndarray,
    vein_signal: np.ndarray,
    *,
    dt_seconds: float,
    lowpass_freq_hz: float,
) -> tuple[HeartbeatAnalysisResult, str]:
    """Detect shared beat boundaries without requiring an arterial mask."""

    for source_name, signal in (
        ("artery", artery_signal),
        ("vein", vein_signal),
    ):
        if np.any(np.isfinite(signal)):
            try:
                heartbeat = run_heartbeat_analysis(
                    signal,
                    dt_seconds=dt_seconds,
                    lowpass_freq_hz=lowpass_freq_hz,
                )
            except ValueError as exc:
                if "No systole peaks detected" not in str(exc):
                    raise
                continue
            return heartbeat, source_name
    return _missing_vessel_heartbeat(artery_signal.size, dt_seconds), "none"


def _missing_vessel_heartbeat(
    signal_length: int,
    dt_seconds: float,
) -> HeartbeatAnalysisResult:
    """Return one full-record cycle when neither vessel can provide timing."""

    if signal_length < 2:
        raise ValueError("At least two frames are required for per-beat analysis.")
    boundaries = np.asarray([0, signal_length - 1], dtype=np.int32)
    missing = np.full(signal_length, np.nan, dtype=np.float32)
    spectral = spectral_heartbeat_analysis(
        missing,
        dt_seconds,
        systole_count=0,
    )
    return HeartbeatAnalysisResult(
        systole=SystoleDetectionResult(
            systole_indexes=boundaries,
            artery_signal_filtered=missing,
            derivative_signal=missing.copy(),
            min_peak_distance=max(1, int(np.floor(0.5 / dt_seconds))),
            min_peak_height=np.float32(np.nan),
        ),
        spectral=spectral,
    )
