"""Arterial waveform analysis, using BloodFlowVelocity/find_systole_index.m."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.find_systole_index import find_systole_index

from .base import CalculationStep as BaseStep


class ArterialWaveformAnalysisStep(BaseStep):
    name = "arterial_waveform_analysis"
    requires = {"retinal_artery_velocity_signal"}
    produces = {
        "retinal_artery_velocity_signal_filtered",
        "retinal_artery_velocity_signal_filtered_perbeat",
        "beat_indices",
        "time_per_beat",
    }

    def _relevant_config(self, ctx):
        return {
            "sampling_freq": ctx.holodoppler_config["sampling_freq"],
            "stride": ctx.holodoppler_config["batch_stride"],
        }

    def slice_interp_beats(self, peaks, sig):
        nbeat = max(0, len(peaks) - 1)

        ninterp = 128 # TODO parametrize

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
        sig = ctx.require("retinal_artery_velocity_signal")
        stride = np.float32(ctx.holodoppler_config["batch_stride"])
        fs = np.float32(ctx.holodoppler_config["sampling_freq"])
        dt = stride / fs

        detection = find_systole_index(sig, dt=dt)
        peaks = detection.systole_indexes
        sig_filtered = detection.artery_signal_filtered

        sig_perbeat = self.slice_interp_beats(peaks, sig_filtered)

        ctx.set("retinal_artery_velocity_signal_filtered_perbeat", sig_perbeat)
        ctx.set("retinal_artery_velocity_signal_filtered", sig_filtered)
        ctx.set("beat_indices", peaks)
        ctx.set(
            "time_per_beat",
            (np.diff(peaks).astype(np.float32) * dt).astype(np.float32),
        ) # TODO parametrize look for params
        ctx.set("beat_detection_min_peak_distance", detection.min_peak_distance)
        ctx.set("beat_detection_min_peak_height", detection.min_peak_height)
