"""Run waveform-velocity calculations for each cardiac cycle."""

from calculations.blood_flow_velocity import run_per_beat_analysis


def analyze_velocity_per_beat(inputs):
    return run_per_beat_analysis(inputs)


__all__ = ["analyze_velocity_per_beat"]
