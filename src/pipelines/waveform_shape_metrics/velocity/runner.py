"""Run and pack per-beat velocity calculations."""

from calculations.blood_flow_velocity import run_per_beat_analysis

from .outputs import pack_velocity_per_beat_outputs


def run_velocity_per_beat_metrics(context):
    result = run_per_beat_analysis(context.per_beat_analysis)
    return result, pack_velocity_per_beat_outputs(result)
