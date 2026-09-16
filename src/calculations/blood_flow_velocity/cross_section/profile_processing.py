"""Compatibility adapter for generic per-beat topology profiles."""
from calculations.topology.profile_interpolation import interpolate_profiles_per_beat


def interpolate_velocity_profiles_per_beat(
    velocity_profiles, cycle_boundary_indexes, *, index_base=0, valid_segments=None,
):
    return interpolate_profiles_per_beat(
        velocity_profiles, cycle_boundary_indexes, index_base=index_base,
        valid_segments=valid_segments,
    )


__all__ = ["interpolate_velocity_profiles_per_beat"]
