"""Build per-beat profile datasets shared by several pipelines."""

import numpy as np

from calculations.math import nanmean_float32
from calculations.topology.profile_interpolation import interpolate_profiles_per_beat
from input_output.writers.h5 import profile_h5_options
from pipeline_engine.base import DatasetValue


def _profile_dataset(
    profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
    spatial_axis: str = "x",
    unit: str = "mm/s",
    valid_segments: np.ndarray | None = None,
) -> DatasetValue:
    if profiles.ndim != 4:
        raise ValueError(
            "profile arrays must have shape "
            "(radius, branch, frame, spatial_sample)."
        )
    profiles_per_beat = interpolate_profiles_per_beat(
        profiles,
        cycle_boundary_indexes,
        index_base=index_base,
        valid_segments=valid_segments,
    )
    return DatasetValue(
        data=profiles_per_beat,
        attrs={
            "unit": unit,
            "dimDesc": [spatial_axis, "time", "beat", "branch", "radius"],
        },
        h5_options=profile_h5_options(profiles_per_beat.shape),
    )


def _temporally_meaned_profile_dataset(profile: DatasetValue) -> DatasetValue:
    """Average an interpolated profile over time within each beat."""
    data = nanmean_float32(np.asarray(profile.data), axis=1)
    attrs = dict(profile.attrs or {})
    dim_desc = list(attrs.get("dimDesc", ()))
    if len(dim_desc) < 2 or dim_desc[1] != "time":
        raise ValueError("profile dataset must have time as its second dimension.")
    del dim_desc[1]
    attrs["dimDesc"] = dim_desc
    attrs["temporal_reduction"] = "mean_over_interpolated_beat_time"
    return DatasetValue(
        data=data,
        attrs=attrs,
        h5_options=profile_h5_options(data.shape),
    )
