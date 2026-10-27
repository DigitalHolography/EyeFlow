"""Shared pipeline-layer packing for per-beat profile datasets."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from calculations.math import nanmean_float32
from calculations.vessel_segments.profiles.per_beat import interpolate_profiles_per_beat
from input_output.writers.h5 import profile_h5_options
from pipeline_engine.base import DatasetValue


def profile_dataset(
    profiles: np.ndarray,
    cycle_boundary_indexes,
    *,
    index_base: int,
    spatial_axis: str = "x",
    unit: str = "mm/s",
    attrs: Mapping[str, object] | None = None,
    valid_segments: np.ndarray | None = None,
) -> DatasetValue:
    """Interpolate profiles per beat and attach their persisted data contract."""

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
    output_attrs = dict(attrs or {})
    output_attrs.update(
        {
            "unit": unit,
            "dimDesc": [spatial_axis, "time", "beat", "branch", "radius"],
        }
    )
    return DatasetValue(
        data=profiles_per_beat,
        attrs=output_attrs,
        h5_options={
            **profile_h5_options(profiles_per_beat.shape),
            "sparse_nan_chunks": True,
        },
    )


def temporally_meaned_profile_dataset(profile: DatasetValue) -> DatasetValue:
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
        h5_options={**profile_h5_options(data.shape), "sparse_nan_chunks": True},
    )


__all__ = [
    "profile_dataset",
    "profile_h5_options",
    "temporally_meaned_profile_dataset",
]
