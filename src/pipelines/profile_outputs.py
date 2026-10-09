"""Shared pipeline-layer packing for per-beat profile datasets."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from calculations.math import nanmean_float32
from calculations.topology.profile_interpolation import interpolate_profiles_per_beat
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
        h5_options=profile_h5_options(profiles_per_beat.shape),
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
        h5_options=profile_h5_options(data.shape),
    )


def profile_h5_options(shape: tuple[int, ...]) -> dict[str, object]:
    """Use lossless compression with chunks aligned to one segment profile."""

    options: dict[str, object] = {
        "compression": "gzip",
        "compression_opts": 4,
        "shuffle": True,
    }
    if len(shape) not in (4, 5) or not all(shape):
        return options

    sample_count, time_count = shape[:2]
    target_elements = (1024 * 1024) // np.dtype(np.float32).itemsize
    if len(shape) == 4:
        options["chunks"] = (sample_count, 1, 1, 1)
        return options
    time_chunk = min(time_count, max(target_elements // sample_count, 1))
    options["chunks"] = (sample_count, time_chunk, 1, 1, 1)
    return options


__all__ = [
    "profile_dataset",
    "profile_h5_options",
    "temporally_meaned_profile_dataset",
]
