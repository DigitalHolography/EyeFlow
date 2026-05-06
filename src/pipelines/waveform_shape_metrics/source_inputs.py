"""Read fixed HD/DV inputs for the waveform-shape metrics pipeline."""

from __future__ import annotations

import numpy as np

from input_output import DOPPLER_VIEW_SCHEMA, HOLODOPPLER_SCHEMA
from input_output.input_access import (
    read_nested_int_setting,
    resolve_required_source_array,
)

from .constants import DOPPLERVIEW_DEFAULT_LOCAL_BACKGROUND_DIST


def dopplerview_cache_from_h5(ctx) -> dict[str, object]:
    return {
        "moment0": _read_required_float_array(
            ctx.hd,
            "HD",
            "moment0",
            HOLODOPPLER_SCHEMA.dataset_path("moment0"),
            dopplerview_moment=True,
        ),
        "moment2": _read_required_float_array(
            ctx.hd,
            "HD",
            "moment2",
            HOLODOPPLER_SCHEMA.dataset_path("moment2"),
            dopplerview_moment=True,
        ),
        "retinal_artery_mask": _read_required_bool_array(
            ctx.dv,
            "DV",
            "retinal artery mask",
            DOPPLER_VIEW_SCHEMA.dataset_path("retinal_artery_mask"),
        ),
        "retinal_vein_mask": _read_required_bool_array(
            ctx.dv,
            "DV",
            "retinal vein mask",
            DOPPLER_VIEW_SCHEMA.dataset_path("retinal_vein_mask"),
        ),
        "retinal_labeled_vessels": _read_optional_int_array(
            ctx.dv,
            DOPPLER_VIEW_SCHEMA.dataset_path("retinal_labeled_vessels"),
        ),
    }


def local_background_dist(ctx) -> int:
    spec = DOPPLER_VIEW_SCHEMA.config_value("local_background_dist")
    value = read_nested_int_setting(
        ctx.dv_config,
        spec.section or "",
        spec.json_key,
        default=int(spec.default or DOPPLERVIEW_DEFAULT_LOCAL_BACKGROUND_DIST),
    )
    return int(value)


def _read_required_float_array(
    source,
    source_name: str,
    logical_name: str,
    path: str,
    *,
    dopplerview_moment: bool = False,
):
    resolved = resolve_required_source_array(
        source,
        source_name=source_name,
        logical_name=logical_name,
        path=path,
    )
    value = (
        _coerce_dopplerview_moment(resolved.value)
        if dopplerview_moment
        else resolved.value
    )
    return np.asarray(value, dtype=np.float32)


def _coerce_dopplerview_moment(value) -> np.ndarray:
    squeezed = np.squeeze(np.asarray(value, dtype=np.float32))
    if squeezed.ndim != 3:
        raise ValueError(
            "Holodoppler moment datasets must become 3-D after squeeze, "
            f"got shape {squeezed.shape}."
        )
    return np.transpose(squeezed, (0, 2, 1))


def _read_required_bool_array(source, source_name: str, logical_name: str, path: str):
    resolved = resolve_required_source_array(
        source,
        source_name=source_name,
        logical_name=logical_name,
        path=path,
    )
    return np.asarray(resolved.value, dtype=bool)


def _read_optional_int_array(source, path: str) -> np.ndarray | None:
    if source.h5file is None:
        return None
    found = source.h5file.get(path)
    if found is None:
        return None
    return np.asarray(found[()], dtype=np.int32)
