"""Resolve the minimal inputs required to detect heartbeat boundaries."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from input_output.schema import HolodopplerTiming


@dataclass(frozen=True)
class HeartbeatInputs:
    moment0: object
    moment2: object
    artery_mask: np.ndarray
    vein_mask: np.ndarray
    optic_disc_center: np.ndarray | None
    optic_disc_width: np.ndarray | None
    optic_disc_height: np.ndarray | None
    timing: HolodopplerTiming
    local_background_dist: int
    index_base: int = 0


def load_heartbeat_inputs(ctx) -> HeartbeatInputs:
    """Load and spatially align only the arrays used for beat detection."""

    ctx.require_inputs("hd", "dv")
    hd = ctx.inputs.hd.as_holodoppler()
    dv = ctx.inputs.dv.as_dopplerview()
    moment0 = hd.moment0_dataset()
    moment2 = hd.moment2_dataset()
    spatial_shape = tuple(int(size) for size in moment0.shape[-2:])
    artery_mask, artery_swapped = _align_mask(
        dv.retinal_artery_mask(),
        spatial_shape,
        "retinal_artery_mask",
    )
    vein_mask, vein_swapped = _align_mask(
        dv.retinal_vein_mask(),
        spatial_shape,
        "retinal_vein_mask",
    )
    swapped = artery_swapped or vein_swapped
    center = _pair(dv.optic_disc_center(), swapped=swapped)
    width, height = _sizes(
        dv.optic_disc_width(),
        dv.optic_disc_height(),
        swapped=swapped,
    )
    return HeartbeatInputs(
        moment0=moment0,
        moment2=moment2,
        artery_mask=artery_mask,
        vein_mask=vein_mask,
        optic_disc_center=center,
        optic_disc_width=width,
        optic_disc_height=height,
        timing=hd.timing(),
        local_background_dist=dv.local_background_dist(),
    )


def _align_mask(array, shape: tuple[int, int], name: str) -> tuple[np.ndarray, bool]:
    value = np.asarray(array)
    if value.shape[-2:] == shape:
        return np.asarray(value, dtype=bool), False
    if value.shape[-2:] == shape[::-1]:
        return np.asarray(np.swapaxes(value, -1, -2), dtype=bool), True
    raise ValueError(f"{name} spatial shape {value.shape[-2:]} does not match {shape}.")


def _pair(value, *, swapped: bool) -> np.ndarray | None:
    if value is None:
        return None
    pair = np.asarray(value, dtype=np.float32).reshape(-1)
    if pair.size < 2:
        return None
    return pair[[1, 0]] if swapped else pair[:2]


def _sizes(width, height, *, swapped: bool):
    return (height, width) if swapped else (width, height)


__all__ = ["HeartbeatInputs", "load_heartbeat_inputs"]
