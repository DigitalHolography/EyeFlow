"""Neutral loading and alignment for vessel-topology inputs."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from calculations.topology import OpticDisc


@dataclass(frozen=True)
class VesselTopologyInputs:
    """The map reference and anatomy needed to prepare vessel topology."""

    moment0: object
    artery_mask: np.ndarray
    vein_mask: np.ndarray
    optic_disc: OpticDisc


def load_vessel_topology_inputs(ctx) -> VesselTopologyInputs:
    """Load one consistently oriented topology input set from HD and DV."""

    ctx.require_inputs("hd", "dv")
    hd = ctx.inputs.hd.as_holodoppler()
    dv = ctx.inputs.dv.as_dopplerview()
    moment0 = hd.moment0_dataset()
    spatial_shape = tuple(int(size) for size in moment0.shape[-2:])
    artery_mask, artery_swapped = align_mask(
        dv.retinal_artery_mask(),
        spatial_shape,
        "retinal_artery_mask",
    )
    vein_mask = apply_mask_alignment(
        dv.retinal_vein_mask(),
        spatial_shape,
        "retinal_vein_mask",
        swapped=artery_swapped,
    )
    optic_disc = dv.optic_disc()
    if artery_swapped:
        optic_disc = optic_disc.transposed()
    if optic_disc.mask is not None and optic_disc.mask.shape != spatial_shape:
        raise ValueError(
            "optic-disc mask must share the DopplerView vessel-mask frame; "
            f"expected {spatial_shape}, got {optic_disc.mask.shape}."
        )
    return VesselTopologyInputs(
        moment0=moment0,
        artery_mask=artery_mask,
        vein_mask=vein_mask,
        optic_disc=optic_disc,
    )


def align_mask(
    array,
    shape: tuple[int, int],
    name: str,
) -> tuple[np.ndarray, bool]:
    value = np.asarray(array)
    if value.shape[-2:] == shape:
        return np.asarray(value, dtype=bool), False
    if value.shape[-2:] == shape[::-1]:
        return np.asarray(np.swapaxes(value, -1, -2), dtype=bool), True
    raise ValueError(f"{name} spatial shape {value.shape[-2:]} does not match {shape}.")


def apply_mask_alignment(
    array,
    shape: tuple[int, int],
    name: str,
    *,
    swapped: bool,
) -> np.ndarray:
    value = np.asarray(array)
    aligned = np.swapaxes(value, -1, -2) if swapped else value
    if aligned.shape[-2:] != shape:
        raise ValueError(
            f"{name} does not share the selected DopplerView spatial orientation; "
            f"aligned shape {aligned.shape[-2:]} does not match {shape}."
        )
    return np.asarray(aligned, dtype=bool)


__all__ = [
    "VesselTopologyInputs",
    "align_mask",
    "apply_mask_alignment",
    "load_vessel_topology_inputs",
]
