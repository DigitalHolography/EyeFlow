"""Segmentation inputs owned by the spatial-gradient pipeline."""

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class SpatialGradientSource:
    retinal_artery_mask: np.ndarray
    retinal_vein_mask: np.ndarray
    optic_disc_mask: np.ndarray | None
    optic_disc_center: object
    optic_disc_width: object
    optic_disc_height: object
    window_size_percentile_kept: float = .95


def load_spatial_gradient_source(ctx, spatial_shape) -> SpatialGradientSource:
    ctx.require_inputs("hd", "dv")
    dv = ctx.inputs.dv.as_dopplerview()
    swapped = False

    def align(value, name):
        nonlocal swapped
        if value is None:
            return None
        array = np.asarray(value, dtype=bool)
        if array.shape == tuple(spatial_shape):
            return array
        if array.shape == tuple(spatial_shape)[::-1]:
            swapped = True
            return array.T
        raise ValueError(f"{name} shape {array.shape} does not match gradient shape {spatial_shape}.")

    artery = align(dv.retinal_artery_mask(), "artery mask")
    vein = align(dv.retinal_vein_mask(), "vein mask")
    # Match the shared input adapters' coordinate-transposition contract.
    labeled = dv.retinal_labeled_vessels()
    if labeled is not None and np.asarray(labeled).shape != tuple(spatial_shape):
        if np.asarray(labeled).shape != tuple(spatial_shape)[::-1]:
            raise ValueError("labeled vessel shape does not match gradient shape.")
        swapped = True
    disc = align(dv.optic_disc_mask(), "optic disc mask")
    center = dv.optic_disc_center()
    width, height = dv.optic_disc_width(), dv.optic_disc_height()
    if swapped:
        if center is not None:
            center = np.asarray(center).copy()
            center.reshape(-1)[:2] = center.reshape(-1)[:2][::-1]
        width, height = height, width
    return SpatialGradientSource(artery, vein, disc, center, width, height)
