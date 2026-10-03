"""Load and align the canonical retinal source-data model."""

from __future__ import annotations

import numpy as np

from app_settings import validate_velocity_estimation_method

from input_output.schema import (
    DopplerViewMetadata,
    ImageMaps,
    RetinalSegmentation,
    RetinalSourceData,
    VesselMasks,
)


def load_vessel_topology_inputs(ctx) -> RetinalSourceData:
    """Load one consistently oriented topology input set from HD and DV."""

    ctx.require_inputs("hd", "dv")
    hd = ctx.inputs.hd.as_holodoppler()
    dv = ctx.inputs.dv.as_dopplerview()
    method = validate_velocity_estimation_method(
        getattr(ctx, "velocity_estimation_method", "doppler_moments")
    )
    return load_retinal_source_data(
        hd,
        dv,
        velocity_estimation_method=method,
    )


def load_retinal_source_data(
    hd,
    dv,
    *,
    velocity_estimation_method: str = "doppler_moments",
) -> RetinalSourceData:
    """Load one canonical source model from typed HD and DV adapters."""

    method = validate_velocity_estimation_method(velocity_estimation_method)
    if method == "frequency_bands":
        band_lf, band_hf = hd.frequency_band_datasets()
        image_maps = ImageMaps(
            moment0=hd.optional_moment0_dataset(),
            moment2=hd.optional_moment2_dataset(),
            band_lf=band_lf,
            band_hf=band_hf,
        )
        spatial_reference = band_lf
    else:
        image_maps = ImageMaps(
            moment0=hd.moment0_dataset(),
            moment2=hd.moment2_dataset(),
        )
        spatial_reference = image_maps.moment0
    spatial_shape = tuple(int(size) for size in spatial_reference.shape[-2:])
    artery_mask, artery_swapped = align_mask(
        dv.retinal_artery_mask(),
        spatial_shape,
        "retinal_artery_mask",
    )
    source_vein_mask = apply_mask_alignment(
        dv.retinal_vein_mask(),
        spatial_shape,
        "retinal_vein_mask",
        swapped=artery_swapped,
    )
    dv_spatial_shape = spatial_shape[::-1] if artery_swapped else spatial_shape
    optic_disc = dv.optic_disc(dv_spatial_shape)
    if artery_swapped:
        optic_disc = optic_disc.transposed()
    if optic_disc.mask is not None and optic_disc.mask.shape != spatial_shape:
        raise ValueError(
            "optic-disc mask must share the DopplerView vessel-mask frame; "
            f"expected {spatial_shape}, got {optic_disc.mask.shape}."
        )
    labeled_vessels = apply_optional_mask_alignment(
        dv.retinal_labeled_vessels(),
        spatial_shape,
        "retinal_labeled_vessels",
        swapped=artery_swapped,
    )
    # A fallback disc means DopplerView did not provide usable optic-disc
    # detection. Venous analysis is not meaningful without that reference,
    # so expose an empty processing mask. Keep the original vein mask only in
    # the shared velocity-background support so arterial inpainting remains
    # exactly as it was before this guard was introduced.
    vein_mask = (
        np.zeros_like(source_vein_mask, dtype=bool)
        if optic_disc.is_fallback
        else source_vein_mask
    )
    return RetinalSourceData(
        image_maps=image_maps,
        segmentation=RetinalSegmentation(
            vessels=VesselMasks(
                artery=artery_mask,
                vein=vein_mask,
                labeled=labeled_vessels,
                velocity_background=(artery_mask | source_vein_mask),
            ),
            optic_disc=optic_disc,
        ),
        holodoppler=hd.metadata(),
        doppler_view=DopplerViewMetadata(
            local_background_dist=dv.local_background_dist(),
            spatial_axes_swapped_to_match_hd=artery_swapped,
        ),
        velocity_estimation_method=method,
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
    return np.asarray(
        apply_spatial_alignment(array, shape, name, swapped=swapped),
        dtype=bool,
    )


def apply_spatial_alignment(
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
    return aligned


def apply_optional_mask_alignment(
    array,
    shape: tuple[int, int],
    name: str,
    *,
    swapped: bool,
) -> np.ndarray | None:
    if array is None:
        return None
    return apply_spatial_alignment(array, shape, name, swapped=swapped)


__all__ = [
    "align_mask",
    "apply_mask_alignment",
    "apply_optional_mask_alignment",
    "apply_spatial_alignment",
    "load_retinal_source_data",
    "load_vessel_topology_inputs",
]
