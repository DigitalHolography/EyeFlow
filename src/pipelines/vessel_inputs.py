"""Load and align the canonical retinal source-data model."""

from __future__ import annotations

import numpy as np

from calculations.topology import OpticDisc

from input_output.spatial_alignment import (
    align_mask,
    apply_mask_alignment,
    apply_optional_mask_alignment,
    apply_spatial_alignment,
)
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
    return load_retinal_source_data(hd, dv)


def load_retinal_source_data(hd, dv) -> RetinalSourceData:
    """Load one canonical source model from typed HD and DV adapters."""

    image_maps = ImageMaps(
        moment0=hd.moment0_dataset(),
        moment2=hd.moment2_dataset(),
    )
    spatial_shape = tuple(int(size) for size in image_maps.moment0.shape[-2:])
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
    measurements = dv.optic_disc_measurements()
    optic_disc = OpticDisc.from_measurements(
        measurements.mask,
        measurements.center,
        measurements.width,
        measurements.height,
        dv_spatial_shape,
    )
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
    )




__all__ = [
    "align_mask",
    "apply_mask_alignment",
    "apply_optional_mask_alignment",
    "apply_spatial_alignment",
    "load_retinal_source_data",
    "load_vessel_topology_inputs",
]
