"""Read and align canonical HD/DV retinal sources."""

from __future__ import annotations

import numpy as np

from app_settings import validate_velocity_estimation_method
from calculations.topology import OpticDisc
from retinal_source_data import (
    DopplerViewMetadata,
    ImageMaps,
    RetinalSegmentation,
    RetinalSourceData,
    VesselMasks,
)
from velocity_calibration import (
    DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
    validate_band_ratio_frequency_scale_hz,
)


def load_retinal_source_data(
    hd,
    dv,
    *,
    velocity_estimation_method: str = "doppler_moments",
    include_spectrograms: bool = True,
    band_ratio_frequency_scale_hz: float = (
        DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ
    ),
) -> RetinalSourceData:
    """Load canonical sources, optionally including spectra for direct callers.

    Runtime velocity/topology consumers omit spectra; the low-frequency
    spectral pipeline validates that optional input independently.
    """

    method = validate_velocity_estimation_method(velocity_estimation_method)
    frequency_scale_hz = validate_band_ratio_frequency_scale_hz(
        band_ratio_frequency_scale_hz
    )
    if method == "frequency_bands":
        band_lf, band_hf = hd.frequency_band_datasets()
        image_maps = ImageMaps(
            moment0=_optional_moment(hd.optional_moment0_dataset),
            moment2=_optional_moment(hd.optional_moment2_dataset),
            band_lf=band_lf,
            band_hf=band_hf,
            spectograms=hd.optional_spectograms_dataset() if include_spectrograms else None,
        )
        spatial_reference = band_lf
    else:
        image_maps = ImageMaps(
            moment0=hd.moment0_dataset(),
            moment2=hd.moment2_dataset(),
            spectograms=hd.optional_spectograms_dataset() if include_spectrograms else None,
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
    metadata = hd.metadata()
    disc_values = dv.optic_disc_measurements()
    optic_disc = OpticDisc.from_measurements(
        mask=disc_values.mask,
        center=disc_values.center,
        width=disc_values.width,
        height=disc_values.height,
        image_shape=dv_spatial_shape,
        pixel_pitch_mm=metadata.pixel_pitch.isotropic_mm,
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
    retinal_vessel_mask = apply_optional_mask_alignment(
        dv.retinal_vessel_mask(),
        spatial_shape,
        "retinal_vessel_segmentation_logits",
        swapped=artery_swapped,
    )
    choroid_vessel_mask = apply_optional_mask_alignment(
        dv.choroid_vessel_mask(),
        spatial_shape,
        "choroid_vessel_mask",
        swapped=artery_swapped,
    )
    return RetinalSourceData(
        image_maps=image_maps,
        segmentation=RetinalSegmentation(
            vessels=VesselMasks(
                artery=artery_mask,
                vein=source_vein_mask,
                labeled=labeled_vessels,
                velocity_background=(artery_mask | source_vein_mask),
                retinal_vessel_mask=retinal_vessel_mask,
                choroid_vessel_mask=choroid_vessel_mask,
            ),
            optic_disc=optic_disc,
        ),
        holodoppler=metadata,
        doppler_view=DopplerViewMetadata(
            local_background_dist=dv.local_background_dist(),
            spatial_axes_swapped_to_match_hd=artery_swapped,
        ),
        velocity_estimation_method=method,
        band_ratio_frequency_scale_hz=frequency_scale_hz,
    )


def align_spectrogram_frame(frame: np.ndarray) -> np.ndarray:
    """Map exported HD (freq, y, x) slices to EF lower-left (freq, x, y).

    HD's spatial Y direction is inverted relative to persisted EF gaps.
    Apply that explicit flip before the axis swap; never infer orientation
    from dimensions or reuse the DV alignment flag for this volume.
    """
    frame = np.asarray(frame)
    if frame.ndim != 3:
        raise ValueError("An HD spectrogram frame must have (freq, y, x) axes.")
    return np.swapaxes(np.flip(frame, axis=-2), -2, -1)


def iter_aligned_spectrogram_frames(spectrograms):
    """Read each lazy HD spectral frame once in the canonical EF orientation."""
    for index in range(spectrograms.shape[0]):
        yield align_spectrogram_frame(spectrograms[index])


def _optional_moment(loader):
    """Malformed inactive moments must not disable the band workflow."""
    try:
        return loader()
    except (KeyError, ValueError, TypeError):
        return None


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
    "align_spectrogram_frame",
    "iter_aligned_spectrogram_frames",
    "apply_mask_alignment",
    "apply_optional_mask_alignment",
    "apply_spatial_alignment",
    "load_retinal_source_data",
]
