"""Chunked retinal vessel-velocity estimation."""

from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass
from time import perf_counter

import numpy as np
from scipy import ndimage as ndi

from calculations.topology import annulus_mask
from utils.logger import Logger

SCRATCH_FRAME_CHUNK_SIZE = 32
SECTION_INNER_RADIUS_FRAC = 0.10
SECTION_OUTER_RADIUS_FRAC = 0.35
DEFAULT_LASER_WAVELENGTH_METERS = 8.52e-7
DEFAULT_NUMERICAL_APERTURE = 0.124
DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ = 1.0
DOPPLER_MOMENTS_METHOD = "doppler_moments"
FREQUENCY_BANDS_METHOD = "frequency_bands"
FREQUENCY_BAND_LF_PATH = "/band_0_3000_9000"
FREQUENCY_BAND_HF_PATH = "/band_1_9000_18000"


@dataclass(frozen=True)
class VelocityEstimatorCacheKey:
    """Run-local identity of every input that affects velocity estimation."""

    velocity_estimation_method: str
    moment0_source: tuple[object, ...] | None
    moment2_source: tuple[object, ...] | None
    band_lf_source: tuple[object, ...] | None
    band_hf_source: tuple[object, ...] | None
    artery_mask: tuple[object, ...]
    vein_mask: tuple[object, ...]
    background_mask: tuple[object, ...]
    optic_disc_center: tuple[object, ...]
    section_inner_radius_frac: float
    section_outer_radius_frac: float
    local_background_dist: int
    laser_wavelength: float
    numerical_aperture: float
    band_ratio_frequency_scale_hz: float | None
    frame_chunk_size: int


def velocity_estimator_cache_key(
    *,
    moment0=None,
    moment2=None,
    band_lf=None,
    band_hf=None,
    velocity_estimation_method: str = DOPPLER_MOMENTS_METHOD,
    artery_mask,
    vein_mask,
    background_mask=None,
    optic_disc_center,
    section_inner_radius_frac: float = SECTION_INNER_RADIUS_FRAC,
    section_outer_radius_frac: float = SECTION_OUTER_RADIUS_FRAC,
    local_background_dist: int,
    laser_wavelength: float = DEFAULT_LASER_WAVELENGTH_METERS,
    numerical_aperture: float = DEFAULT_NUMERICAL_APERTURE,
    band_ratio_frequency_scale_hz: float = (
        DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ
    ),
) -> VelocityEstimatorCacheKey:
    """Build a conservative key for run-scoped estimator-result reuse."""

    method, first_volume, second_volume = _active_velocity_volumes(
        velocity_estimation_method=velocity_estimation_method,
        moment0=moment0,
        moment2=moment2,
        band_lf=band_lf,
        band_hf=band_hf,
    )
    _validate_matching_volumes(method, first_volume, second_volume)
    artery = np.asarray(artery_mask, dtype=bool)
    vein = np.asarray(vein_mask, dtype=bool)
    background = (
        artery | vein
        if background_mask is None
        else np.asarray(background_mask, dtype=bool)
    )
    return VelocityEstimatorCacheKey(
        velocity_estimation_method=method,
        moment0_source=(
            _volume_source_key(first_volume)
            if method == DOPPLER_MOMENTS_METHOD
            else None
        ),
        moment2_source=(
            _volume_source_key(second_volume)
            if method == DOPPLER_MOMENTS_METHOD
            else None
        ),
        band_lf_source=(
            _volume_source_key(first_volume)
            if method == FREQUENCY_BANDS_METHOD
            else None
        ),
        band_hf_source=(
            _volume_source_key(second_volume)
            if method == FREQUENCY_BANDS_METHOD
            else None
        ),
        artery_mask=_array_value_key(artery, dtype=bool),
        vein_mask=_array_value_key(vein, dtype=bool),
        background_mask=_array_value_key(background, dtype=bool),
        optic_disc_center=_array_value_key(
            optic_disc_center,
            dtype=np.float32,
        ),
        section_inner_radius_frac=float(section_inner_radius_frac),
        section_outer_radius_frac=float(section_outer_radius_frac),
        local_background_dist=int(local_background_dist),
        laser_wavelength=float(laser_wavelength),
        numerical_aperture=float(numerical_aperture),
        band_ratio_frequency_scale_hz=(
            _validate_band_ratio_frequency_scale(
                band_ratio_frequency_scale_hz
            )
            if method == FREQUENCY_BANDS_METHOD
            else None
        ),
        frame_chunk_size=SCRATCH_FRAME_CHUNK_SIZE,
    )


def _velocity_from_delta_frequency(
    delta_frequency,
    laser_wavelength: float = DEFAULT_LASER_WAVELENGTH_METERS,
    numerical_aperture: float = DEFAULT_NUMERICAL_APERTURE,
) -> np.ndarray:
    """Convert a Doppler-frequency shift in Hz to velocity in mm/s."""
    delta_frequency = np.asarray(delta_frequency, dtype=np.float32)
    return (
        np.float32(1e3) * laser_wavelength * delta_frequency / numerical_aperture
    ).astype(np.float32, copy=False)


def run_chunked_velocity_estimator(
    *,
    moment0=None,
    moment2=None,
    band_lf=None,
    band_hf=None,
    velocity_estimation_method: str = DOPPLER_MOMENTS_METHOD,
    artery_mask,
    vein_mask,
    background_mask=None,
    optic_disc_center,
    section_inner_radius_frac: float = SECTION_INNER_RADIUS_FRAC,
    section_outer_radius_frac: float = SECTION_OUTER_RADIUS_FRAC,
    local_background_dist: int,
    scratch_h5,
    laser_wavelength: float = DEFAULT_LASER_WAVELENGTH_METERS,
    numerical_aperture: float = DEFAULT_NUMERICAL_APERTURE,
    band_ratio_frequency_scale_hz: float = (
        DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ
    ),
    retain_velocity_video: bool = True,
    velocity_video_output=None,
) -> dict[str, object]:
    """Estimate velocity into scratch datasets without materializing full videos."""

    method, first_volume, second_volume = _active_velocity_volumes(
        velocity_estimation_method=velocity_estimation_method,
        moment0=moment0,
        moment2=moment2,
        band_lf=band_lf,
        band_hf=band_hf,
    )
    _validate_matching_volumes(method, first_volume, second_volume)
    frequency_scale_hz = (
        _validate_band_ratio_frequency_scale(
            band_ratio_frequency_scale_hz
        )
        if method == FREQUENCY_BANDS_METHOD
        else None
    )
    frame_count, height, width = (int(size) for size in first_volume.shape)
    artery = np.asarray(artery_mask, dtype=bool)
    vein = np.asarray(vein_mask, dtype=bool)
    background = (
        artery | vein
        if background_mask is None
        else np.asarray(background_mask, dtype=bool)
    )
    if (
        artery.shape != (height, width)
        or vein.shape != (height, width)
        or background.shape != (height, width)
    ):
        raise ValueError(
            "Velocity masks must match the active HoloDoppler estimator "
            f"spatial shape {(height, width)}."
        )

    if method == DOPPLER_MOMENTS_METHOD:
        Logger.log("Velocity estimator uses raw HD moments.")
    else:
        Logger.log(
            "Velocity estimator converts the HoloDoppler band ratio "
            f"{FREQUENCY_BAND_HF_PATH} / {FREQUENCY_BAND_LF_PATH} to RMS "
            f"frequency using {frequency_scale_hz:g} Hz per ratio unit."
        )
    Logger.log(
        f"Velocity estimator uses {SCRATCH_FRAME_CHUNK_SIZE}-frame batched "
        "inpainting and summary-only frequency intermediates."
    )

    velocity_dataset = _velocity_video_storage(
        scratch_h5,
        (frame_count, height, width),
        retain_velocity_video=retain_velocity_video,
        velocity_video_output=velocity_video_output,
    )
    disk, inpaint = _skimage_dependencies()
    inpaint_mask = _dilated_mask(background, disk(int(local_background_dist)))
    section_mask = annulus_mask(
        (height, width),
        optic_disc_center,
        section_inner_radius_frac,
        section_outer_radius_frac,
    )
    artery_section = section_mask & artery
    vein_section = section_mask & vein

    averages = {
        name: np.zeros((height, width), dtype=np.float64)
        for name in ("background", "velocity", "fRMS", "fRMS_bkg", "deltafRMS")
    }
    signals = {
        name: np.full(frame_count, np.nan, dtype=np.float32)
        for name in (
            "artery_velocity",
            "vein_velocity",
            "artery_fRMS",
            "vein_fRMS",
            "artery_fRMS_bkg",
            "vein_fRMS_bkg",
            "vessel_fRMS_bkg",
            "artery_deltafRMS",
            "vein_deltafRMS",
        )
    }

    estimation_started = perf_counter()
    chunk_count = max(1, (frame_count + SCRATCH_FRAME_CHUNK_SIZE - 1) // SCRATCH_FRAME_CHUNK_SIZE)
    for chunk_index, start in enumerate(range(0, frame_count, SCRATCH_FRAME_CHUNK_SIZE), start=1):
        stop = min(start + SCRATCH_FRAME_CHUNK_SIZE, frame_count)
        frame_slice = slice(start, stop)
        if method == DOPPLER_MOMENTS_METHOD:
            background_image = _read_volume_chunk(first_volume, frame_slice)
            moment2_chunk = _read_volume_chunk(second_volume, frame_slice)
            mean_m0 = np.mean(
                background_image,
                axis=(-1, -2),
                keepdims=True,
                dtype=np.float32,
            )
            f_rms = np.sqrt(
                np.divide(
                    moment2_chunk,
                    mean_m0,
                    out=np.zeros_like(moment2_chunk, dtype=np.float32),
                    where=mean_m0 != 0,
                )
            ).astype(np.float32, copy=False)
        else:
            background_image = _read_validated_band_chunk(
                first_volume,
                frame_slice,
                FREQUENCY_BAND_LF_PATH,
            )
            high_frequency = _read_validated_band_chunk(
                second_volume,
                frame_slice,
                FREQUENCY_BAND_HF_PATH,
            )
            f_rms = _band_ratio_to_frequency(
                _safe_band_ratio(
                    high_frequency,
                    background_image,
                    frame_slice=frame_slice,
                ),
                frequency_scale_hz=frequency_scale_hz,
                frame_slice=frame_slice,
            )
        f_rms_background = _inpaint_frame_batch(
            f_rms,
            inpaint_mask,
            inpaint,
        )
        delta = _signed_rms_difference(f_rms, f_rms_background)
        velocity = _velocity_from_delta_frequency(
            delta,
            laser_wavelength=laser_wavelength,
            numerical_aperture=numerical_aperture,
        )

        if velocity_dataset is not None:
            velocity_dataset[frame_slice] = velocity
        averages["background"] += np.sum(
            background_image,
            axis=0,
            dtype=np.float64,
        )
        averages["velocity"] += np.sum(velocity, axis=0, dtype=np.float64)
        averages["fRMS"] += np.sum(f_rms, axis=0, dtype=np.float64)
        averages["fRMS_bkg"] += np.sum(f_rms_background, axis=0, dtype=np.float64)
        averages["deltafRMS"] += np.sum(delta, axis=0, dtype=np.float64)
        signals["artery_velocity"][frame_slice] = _masked_signal(
            velocity,
            artery_section,
        )
        signals["vein_velocity"][frame_slice] = _masked_signal(
            velocity,
            vein_section,
        )
        signals["artery_fRMS"][frame_slice] = _masked_signal(f_rms, artery_section)
        signals["vein_fRMS"][frame_slice] = _masked_signal(f_rms, vein_section)
        signals["artery_fRMS_bkg"][frame_slice] = _masked_signal(
            f_rms_background,
            artery_section,
        )
        signals["vein_fRMS_bkg"][frame_slice] = _masked_signal(
            f_rms_background,
            vein_section,
        )
        signals["vessel_fRMS_bkg"][frame_slice] = _masked_signal(
            f_rms_background,
            artery_section | vein_section,
        )
        signals["artery_deltafRMS"][frame_slice] = _masked_signal(
            delta,
            artery_section,
        )
        signals["vein_deltafRMS"][frame_slice] = _masked_signal(
            delta,
            vein_section,
        )
        if chunk_index == chunk_count or chunk_index % 10 == 0:
            Logger.log(
                f"Velocity estimation completed chunk {chunk_index}/{chunk_count} "
                f"({stop}/{frame_count} frames)."
            )

    Logger.log(
        f"Completed chunked velocity estimation in {perf_counter() - estimation_started:.1f}s."
    )

    divisor = np.float64(max(frame_count, 1))
    background_average = (averages["background"] / divisor).astype(np.float32)
    return {
        "fRMS": None,
        "fRMS_bkg": None,
        "deltafRMS": None,
        "velocity_map": velocity_dataset,
        "retinal_vessel_velocity": velocity_dataset,
        # Preserve the legacy key consumed by existing figures. In band mode,
        # it intentionally contains the LF temporal mean used as the display
        # background rather than requiring moment0 only for visualization.
        "moment0_avg": background_average,
        "background_image_avg": background_average,
        "velocity_map_avg": (averages["velocity"] / divisor).astype(np.float32),
        "fRMS_avg": (averages["fRMS"] / divisor).astype(np.float32),
        "fRMS_bkg_avg": (averages["fRMS_bkg"] / divisor).astype(np.float32),
        "deltafRMS_avg": (averages["deltafRMS"] / divisor).astype(np.float32),
        "velocity_section_mask": section_mask,
        "velocity_section_geometry": "optic_disc_centered_frame_fraction",
        "retinal_artery_velocity_signal": signals["artery_velocity"],
        "retinal_vein_velocity_signal": signals["vein_velocity"],
        "retinal_artery_fRMS_signal": signals["artery_fRMS"],
        "retinal_vein_fRMS_signal": signals["vein_fRMS"],
        "retinal_artery_fRMS_bkg_signal": signals["artery_fRMS_bkg"],
        "retinal_vein_fRMS_bkg_signal": signals["vein_fRMS_bkg"],
        "retinal_vessel_fRMS_bkg_signal": signals["vessel_fRMS_bkg"],
        "retinal_artery_deltafRMS_signal": signals["artery_deltafRMS"],
        "retinal_vein_deltafRMS_signal": signals["vein_deltafRMS"],
        "velocity_estimation_method": method,
        "velocity_quantity": "physical_velocity",
        "velocity_unit": "mm/s",
        "band_ratio_frequency_scale_hz": frequency_scale_hz,
        "band_lf_source_path": (
            FREQUENCY_BAND_LF_PATH if method == FREQUENCY_BANDS_METHOD else None
        ),
        "band_hf_source_path": (
            FREQUENCY_BAND_HF_PATH if method == FREQUENCY_BANDS_METHOD else None
        ),
    }


def _velocity_video_storage(
    scratch_h5,
    shape: tuple[int, int, int],
    *,
    retain_velocity_video: bool,
    velocity_video_output,
):
    group = scratch_h5.require_group("waveform")
    if not retain_velocity_video:
        if velocity_video_output is not None:
            raise ValueError(
                "velocity_video_output requires retain_velocity_video=True."
            )
        return None
    if velocity_video_output is not None:
        if tuple(velocity_video_output.shape) != shape:
            raise ValueError(
                "velocity_video_output must match the active estimator volume shape."
            )
        if np.dtype(velocity_video_output.dtype) != np.dtype(np.float32):
            raise ValueError("velocity_video_output must have dtype float32.")
        return velocity_video_output

    return group.create_dataset(
        "velocity",
        shape=shape,
        dtype=np.float32,
        chunks=(
            min(64, shape[0]),
            min(32, shape[1]),
            min(32, shape[2]),
        ),
        compression=None,
    )


def _active_velocity_volumes(
    *,
    velocity_estimation_method: str,
    moment0,
    moment2,
    band_lf,
    band_hf,
) -> tuple[str, object, object]:
    method = str(velocity_estimation_method)
    if method == DOPPLER_MOMENTS_METHOD:
        missing = [
            name
            for name, value in (("moment0", moment0), ("moment2", moment2))
            if value is None
        ]
        if missing:
            raise ValueError(
                "velocity_estimation_method='doppler_moments' requires "
                f"{', '.join(missing)}."
            )
        return method, moment0, moment2
    if method == FREQUENCY_BANDS_METHOD:
        missing = [
            path
            for path, value in (
                (FREQUENCY_BAND_LF_PATH, band_lf),
                (FREQUENCY_BAND_HF_PATH, band_hf),
            )
            if value is None
        ]
        if missing:
            raise ValueError(
                "velocity_estimation_method='frequency_bands' requires "
                f"HoloDoppler datasets {', '.join(missing)}."
            )
        return method, band_lf, band_hf
    raise ValueError(
        "velocity_estimation_method must be 'doppler_moments' or "
        f"'frequency_bands', got {velocity_estimation_method!r}."
    )


def _validate_matching_volumes(method: str, first_volume, second_volume) -> None:
    if method == DOPPLER_MOMENTS_METHOD:
        first_name, second_name = "moment0", "moment2"
    else:
        first_name, second_name = FREQUENCY_BAND_LF_PATH, FREQUENCY_BAND_HF_PATH
    for name, volume in (
        (first_name, first_volume),
        (second_name, second_volume),
    ):
        shape = tuple(int(size) for size in getattr(volume, "shape", ()))
        if len(shape) != 3:
            raise ValueError(
                f"{name} must be a 3-D (frame, y, x) dataset, got shape {shape}."
            )
        dtype = getattr(volume, "dtype", None)
        if dtype is None or not np.issubdtype(np.dtype(dtype), np.number):
            raise TypeError(f"{name} must be numeric, got dtype {dtype}.")
    if tuple(first_volume.shape) != tuple(second_volume.shape):
        raise ValueError(
            f"{first_name} and {second_name} must have identical "
            f"(frame, y, x) shapes, got {first_volume.shape} and "
            f"{second_volume.shape}."
        )


def _volume_source_key(value) -> tuple[object, ...]:
    shape = tuple(int(size) for size in value.shape)
    dtype = np.dtype(value.dtype).str
    try:
        filename = value.file.filename
        dataset_name = value.name
    except (AttributeError, RuntimeError, ValueError):
        filename = None
        dataset_name = None
    if filename is not None and dataset_name is not None:
        return (
            "hdf5",
            os.path.normcase(os.path.abspath(str(filename))),
            str(dataset_name),
            shape,
            dtype,
        )
    return ("run_object", id(value), shape, dtype)


def _optional_array_value_key(value, *, dtype) -> tuple[object, ...] | None:
    return None if value is None else _array_value_key(value, dtype=dtype)


def _array_value_key(value, *, dtype) -> tuple[object, ...]:
    array = np.ascontiguousarray(np.asarray(value, dtype=dtype))
    digest = hashlib.blake2b(
        array.view(np.uint8),
        digest_size=16,
    ).hexdigest()
    return (array.shape, array.dtype.str, digest)


def _read_volume_chunk(volume, frame_slice: slice) -> np.ndarray:
    return np.asarray(volume[frame_slice], dtype=np.float32)


def _read_validated_band_chunk(
    volume,
    frame_slice: slice,
    dataset_path: str,
) -> np.ndarray:
    values = _read_volume_chunk(volume, frame_slice)
    if not np.all(np.isfinite(values)):
        raise ValueError(
            f"HoloDoppler frequency-band dataset {dataset_path} contains "
            f"non-finite values in frames {frame_slice.start}:{frame_slice.stop}."
        )
    if np.any(values < 0.0):
        raise ValueError(
            f"HoloDoppler frequency-band dataset {dataset_path} contains "
            f"negative values in frames {frame_slice.start}:{frame_slice.stop}."
        )
    return values


def _safe_band_ratio(
    high_frequency: np.ndarray,
    low_frequency: np.ndarray,
    *,
    frame_slice: slice,
) -> np.ndarray:
    """Return HF/LF with exact-zero LF mapped to zero and no epsilon bias."""

    high64 = np.asarray(high_frequency, dtype=np.float64)
    low64 = np.asarray(low_frequency, dtype=np.float64)
    ratio64 = np.divide(
        high64,
        low64,
        out=np.zeros_like(high64),
        where=low64 != 0.0,
    )
    float32_limit = np.finfo(np.float32).max
    if not np.all(np.isfinite(ratio64)) or np.any(ratio64 > float32_limit):
        raise ValueError(
            "The HoloDoppler HF/LF ratio is outside the finite float32 "
            f"range in frames {frame_slice.start}:{frame_slice.stop}; no "
            "infinite RMS-frequency estimate was produced."
        )
    return ratio64.astype(np.float32, copy=False)


def _validate_band_ratio_frequency_scale(value: float) -> float:
    """Return a finite, positive Hz-per-ratio calibration factor."""

    scale_hz = float(value)
    if not np.isfinite(scale_hz) or scale_hz <= 0.0:
        raise ValueError(
            "band_ratio_frequency_scale_hz must be a finite positive value "
            "in Hz per ratio unit."
        )
    return scale_hz


def _band_ratio_to_frequency(
    ratio: np.ndarray,
    *,
    frequency_scale_hz: float,
    frame_slice: slice,
) -> np.ndarray:
    """Convert the dimensionless HF/LF ratio to an RMS frequency in Hz."""

    frequency64 = np.asarray(ratio, dtype=np.float64) * frequency_scale_hz
    float32_limit = np.finfo(np.float32).max
    if not np.all(np.isfinite(frequency64)) or np.any(
        frequency64 > float32_limit
    ):
        raise ValueError(
            "The calibrated HoloDoppler band-ratio frequency is outside the "
            "finite float32 range in frames "
            f"{frame_slice.start}:{frame_slice.stop}."
        )
    return frequency64.astype(np.float32, copy=False)


def _inpaint_frame_batch(
    frames: np.ndarray,
    mask: np.ndarray,
    inpaint,
) -> np.ndarray:
    """Inpaint all frames as channels so the sparse system is built once."""
    source = np.asarray(frames, dtype=np.float32)
    channels_last = np.moveaxis(source, 0, -1)
    result = inpaint.inpaint_biharmonic(
        channels_last,
        np.asarray(mask, dtype=bool),
        channel_axis=-1,
    )
    inpainted = np.moveaxis(result, -1, 0).astype(np.float32, copy=False)
    square_safe_limit = np.float32(np.sqrt(np.finfo(np.float32).max))
    valid_frames = np.all(np.isfinite(inpainted), axis=(-1, -2)) & np.all(
        np.abs(inpainted) <= square_safe_limit,
        axis=(-1, -2),
    )
    if np.all(valid_frames):
        return inpainted

    output = inpainted.copy()
    output[~valid_frames] = _bounded_inpaint_result(
        inpainted[~valid_frames],
        source[~valid_frames],
        mask,
    )
    return output


def _bounded_inpaint_result(
    inpainted: np.ndarray,
    source: np.ndarray,
    mask: np.ndarray,
) -> np.ndarray:
    """Keep each inpainted frame within its finite background range."""
    background = np.asarray(source, dtype=np.float32)[:, ~np.asarray(mask, dtype=bool)]
    finite = np.isfinite(background)
    if background.shape[1] == 0 or np.any(~np.any(finite, axis=1)):
        raise ValueError("Background inpainting requires a finite unmasked pixel per frame.")

    finite_background = np.where(finite, background, np.nan)
    lower = np.nanmin(finite_background, axis=1)
    upper = np.nanmax(finite_background, axis=1)
    fallback = np.nanmean(finite_background, axis=1, dtype=np.float32)
    bounded = np.where(
        np.isfinite(inpainted),
        inpainted,
        fallback[:, None, None],
    )
    return np.clip(
        bounded,
        lower[:, None, None],
        upper[:, None, None],
    ).astype(np.float32, copy=False)


def _signed_rms_difference(
    f_rms: np.ndarray,
    f_rms_background: np.ndarray,
) -> np.ndarray:
    """Return signed root-square difference without float32 square overflow.
    Reason: Encountered overflow.
    """
    foreground64 = np.asarray(f_rms, dtype=np.float64)
    background64 = np.asarray(f_rms_background, dtype=np.float64)
    with np.errstate(invalid="ignore"):
        difference64 = np.square(foreground64) - np.square(background64)
        delta64 = np.sign(difference64) * np.sqrt(np.abs(difference64))
    return delta64.astype(np.float32, copy=False)


def _skimage_dependencies():
    try:
        from skimage.morphology import disk
        from skimage.restoration import inpaint
    except ModuleNotFoundError as exc:
        raise ImportError(
            "Retinal velocity estimation requires scikit-image."
        ) from exc
    return disk, inpaint


def _dilated_mask(vessel_mask: np.ndarray, footprint: np.ndarray) -> np.ndarray:
    mask = np.asarray(vessel_mask, dtype=bool)
    if mask.ndim != 2:
        raise ValueError(f"vessel_mask must be 2-D for dilation, got {mask.shape}.")
    return ndi.binary_dilation(mask, structure=np.asarray(footprint, dtype=bool))


def _masked_signal(velocity_map: np.ndarray, mask: np.ndarray) -> np.ndarray:
    selected = velocity_map[:, np.asarray(mask, dtype=bool)]
    if not np.any(np.isfinite(selected)):
        return np.full((velocity_map.shape[0],), np.nan, dtype=np.float32)
    return np.nanmean(selected, axis=1, dtype=np.float64).astype(
        np.float32,
        copy=False,
    )
