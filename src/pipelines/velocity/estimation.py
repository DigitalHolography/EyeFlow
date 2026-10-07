"""Chunked retinal vessel-velocity estimation."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass, field
from functools import singledispatch
from time import perf_counter

import numpy as np
from scipy import ndimage as ndi

from calculations.topology import annulus_mask
from input_output.schema import HD_BAND_HF_PATH, HD_BAND_LF_PATH, ImageMaps
from utils.logger import Logger
from velocity_calibration import (
    BAND_LF_LOW_RELATIVE_THRESHOLD,
    DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
    DEFAULT_LASER_WAVELENGTH_METERS,
    DEFAULT_NUMERICAL_APERTURE,
    physical_velocity_provenance,
    validate_band_ratio_frequency_scale_hz,
)

from .models import (
    RetinalVelocityData,
    RetinalVelocityMaps,
    VesselVelocitySignals,
)

SCRATCH_FRAME_CHUNK_SIZE = 32
SECTION_INNER_RADIUS_FRAC = 0.10
SECTION_OUTER_RADIUS_FRAC = 0.35
DOPPLER_MOMENTS_METHOD = "doppler_moments"
FREQUENCY_BANDS_METHOD = "frequency_bands"
FREQUENCY_BAND_LF_PATH = f"/{HD_BAND_LF_PATH}"
FREQUENCY_BAND_HF_PATH = f"/{HD_BAND_HF_PATH}"


@dataclass(frozen=True, slots=True)
class VelocityEstimatorInputs:
    """Validated source volumes shared by every velocity estimator."""

    method: str
    first_volume: object
    second_volume: object


@dataclass(frozen=True, slots=True)
class DopplerMomentsEstimatorInputs(VelocityEstimatorInputs):
    """Moment volumes used by the Doppler-moments estimator."""


@dataclass(frozen=True, slots=True)
class FrequencyBandsEstimatorInputs(VelocityEstimatorInputs):
    """Band volumes and calibration used by the frequency-band estimator."""

    frequency_scale_hz: float


@dataclass(frozen=True, slots=True)
class VelocityEstimatorChunk:
    """Standard method-independent calculation returned for one frame slice."""

    chunk_index: int
    chunk_count: int
    frame_slice: slice
    background_image: np.ndarray
    rms_frequency: np.ndarray
    diagnostic_counts: Mapping[str, int] = field(default_factory=dict)


@dataclass(slots=True)
class _VelocityResultAccumulator:
    """Accumulate retained maps and vessel signals across estimator chunks."""

    frame_count: int
    section_mask: np.ndarray
    artery_section: np.ndarray
    vein_section: np.ndarray
    velocity_video: np.ndarray | None
    averages: dict[str, np.ndarray]
    signals: dict[str, np.ndarray]

    @classmethod
    def create(
        cls,
        *,
        frame_count: int,
        spatial_shape: tuple[int, int],
        section_mask: np.ndarray,
        artery_mask: np.ndarray,
        vein_mask: np.ndarray,
        retain_velocity_video: bool,
    ) -> _VelocityResultAccumulator:
        averages = {
            name: np.zeros(spatial_shape, dtype=np.float64)
            for name in (
                "background",
                "velocity",
                "velocity_unmasked",
                "frms",
                "frms_background",
                "delta_frms",
            )
        }
        signals = {
            name: np.full(frame_count, np.nan, dtype=np.float32)
            for name in (
                "artery_velocity",
                "vein_velocity",
                "artery_frms",
                "vein_frms",
                "artery_frms_background",
                "vein_frms_background",
                "vessel_frms_background",
                "artery_delta_frms",
                "vein_delta_frms",
            )
        }
        return cls(
            frame_count=frame_count,
            section_mask=section_mask,
            artery_section=section_mask & artery_mask,
            vein_section=section_mask & vein_mask,
            velocity_video=(
                np.empty((frame_count, *spatial_shape), dtype=np.float32)
                if retain_velocity_video
                else None
            ),
            averages=averages,
            signals=signals,
        )

    def update(
        self,
        estimator_chunk: VelocityEstimatorChunk,
        *,
        rms_frequency_background: np.ndarray,
        delta_rms_frequency: np.ndarray,
        velocity: np.ndarray,
        velocity_unmasked: np.ndarray,
    ) -> None:
        frame_slice = estimator_chunk.frame_slice
        if self.velocity_video is not None:
            self.velocity_video[frame_slice] = velocity

        for name, values in (
            ("background", estimator_chunk.background_image),
            ("velocity", velocity),
            ("velocity_unmasked", velocity_unmasked),
            ("frms", estimator_chunk.rms_frequency),
            ("frms_background", rms_frequency_background),
            ("delta_frms", delta_rms_frequency),
        ):
            self.averages[name] += np.sum(values, axis=0, dtype=np.float64)

        for vessel, section in (
            ("artery", self.artery_section),
            ("vein", self.vein_section),
        ):
            self.signals[f"{vessel}_velocity"][frame_slice] = _masked_signal(
                velocity,
                section,
            )
            self.signals[f"{vessel}_frms"][frame_slice] = _masked_signal(
                estimator_chunk.rms_frequency,
                section,
            )
            self.signals[f"{vessel}_frms_background"][frame_slice] = _masked_signal(
                rms_frequency_background,
                section,
            )
            self.signals[f"{vessel}_delta_frms"][frame_slice] = _masked_signal(
                delta_rms_frequency,
                section,
            )
        self.signals["vessel_frms_background"][frame_slice] = _masked_signal(
            rms_frequency_background,
            self.artery_section | self.vein_section,
        )

    def build_result(
        self,
        *,
        provenance: Mapping[str, object],
    ) -> RetinalVelocityData:
        return RetinalVelocityData(
            maps=RetinalVelocityMaps(
                velocity=self.velocity_video,
                # In band mode this is the LF mean used as the display background.
                moment0_average=self._average("background"),
                velocity_average=self._average("velocity_unmasked"),
                velocity_average_masked=self._average("velocity"),
                frms_average=self._average("frms"),
                frms_background_average=self._average("frms_background"),
                delta_frms_average=self._average("delta_frms"),
                section_mask=self.section_mask,
            ),
            artery=self._vessel_signals("artery"),
            vein=self._vessel_signals("vein"),
            vessel_frms_background=self.signals["vessel_frms_background"],
            provenance=provenance,
        )

    def _average(self, name: str) -> np.ndarray:
        divisor = np.float64(max(self.frame_count, 1))
        return (self.averages[name] / divisor).astype(np.float32)

    def _vessel_signals(self, vessel: str) -> VesselVelocitySignals:
        return VesselVelocitySignals(
            velocity=self.signals[f"{vessel}_velocity"],
            frms=self.signals[f"{vessel}_frms"],
            frms_background=self.signals[f"{vessel}_frms_background"],
            delta_frms=self.signals[f"{vessel}_delta_frms"],
        )


def resolve_velocity_estimator_inputs(
    image_maps: ImageMaps,
    *,
    velocity_estimation_method: str,
    band_ratio_frequency_scale_hz: float = (
        DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ
    ),
) -> VelocityEstimatorInputs:
    """Resolve and validate the method-specific estimator input model."""

    method = str(velocity_estimation_method)
    match method:
        case "doppler_moments":
            first_name, first_volume = "moment0", image_maps.moment0
            second_name, second_volume = "moment2", image_maps.moment2
            missing_prefix = ""
        case "frequency_bands":
            first_name, first_volume = FREQUENCY_BAND_LF_PATH, image_maps.band_lf
            second_name, second_volume = FREQUENCY_BAND_HF_PATH, image_maps.band_hf
            missing_prefix = "HoloDoppler datasets "
        case _:
            raise ValueError(
                "velocity_estimation_method must be 'doppler_moments' or "
                f"'frequency_bands', got {velocity_estimation_method!r}."
            )

    missing = [
        name
        for name, volume in (
            (first_name, first_volume),
            (second_name, second_volume),
        )
        if volume is None
    ]
    if missing:
        raise ValueError(
            f"velocity_estimation_method={method!r} requires "
            f"{missing_prefix}{', '.join(missing)}."
        )
    _validate_matching_volumes(method, first_volume, second_volume)
    if method == FREQUENCY_BANDS_METHOD:
        frequency_scale_hz = validate_band_ratio_frequency_scale_hz(
            band_ratio_frequency_scale_hz
        )
        inputs: VelocityEstimatorInputs = FrequencyBandsEstimatorInputs(
            method=method,
            first_volume=first_volume,
            second_volume=second_volume,
            frequency_scale_hz=frequency_scale_hz,
        )
        Logger.log(
            "Velocity estimator converts the HoloDoppler band ratio "
            f"{FREQUENCY_BAND_HF_PATH} / {FREQUENCY_BAND_LF_PATH} to RMS "
            f"frequency using {frequency_scale_hz:g} Hz per ratio unit."
        )
    else:
        inputs = DopplerMomentsEstimatorInputs(
            method=method,
            first_volume=first_volume,
            second_volume=second_volume,
        )
        Logger.log("Velocity estimator uses raw HD moments.")
    return inputs


@singledispatch
def iter_velocity_estimator_chunks(
    inputs: VelocityEstimatorInputs,
    *,
    vessel_mask: np.ndarray,
    neighborhood_mask: np.ndarray,
) -> Iterator[VelocityEstimatorChunk]:
    """Yield normalized slice calculations for a resolved estimator type."""

    raise TypeError(f"Unsupported velocity estimator inputs: {type(inputs).__name__}.")


@singledispatch
def _build_velocity_provenance(
    inputs: VelocityEstimatorInputs,
    *,
    diagnostic_counts: Mapping[str, int],
    laser_wavelength_m: float,
    numerical_aperture: float,
) -> dict[str, object]:
    """Build complete provenance for a resolved estimator type."""

    raise TypeError(f"Unsupported velocity estimator inputs: {type(inputs).__name__}.")


@_build_velocity_provenance.register
def _build_doppler_moments_provenance(
    inputs: DopplerMomentsEstimatorInputs,
    *,
    diagnostic_counts: Mapping[str, int],
    laser_wavelength_m: float,
    numerical_aperture: float,
) -> dict[str, object]:
    del diagnostic_counts
    return physical_velocity_provenance(
        velocity_estimation_method=inputs.method,
        laser_wavelength_m=laser_wavelength_m,
        numerical_aperture=numerical_aperture,
    )


@_build_velocity_provenance.register
def _build_frequency_bands_provenance(
    inputs: FrequencyBandsEstimatorInputs,
    *,
    diagnostic_counts: Mapping[str, int],
    laser_wavelength_m: float,
    numerical_aperture: float,
) -> dict[str, object]:
    quality_counts = _empty_band_quality_counts()
    quality_counts.update(diagnostic_counts)
    Logger.log(
        "Frequency-band LF quality counts: "
        f"zero={quality_counts['band_lf_zero_sample_count']}, "
        f"near_zero={quality_counts['band_lf_near_zero_sample_count']}."
    )
    return {
        **physical_velocity_provenance(
            velocity_estimation_method=inputs.method,
            band_ratio_frequency_scale_hz=inputs.frequency_scale_hz,
            laser_wavelength_m=laser_wavelength_m,
            numerical_aperture=numerical_aperture,
        ),
        "band_lf_source_path": FREQUENCY_BAND_LF_PATH,
        "band_hf_source_path": FREQUENCY_BAND_HF_PATH,
        **quality_counts,
    }


@iter_velocity_estimator_chunks.register
def _iter_doppler_moment_chunks(
    inputs: DopplerMomentsEstimatorInputs,
    *,
    vessel_mask: np.ndarray,
    neighborhood_mask: np.ndarray,
) -> Iterator[VelocityEstimatorChunk]:
    del vessel_mask, neighborhood_mask
    for chunk_index, chunk_count, frame_slice in _estimator_chunk_slices(inputs):
        background_image = _read_volume_chunk(inputs.first_volume, frame_slice)
        moment2_chunk = _read_volume_chunk(inputs.second_volume, frame_slice)
        mean_m0 = np.mean(
            background_image,
            axis=(-1, -2),
            keepdims=True,
            dtype=np.float32,
        )
        rms_frequency = np.sqrt(
            np.divide(
                moment2_chunk,
                mean_m0,
                out=np.zeros_like(moment2_chunk, dtype=np.float32),
                where=mean_m0 != 0,
            )
        ).astype(np.float32, copy=False)
        yield VelocityEstimatorChunk(
            chunk_index=chunk_index,
            chunk_count=chunk_count,
            frame_slice=frame_slice,
            background_image=background_image,
            rms_frequency=rms_frequency,
        )


@iter_velocity_estimator_chunks.register
def _iter_frequency_band_chunks(
    inputs: FrequencyBandsEstimatorInputs,
    *,
    vessel_mask: np.ndarray,
    neighborhood_mask: np.ndarray,
) -> Iterator[VelocityEstimatorChunk]:
    for chunk_index, chunk_count, frame_slice in _estimator_chunk_slices(inputs):
        low_frequency = _read_validated_band_chunk(
            inputs.first_volume,
            frame_slice,
            FREQUENCY_BAND_LF_PATH,
        )
        high_frequency = _read_validated_band_chunk(
            inputs.second_volume,
            frame_slice,
            FREQUENCY_BAND_HF_PATH,
        )
        yield VelocityEstimatorChunk(
            chunk_index=chunk_index,
            chunk_count=chunk_count,
            frame_slice=frame_slice,
            background_image=low_frequency,
            rms_frequency=_band_ratio_to_frequency(
                _safe_band_ratio(
                    high_frequency,
                    low_frequency,
                    frame_slice=frame_slice,
                ),
                frequency_scale_hz=inputs.frequency_scale_hz,
                frame_slice=frame_slice,
            ),
            diagnostic_counts=_band_quality_counts_for_chunk(
                low_frequency,
                vessel_mask=vessel_mask,
                neighborhood_mask=neighborhood_mask,
            ),
        )


def _estimator_chunk_slices(
    inputs: VelocityEstimatorInputs,
) -> Iterator[tuple[int, int, slice]]:
    frame_count = int(inputs.first_volume.shape[0])
    chunk_count = max(
        1,
        (frame_count + SCRATCH_FRAME_CHUNK_SIZE - 1) // SCRATCH_FRAME_CHUNK_SIZE,
    )
    for chunk_index, start in enumerate(
        range(0, frame_count, SCRATCH_FRAME_CHUNK_SIZE),
        start=1,
    ):
        yield (
            chunk_index,
            chunk_count,
            slice(start, min(start + SCRATCH_FRAME_CHUNK_SIZE, frame_count)),
        )


def _doppler_frequency_to_velocity_mm_s(
    frequency_hz,
    laser_wavelength_m: float = DEFAULT_LASER_WAVELENGTH_METERS,
    numerical_aperture: float = DEFAULT_NUMERICAL_APERTURE,
) -> np.ndarray:
    """Apply v = 2 * wavelength * Doppler frequency / NA and return mm/s."""
    frequency_hz = np.asarray(frequency_hz, dtype=np.float32)
    return (
        np.float32(2e3)
        * laser_wavelength_m
        * frequency_hz
        / numerical_aperture
    ).astype(np.float32, copy=False)


def estimate_retinal_velocity(
    *,
    image_maps: ImageMaps,
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
    retain_velocity_video: bool = True,
) -> RetinalVelocityData:
    """Estimate velocity in chunks, retaining its video only when requested."""

    estimator_inputs = resolve_velocity_estimator_inputs(
        image_maps,
        velocity_estimation_method=velocity_estimation_method,
        band_ratio_frequency_scale_hz=band_ratio_frequency_scale_hz,
    )
    first_volume = estimator_inputs.first_volume
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

    Logger.log(
        f"Velocity estimator uses {SCRATCH_FRAME_CHUNK_SIZE}-frame batched "
        "inpainting and summary-only frequency intermediates."
    )

    disk, inpaint = _skimage_dependencies()
    inpaint_mask = _dilated_mask(background, disk(int(local_background_dist)))
    section_mask = annulus_mask(
        (height, width),
        optic_disc_center,
        section_inner_radius_frac,
        section_outer_radius_frac,
    )
    accumulator = _VelocityResultAccumulator.create(
        frame_count=frame_count,
        spatial_shape=(height, width),
        section_mask=section_mask,
        artery_mask=artery,
        vein_mask=vein,
        retain_velocity_video=retain_velocity_video,
    )
    diagnostic_counts: dict[str, int] = {}

    estimation_started = perf_counter()
    for estimator_chunk in iter_velocity_estimator_chunks(
        estimator_inputs,
        vessel_mask=(artery | vein),
        neighborhood_mask=~inpaint_mask,
    ):
        frame_slice = estimator_chunk.frame_slice
        f_rms = estimator_chunk.rms_frequency
        for name, count in estimator_chunk.diagnostic_counts.items():
            diagnostic_counts[name] = diagnostic_counts.get(name, 0) + count
        f_rms_background = _inpaint_frame_batch(
            f_rms,
            inpaint_mask,
            inpaint,
        )
        delta = _signed_rms_difference(f_rms, f_rms_background)
        velocity = _doppler_frequency_to_velocity_mm_s(
            delta,
            laser_wavelength_m=laser_wavelength,
            numerical_aperture=numerical_aperture,
        )
        velocity_unmasked = _doppler_frequency_to_velocity_mm_s(
            f_rms,
            laser_wavelength_m=laser_wavelength,
            numerical_aperture=numerical_aperture,
        )

        accumulator.update(
            estimator_chunk,
            rms_frequency_background=f_rms_background,
            delta_rms_frequency=delta,
            velocity=velocity,
            velocity_unmasked=velocity_unmasked,
        )
        if (
            estimator_chunk.chunk_index == estimator_chunk.chunk_count
            or estimator_chunk.chunk_index % 10 == 0
        ):
            Logger.log(
                "Velocity estimation completed chunk "
                f"{estimator_chunk.chunk_index}/{estimator_chunk.chunk_count} "
                f"({frame_slice.stop}/{frame_count} frames)."
            )

    Logger.log(
        f"Completed chunked velocity estimation in {perf_counter() - estimation_started:.1f}s."
    )

    provenance = _build_velocity_provenance(
        estimator_inputs,
        diagnostic_counts=diagnostic_counts,
        laser_wavelength_m=laser_wavelength,
        numerical_aperture=numerical_aperture,
    )
    return accumulator.build_result(provenance=provenance)


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


def _empty_band_quality_counts() -> dict[str, int | float]:
    return {
        "band_lf_low_relative_threshold": BAND_LF_LOW_RELATIVE_THRESHOLD,
        "band_lf_zero_sample_count": 0,
        "band_lf_near_zero_sample_count": 0,
        "band_lf_vessel_zero_sample_count": 0,
        "band_lf_vessel_near_zero_sample_count": 0,
        "band_lf_neighborhood_zero_sample_count": 0,
        "band_lf_neighborhood_near_zero_sample_count": 0,
    }


def _update_band_quality_counts(
    counts: dict[str, int | float],
    low_frequency: np.ndarray,
    *,
    vessel_mask: np.ndarray,
    neighborhood_mask: np.ndarray,
) -> None:
    """Accumulate zero and frame-relative near-zero LF quality counts."""

    low = np.asarray(low_frequency, dtype=np.float32)
    zero = low == 0.0
    frame_max = np.max(low, axis=(-1, -2), keepdims=True)
    near_zero = (
        (low > 0.0)
        & (low <= frame_max * np.float32(BAND_LF_LOW_RELATIVE_THRESHOLD))
    )
    vessel = np.asarray(vessel_mask, dtype=bool)[None, :, :]
    neighborhood = np.asarray(neighborhood_mask, dtype=bool)[None, :, :]
    counts["band_lf_zero_sample_count"] += int(np.count_nonzero(zero))
    counts["band_lf_near_zero_sample_count"] += int(np.count_nonzero(near_zero))
    counts["band_lf_vessel_zero_sample_count"] += int(
        np.count_nonzero(zero & vessel)
    )
    counts["band_lf_vessel_near_zero_sample_count"] += int(
        np.count_nonzero(near_zero & vessel)
    )
    counts["band_lf_neighborhood_zero_sample_count"] += int(
        np.count_nonzero(zero & neighborhood)
    )
    counts["band_lf_neighborhood_near_zero_sample_count"] += int(
        np.count_nonzero(near_zero & neighborhood)
    )


def _band_quality_counts_for_chunk(
    low_frequency: np.ndarray,
    *,
    vessel_mask: np.ndarray,
    neighborhood_mask: np.ndarray,
) -> dict[str, int]:
    counts = _empty_band_quality_counts()
    _update_band_quality_counts(
        counts,
        low_frequency,
        vessel_mask=vessel_mask,
        neighborhood_mask=neighborhood_mask,
    )
    return {
        name: int(value)
        for name, value in counts.items()
        if name.endswith("_count")
    }


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
