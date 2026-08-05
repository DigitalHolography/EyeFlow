"""Copied from DopplerView pipeline/steps/vessel_velocity_estimator.py."""

from __future__ import annotations

import os
from functools import partial
from time import perf_counter

import numpy as np
from scipy import ndimage as ndi

from calculations.blood_flow_velocity.cross_section.segment_geometry import annulus_mask
from runtime_limits import cap_parallel_jobs
from utils.logger import Logger

from .flat_field import (
    corrected_flat_field_chunk,
    fit_flat_field_parameters,
)


SCRATCH_FRAME_CHUNK_SIZE = 8
SECTION_INNER_RADIUS_FRAC = 0.10
SECTION_OUTER_RADIUS_FRAC = 0.35
DOPPLER_WAVELENGTH_METERS = 8.52e-7
DOPPLER_ANGLE_RADIANS = 0.25


def _velocity_scale_mm_per_s_per_hz() -> np.float32:
    return np.float32(
        1000.0 * 2.0 * DOPPLER_WAVELENGTH_METERS / DOPPLER_ANGLE_RADIANS
    )


class VesselVelocityEstimatorStep:

    def _relevant_config(self, ctx):
        return {
            "LocalBackgroundDist": ctx.dv_config_value(
                "VelocityEstimation",
                "LocalBackgroundDist",
            ),
        }

    def run(self, ctx):

        # ---- Requires ----
        moment0 = np.asarray(ctx.require("moment0"), dtype=np.float32)
        moment2 = np.asarray(ctx.require("moment2"), dtype=np.float32)

        artery_mask = ctx.require("retinal_artery_mask")
        vein_mask = ctx.require("retinal_vein_mask")
        vessel_mask = artery_mask | vein_mask

        # Compute fRMS
        mean_m0 = np.mean(moment0, axis=(-1, -2), keepdims=True, dtype=np.float32)
        fRMS = np.sqrt(moment2 / mean_m0).astype(np.float32, copy=False)

        # Inpaint fRMS to estimate background
        local_background_dist = ctx.dv_config_value(
            "VelocityEstimation",
            "LocalBackgroundDist",
        )
        disk, inpaint = _skimage_dependencies()
        mask = _dilated_mask(vessel_mask, disk(local_background_dist))

        n_jobs = cap_parallel_jobs(_cpu_count())

        print(f"    - Inpainting fRMS with {n_jobs} parallel jobs")

        def _inpaint_frame(frame, mask):
            return inpaint.inpaint_biharmonic(frame, mask).astype(np.float32)

        fRMSbkg = _run_in_parallel(
            partial(_inpaint_frame, mask=mask), fRMS, n_jobs=n_jobs, chunking=False
        )
        ctx.set("fRMS", fRMS)
        ctx.set("fRMS_bkg", fRMSbkg)

        # fRMSbkg = np.stack(np.array([inpaint.inpaint_biharmonic(frame, mask) for frame in fRMS]), axis=0)

        # Velocity estimation
        A = (fRMS**2 - fRMSbkg**2).astype(np.float32, copy=False)
        deltafRMS = (np.sign(A) * np.sqrt(np.abs(A))).astype(np.float32, copy=False)
        ctx.set("deltafRMS", deltafRMS)

        velocity_scale = _velocity_scale_mm_per_s_per_hz()
        velocity_map = (velocity_scale * deltafRMS).astype(
            np.float32,
            copy=False,
        )  # mm/s

        ctx.set("velocity_map", velocity_map)

        # num_bins = 256  # for 8-bit grayscale
        # hist_matrix = np.zeros((velocity_map.shape[2], num_bins))
        # v_range = (velocity_map.min(),velocity_map.max())

        # for i in range(velocity_map.shape[2]):
        #     masked_pixels = velocity_map[:,:,i][mask]  # select only pixels under mask
        #     hist, _ = np.histogram(masked_pixels, bins=num_bins, range=v_range)
        #     hist_matrix[i,:] = hist

        # ctx.set("hist_matrix", hist_matrix)
        ctx.set("velocity_map_avg", np.mean(velocity_map, axis=0, dtype=np.float32))
        ctx.set("fRMS_avg", np.mean(fRMS, axis=0, dtype=np.float32))
        ctx.set("fRMS_bkg_avg", np.mean(fRMSbkg, axis=0, dtype=np.float32))

        sz = velocity_map.shape
        optic_disc_center = getattr(ctx, "cache", {}).get("optic_disc_center")
        section_mask = annulus_mask(
            (sz[-2], sz[-1]),
            optic_disc_center,
            SECTION_INNER_RADIUS_FRAC,
            SECTION_OUTER_RADIUS_FRAC,
        )
        ctx.set("velocity_section_mask", section_mask)

        artery_sig = _masked_signal(velocity_map, section_mask & artery_mask)

        vein_sig = _masked_signal(velocity_map, section_mask & vein_mask)
        ctx.set("retinal_vessel_velocity", velocity_map)
        ctx.set("retinal_artery_velocity_signal", artery_sig)
        ctx.set("retinal_vein_velocity_signal", vein_sig)


def run_chunked_velocity_estimator(
    *,
    moment0,
    moment2,
    moment0_flat_field_source: str = "dopplerview_recomputed_from_raw",
    moment2_flat_field_source: str = "dopplerview_recomputed_from_raw",
    artery_mask,
    vein_mask,
    optic_disc_center=None,
    optic_disc_width=None,
    optic_disc_height=None,
    section_inner_radius_frac: float = SECTION_INNER_RADIUS_FRAC,
    section_outer_radius_frac: float = SECTION_OUTER_RADIUS_FRAC,
    local_background_dist: int,
    flat_field_gaussian_ratio: float,
    flat_field_border: float,
    scratch_h5,
) -> dict[str, object]:
    """Estimate velocity into scratch datasets without materializing full videos."""

    if tuple(moment0.shape) != tuple(moment2.shape) or len(moment0.shape) != 3:
        raise ValueError(
            "moment0 and moment2 must be matching 3-D datasets, got "
            f"{moment0.shape} and {moment2.shape}."
        )
    frame_count, height, width = (int(size) for size in moment0.shape)
    artery = np.asarray(artery_mask, dtype=bool)
    vein = np.asarray(vein_mask, dtype=bool)
    if artery.shape != (height, width) or vein.shape != (height, width):
        raise ValueError("Velocity masks must match the HD moment spatial shape.")

    gaussian_width = float(flat_field_gaussian_ratio) * float(width)
    preparation_started = perf_counter()
    moment0_params = _flat_field_parameters_if_needed(
        moment0, gaussian_width, flat_field_border, "moment0", moment0_flat_field_source
    )
    moment2_params = _flat_field_parameters_if_needed(
        moment2, gaussian_width, flat_field_border, "moment2", moment2_flat_field_source
    )
    Logger.log(
        "Completed DopplerView flat-field preparation "
        f"in {perf_counter() - preparation_started:.1f}s "
        f"(moment0={moment0_flat_field_source}, moment2={moment2_flat_field_source})."
    )

    group = scratch_h5.require_group("waveform")
    datasets = {
        name: group.create_dataset(
            name,
            shape=(frame_count, height, width),
            dtype=np.float32,
            chunks=(
                min(SCRATCH_FRAME_CHUNK_SIZE, frame_count),
                min(128, height),
                min(128, width),
            ),
            compression="lzf",
        )
        for name in ("fRMS", "fRMS_bkg", "deltafRMS", "velocity")
    }
    vessel_mask = artery | vein
    disk, inpaint = _skimage_dependencies()
    inpaint_mask = _dilated_mask(vessel_mask, disk(int(local_background_dist)))
    section_mask = annulus_mask(
        (height, width),
        optic_disc_center,
        section_inner_radius_frac,
        section_outer_radius_frac,
        optic_disc_width=optic_disc_width,
        optic_disc_height=optic_disc_height,
        optic_disc_boundary_radius_frac=section_inner_radius_frac,
    )
    artery_section = section_mask & artery
    vein_section = section_mask & vein

    averages = {
        name: np.zeros((height, width), dtype=np.float64)
        for name in ("moment0", "velocity", "fRMS", "fRMS_bkg")
    }
    artery_signal = np.full(frame_count, np.nan, dtype=np.float32)
    vein_signal = np.full(frame_count, np.nan, dtype=np.float32)
    velocity_scale = _velocity_scale_mm_per_s_per_hz()
    n_jobs = min(cap_parallel_jobs(_cpu_count()), SCRATCH_FRAME_CHUNK_SIZE)

    estimation_started = perf_counter()
    chunk_count = max(1, (frame_count + SCRATCH_FRAME_CHUNK_SIZE - 1) // SCRATCH_FRAME_CHUNK_SIZE)
    for chunk_index, start in enumerate(range(0, frame_count, SCRATCH_FRAME_CHUNK_SIZE), start=1):
        stop = min(start + SCRATCH_FRAME_CHUNK_SIZE, frame_count)
        frame_slice = slice(start, stop)
        m0 = _read_moment_chunk(moment0, frame_slice, moment0_params)
        m2 = _read_moment_chunk(moment2, frame_slice, moment2_params)
        mean_m0 = np.mean(m0, axis=(-1, -2), keepdims=True, dtype=np.float32)
        f_rms = np.sqrt(
            np.divide(
                m2,
                mean_m0,
                out=np.zeros_like(m2, dtype=np.float32),
                where=mean_m0 != 0,
            )
        ).astype(np.float32, copy=False)
        f_background = _run_in_parallel(
            lambda frame: inpaint.inpaint_biharmonic(
                frame,
                inpaint_mask,
            ).astype(np.float32),
            f_rms,
            n_jobs=n_jobs,
            chunking=False,
        )
        difference = (f_rms**2 - f_background**2).astype(np.float32, copy=False)
        delta = (np.sign(difference) * np.sqrt(np.abs(difference))).astype(
            np.float32,
            copy=False,
        )
        velocity = (velocity_scale * delta).astype(np.float32, copy=False)

        datasets["fRMS"][frame_slice] = f_rms
        datasets["fRMS_bkg"][frame_slice] = f_background
        datasets["deltafRMS"][frame_slice] = delta
        datasets["velocity"][frame_slice] = velocity
        averages["moment0"] += np.sum(m0, axis=0, dtype=np.float64)
        averages["velocity"] += np.sum(velocity, axis=0, dtype=np.float64)
        averages["fRMS"] += np.sum(f_rms, axis=0, dtype=np.float64)
        averages["fRMS_bkg"] += np.sum(f_background, axis=0, dtype=np.float64)
        artery_signal[frame_slice] = _masked_signal(velocity, artery_section)
        vein_signal[frame_slice] = _masked_signal(velocity, vein_section)
        if chunk_index == chunk_count or chunk_index % 10 == 0:
            Logger.log(
                f"Velocity estimation completed chunk {chunk_index}/{chunk_count} "
                f"({stop}/{frame_count} frames)."
            )

    Logger.log(
        f"Completed chunked velocity estimation in {perf_counter() - estimation_started:.1f}s."
    )

    divisor = np.float64(max(frame_count, 1))
    return {
        "fRMS": datasets["fRMS"],
        "fRMS_bkg": datasets["fRMS_bkg"],
        "deltafRMS": datasets["deltafRMS"],
        "velocity_map": datasets["velocity"],
        "retinal_vessel_velocity": datasets["velocity"],
        "moment0_avg": (averages["moment0"] / divisor).astype(np.float32),
        "velocity_map_avg": (averages["velocity"] / divisor).astype(np.float32),
        "fRMS_avg": (averages["fRMS"] / divisor).astype(np.float32),
        "fRMS_bkg_avg": (averages["fRMS_bkg"] / divisor).astype(np.float32),
        "velocity_section_mask": section_mask,
        "velocity_section_geometry": (
            "optic_disc_relative"
            if _has_optic_disc_geometry(optic_disc_width, optic_disc_height)
            else "frame_relative_fallback"
        ),
        "retinal_artery_velocity_signal": artery_signal,
        "retinal_vein_velocity_signal": vein_signal,
        "moment0_flat_field_source": moment0_flat_field_source,
        "moment2_flat_field_source": moment2_flat_field_source,
    }


def _has_optic_disc_geometry(width, height) -> bool:
    for value in (width, height):
        if value is None:
            return False
        array = np.asarray(value, dtype=np.float32).reshape(-1)
        if array.size == 0 or not np.isfinite(array[0]) or array[0] <= 0:
            return False
    return True


def _flat_field_parameters(
    volume,
    gaussian_width: float,
    border_amount: float,
    name: str,
):
    Logger.log(f"Calculating DopplerView flat field from raw HD {name}.")
    return fit_flat_field_parameters(
        volume,
        gaussian_width=gaussian_width,
        border_amount=border_amount,
        frame_chunk_size=SCRATCH_FRAME_CHUNK_SIZE,
    )


def _flat_field_parameters_if_needed(
    volume,
    gaussian_width: float,
    border_amount: float,
    name: str,
    source: str,
):
    if source == "holodoppler_precomputed_flat_field":
        Logger.log(f"Reusing precomputed HD {name} flat field.")
        return None
    return _flat_field_parameters(volume, gaussian_width, border_amount, name)


def _read_moment_chunk(volume, frame_slice: slice, parameters) -> np.ndarray:
    if parameters is None:
        return np.asarray(volume[frame_slice], dtype=np.float32)
    return corrected_flat_field_chunk(volume, frame_slice, parameters)


def _skimage_dependencies():
    try:
        from skimage.morphology import disk
        from skimage.restoration import inpaint
    except ModuleNotFoundError as exc:
        raise ImportError(
            "DopplerView velocity estimation requires scikit-image."
        ) from exc
    return disk, inpaint


def _dilated_mask(vessel_mask: np.ndarray, footprint: np.ndarray) -> np.ndarray:
    mask = np.asarray(vessel_mask, dtype=bool)
    if mask.ndim != 2:
        raise ValueError(f"vessel_mask must be 2-D for dilation, got {mask.shape}.")
    return ndi.binary_dilation(mask, structure=np.asarray(footprint, dtype=bool))


def _cpu_count() -> int:
    return os.cpu_count() or 1


def _masked_signal(velocity_map: np.ndarray, mask: np.ndarray) -> np.ndarray:
    selected = velocity_map[:, np.asarray(mask, dtype=bool)]
    if not np.any(np.isfinite(selected)):
        return np.full((velocity_map.shape[0],), np.nan, dtype=np.float32)
    return np.nanmean(selected, axis=1, dtype=np.float32).astype(
        np.float32,
        copy=False,
    )


def _run_in_parallel(func, iterable, n_jobs=-1, chunking=False):
    try:
        import joblib
    except ModuleNotFoundError:
        return np.stack([func(item) for item in iterable], axis=0).astype(
            np.float32,
            copy=False,
        )
    if n_jobs == -1:
        n_jobs = _cpu_count()
    n_jobs = cap_parallel_jobs(n_jobs)
    results = joblib.Parallel(n_jobs=n_jobs, backend="threading")(
        joblib.delayed(func)(item) for item in iterable
    )
    return np.stack(results, axis=0).astype(np.float32, copy=False)
