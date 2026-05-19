"""Copied from DopplerView pipeline/steps/vessel_velocity_estimator.py."""

from __future__ import annotations

from functools import partial

import numpy as np
from scipy import ndimage as ndi

from calculations.math import butter_lowpass_filtfilt
from runtime_limits import cap_parallel_jobs

from ._masks import elliptical_mask


class VesselVelocityEstimatorStep:

    def _relevant_config(self, ctx):
        return {
            "LocalBackgroundDist": ctx.dv_config_value(
                "VelocityEstimation",
                "LocalBackgroundDist",
            ),
            "FilterSignals": ctx.dv_config_value(
                "PulseAnalysis",
                "FilterSignals",
                True,
            ),
            "LowpassFreqHz": ctx.dv_config_value(
                "PulseAnalysis",
                "LowpassFreqHz",
                15.0,
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

        # fRMSbkg = np.stack(np.array([inpaint.inpaint_biharmonic(frame, mask) for frame in fRMS]), axis=0)

        # Velocity estimation
        A = (fRMS**2 - fRMSbkg**2).astype(np.float32, copy=False)
        deltafRMS = (np.sign(A) * np.sqrt(np.abs(A))).astype(np.float32, copy=False)

        velocity_scale = np.float32(2 * 852e-9 / np.sin(0.25) * 1e6)
        velocity_map = (velocity_scale * deltafRMS).astype(np.float32)  # mm/s

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

        section_mask = elliptical_mask(sz[-2], sz[-1], 0.5) & (~(elliptical_mask(sz[-2], sz[-1], 0.2)))

        artery_sig = _masked_signal(velocity_map, section_mask & artery_mask)

        vein_sig = _masked_signal(velocity_map, section_mask & vein_mask)
        if ctx.dv_config_value("PulseAnalysis", "FilterSignals", True):
            dt_seconds = _dt_seconds(ctx)
            lowpass_freq_hz = float(
                ctx.dv_config_value("PulseAnalysis", "LowpassFreqHz", 15.0)
            )
            artery_sig = butter_lowpass_filtfilt(
                artery_sig,
                dt_seconds=dt_seconds,
                lowpass_freq_hz=lowpass_freq_hz,
                order=4,
            )
            vein_sig = butter_lowpass_filtfilt(
                vein_sig,
                dt_seconds=dt_seconds,
                lowpass_freq_hz=lowpass_freq_hz,
                order=4,
            )

        ctx.set("retinal_vessel_velocity", velocity_map)
        ctx.set("retinal_artery_velocity_signal", artery_sig)
        ctx.set("retinal_vein_velocity_signal", vein_sig)


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
    try:
        import joblib
    except ModuleNotFoundError:
        return 1
    return joblib.cpu_count()


def _masked_signal(velocity_map: np.ndarray, mask: np.ndarray) -> np.ndarray:
    count = np.count_nonzero(mask)
    if count == 0:
        return np.full((velocity_map.shape[0],), np.nan, dtype=np.float32)
    total = np.sum(velocity_map * mask, axis=(-2, -1), dtype=np.float32)
    return (total / np.float32(count)).astype(np.float32, copy=False)


def _dt_seconds(ctx) -> float:
    fs = float(ctx.hd_config_value("sampling_freq"))
    stride = float(ctx.hd_config_value("batch_stride"))
    if fs <= 0:
        raise ValueError("sampling_freq must be positive for velocity filtering.")
    return stride / fs


def _run_in_parallel(func, iterable, n_jobs=-1, chunking=False):
    try:
        import joblib
    except ModuleNotFoundError:
        return np.stack([func(item) for item in iterable], axis=0).astype(np.float32)
    if n_jobs == -1:
        n_jobs = joblib.cpu_count()
    n_jobs = cap_parallel_jobs(n_jobs)
    results = joblib.Parallel(n_jobs=n_jobs, backend="threading")(
        joblib.delayed(func)(item) for item in iterable
    )
    return np.stack(results, axis=0).astype(np.float32)
