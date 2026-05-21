"""Small MATLAB-legacy pulse-analysis proof of concept."""

from __future__ import annotations

import numpy as np

from calculations.blood_flow_velocity.context_builders.signal import find_systole_index
from calculations.blood_flow_velocity.signal_analysis.signal.per_beat_signal import (
    per_beat_signal_analysis,
)
from calculations.math import butter_lowpass_filtfilt
from input_output.schema import DopplerViewSource, HD_MOMENT0_PATH, HD_MOMENT2_PATH
from pipeline_engine.imports import (
    ProcessResult,
    pipeline,
    resolve_holodoppler_timing,
    with_attrs,
)

MOMENT2_MATLAB_SCALE = np.float32(1e-6)
VELOCITY_SCALE_MM_S = np.float32(1000 * 1000 * 2 * 8.52e-7 / 0.25)
LOWPASS_FREQ_HZ = np.float32(15.0)
ANNULUS_INNER_RADIUS = np.float32(0.10)
ANNULUS_OUTER_RADIUS = np.float32(0.35)
BAND_LIMITED_SIGNAL_HARMONIC_COUNT = 13
DV_ARTERY_MASK_PATH = "segmentation/Retina/artery_mask"
DV_VEIN_MASK_PATH = "segmentation/Retina/vein_mask"


@pipeline(
    name="matlab_pulse_poc",
    description=(
        "POC: compute MATLAB-like retinal velocity and heartbeat metrics from HD "
        "moments plus DopplerView artery/vein masks."
    ),
    requires=["numpy", "h5py", "scipy"],
    dag_produces=["matlab_pulse_poc"],
    input_slot="both",
)
def run(ctx) -> ProcessResult:
    ctx.require_inputs("hd", "dv")

    timing = resolve_holodoppler_timing(ctx)
    dv_source = DopplerViewSource.from_context(ctx)
    artery_mask = dv_source.retinal_artery_mask()
    vein_mask = dv_source.retinal_vein_mask()
    vessel_mask = artery_mask | vein_mask
    section_mask = _annulus_mask(artery_mask.shape)

    artery_section = artery_mask & section_mask
    vein_section = vein_mask & section_mask
    background_section = section_mask & ~vessel_mask
    _validate_masks(artery_section, vein_section, background_section)

    moment0 = ctx.sources.hd.dataset(HD_MOMENT0_PATH)
    moment2 = ctx.sources.hd.dataset(HD_MOMENT2_PATH)
    artery_signal, vein_signal = _velocity_signals(
        moment0=moment0,
        moment2=moment2,
        artery_mask=artery_section,
        vein_mask=vein_section,
        background_mask=background_section,
    )

    artery_filtered = butter_lowpass_filtfilt(
        artery_signal,
        dt_seconds=np.float32(timing.dt_seconds),
        lowpass_freq_hz=LOWPASS_FREQ_HZ,
    )
    vein_filtered = butter_lowpass_filtfilt(
        vein_signal,
        dt_seconds=np.float32(timing.dt_seconds),
        lowpass_freq_hz=LOWPASS_FREQ_HZ,
    )

    detection = find_systole_index(
        artery_signal,
        dt=np.float32(timing.dt_seconds),
        lowpass_freq_hz=LOWPASS_FREQ_HZ,
    )
    systole_zero_based = detection.systole_indexes.astype(np.int32, copy=False)
    beat_period_idx = np.diff(systole_zero_based).astype(np.int32, copy=False)
    beat_period_seconds = (
        beat_period_idx.astype(np.float32) * np.float32(timing.dt_seconds)
    )
    artery_per_beat = _per_beat_waveforms(artery_filtered, systole_zero_based)
    vein_per_beat = _per_beat_waveforms(vein_filtered, systole_zero_based)
    artery_vti_per_beat = _vti_per_beat(
        artery_per_beat.velocity_signal_per_beat,
        np.float32(timing.dt_seconds),
    )
    vein_vti_per_beat = _vti_per_beat(
        vein_per_beat.velocity_signal_per_beat,
        np.float32(timing.dt_seconds),
    )
    artery_stats = _waveform_stats(artery_filtered, systole_zero_based)
    vein_stats = _waveform_stats(vein_filtered, systole_zero_based)

    metrics = _legacy_global_metrics(
        "Artery",
        filtered_signal=artery_filtered,
        raw_signal=artery_signal,
        stats=artery_stats,
    )
    metrics.update(
        _legacy_global_metrics(
            "Vein",
            filtered_signal=vein_filtered,
            raw_signal=vein_signal,
            stats=vein_stats,
        )
    )
    metrics.update(
        _legacy_arterial_waveform_metrics(
            systole_zero_based=systole_zero_based,
            beat_period_idx=beat_period_idx,
            beat_period_seconds=beat_period_seconds,
            heart_rate_bpm=_heart_rate_bpm(beat_period_seconds),
        )
    )
    metrics.update(
        _legacy_per_beat_metrics(
            "Artery",
            per_beat=artery_per_beat,
            vti_per_beat=artery_vti_per_beat,
            beat_period_idx=beat_period_idx,
            beat_period_seconds=beat_period_seconds,
        )
    )
    metrics.update(
        _legacy_per_beat_metrics(
            "Vein",
            per_beat=vein_per_beat,
            vti_per_beat=vein_vti_per_beat,
            beat_period_idx=beat_period_idx,
            beat_period_seconds=beat_period_seconds,
        )
    )
    metrics.update(
        {
            "poc/matlab_pulse/artery_velocity_signal": with_attrs(
                artery_filtered,
                {"unit": "mm/s", "dimDesc": ["frame"]},
            ),
            "poc/matlab_pulse/vein_velocity_signal": with_attrs(
                vein_filtered,
                {"unit": "mm/s", "dimDesc": ["frame"]},
            ),
            "poc/matlab_pulse/artery_velocity_signal_raw": with_attrs(
                artery_signal,
                {"unit": "mm/s", "dimDesc": ["frame"]},
            ),
            "poc/matlab_pulse/vein_velocity_signal_raw": with_attrs(
                vein_signal,
                {"unit": "mm/s", "dimDesc": ["frame"]},
            ),
            "poc/matlab_pulse/systolic_acceleration_peak_indexes_zero_based": (
                systole_zero_based,
                {"unit": "frame", "index_base": np.int32(0), "dimDesc": ["beat"]},
            ),
            "poc/matlab_pulse/systolic_acceleration_peak_indexes_matlab": (
                systole_zero_based + np.int32(1),
                {"unit": "frame", "index_base": np.int32(1), "dimDesc": ["beat"]},
            ),
            "poc/matlab_pulse/beat_period_seconds": with_attrs(
                beat_period_seconds,
                {"unit": "s", "dimDesc": ["beat_interval"]},
            ),
            "poc/matlab_pulse/artery_vti_per_beat": with_attrs(
                artery_vti_per_beat,
                {"unit": "mm", "dimDesc": ["beat"]},
            ),
            "poc/matlab_pulse/vein_vti_per_beat": with_attrs(
                vein_vti_per_beat,
                {"unit": "mm", "dimDesc": ["beat"]},
            ),
            "poc/matlab_pulse/artery_vti_mean": with_attrs(
                _nanmean_or_nan(artery_vti_per_beat),
                {"unit": "mm"},
            ),
            "poc/matlab_pulse/vein_vti_mean": with_attrs(
                _nanmean_or_nan(vein_vti_per_beat),
                {"unit": "mm"},
            ),
            "poc/matlab_pulse/artery_vti_std": with_attrs(
                _nanstd_or_nan(artery_vti_per_beat),
                {"unit": "mm"},
            ),
            "poc/matlab_pulse/vein_vti_std": with_attrs(
                _nanstd_or_nan(vein_vti_per_beat),
                {"unit": "mm"},
            ),
            "poc/matlab_pulse/heart_rate_bpm": with_attrs(
                _heart_rate_bpm(beat_period_seconds),
                {"unit": "bpm"},
            ),
            "poc/matlab_pulse/artery_resistivity_index": with_attrs(
                artery_stats["resistivity_index"],
                {"matlab_source": "ArterialResistivityIndex.m"},
            ),
            "poc/matlab_pulse/vein_resistivity_index": with_attrs(
                vein_stats["resistivity_index"],
                {"matlab_source": "ArterialResistivityIndex.m"},
            ),
            "poc/matlab_pulse/artery_mean_velocity": with_attrs(
                artery_stats["mean"],
                {"unit": "mm/s"},
            ),
            "poc/matlab_pulse/vein_mean_velocity": with_attrs(
                vein_stats["mean"],
                {"unit": "mm/s"},
            ),
            "poc/matlab_pulse/artery_mask_pixel_count": with_attrs(
                np.int32(np.count_nonzero(artery_section)),
                {"unit": "pixels"},
            ),
            "poc/matlab_pulse/vein_mask_pixel_count": with_attrs(
                np.int32(np.count_nonzero(vein_section)),
                {"unit": "pixels"},
            ),
        }
    )

    attrs = {
        "matlab_reference_step": "BloodFlowVelocity/pulseAnalysis.m",
        "poc_background_method": (
            "mean fRMS over non-vessel pixels in fixed elliptical annulus"
        ),
        "poc_frame_transform": "transpose HD y/x frames before applying DV masks",
        "poc_moment2_scale": float(MOMENT2_MATLAB_SCALE),
        "poc_velocity_scale_mm_s": float(VELOCITY_SCALE_MM_S),
        "poc_annulus_inner_radius": float(ANNULUS_INNER_RADIUS),
        "poc_annulus_outer_radius": float(ANNULUS_OUTER_RADIUS),
        "poc_lowpass_freq_hz": float(LOWPASS_FREQ_HZ),
        "poc_band_limited_signal_harmonic_count": int(
            BAND_LIMITED_SIGNAL_HARMONIC_COUNT
        ),
        "poc_vti_method": (
            "perBeatAnalysis.m style: sum(interpolated filtered velocity beat) * dt"
        ),
        "poc_hd_moment0_path": HD_MOMENT0_PATH,
        "poc_hd_moment2_path": HD_MOMENT2_PATH,
        "poc_dv_artery_mask_path": DV_ARTERY_MASK_PATH,
        "poc_dv_vein_mask_path": DV_VEIN_MASK_PATH,
        "poc_sampling_freq": float(timing.sampling_freq),
        "poc_batch_stride": float(timing.batch_stride),
        "poc_dt_seconds": float(timing.dt_seconds),
    }
    ctx.set_var("matlab_pulse_poc", {"systole_indexes": systole_zero_based})
    return ProcessResult(metrics=metrics, attrs=attrs)


def _velocity_signals(
    *,
    moment0,
    moment2,
    artery_mask: np.ndarray,
    vein_mask: np.ndarray,
    background_mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    if moment0.shape != moment2.shape:
        raise ValueError(
            f"moment0 and moment2 shapes differ: {moment0.shape} != {moment2.shape}"
        )
    if len(moment0.shape) != 3:
        raise ValueError(f"HD moments must be 3-D, got shape {moment0.shape}.")

    n_frames = int(moment0.shape[0])
    artery_signal = np.empty(n_frames, dtype=np.float32)
    vein_signal = np.empty(n_frames, dtype=np.float32)

    for frame_index in range(n_frames):
        m0_frame = _hd_frame_for_dv_mask(moment0[frame_index], artery_mask.shape)
        m2_frame = _hd_frame_for_dv_mask(moment2[frame_index], artery_mask.shape)
        velocity_frame = _velocity_frame(m0_frame, m2_frame, background_mask)
        artery_signal[frame_index] = np.nanmean(velocity_frame[artery_mask])
        vein_signal[frame_index] = np.nanmean(velocity_frame[vein_mask])

    return artery_signal, vein_signal


def _hd_frame_for_dv_mask(frame, mask_shape: tuple[int, int]) -> np.ndarray:
    value = np.asarray(frame, dtype=np.float32).T
    if value.shape != mask_shape:
        raise ValueError(
            "Transposed HD frame shape does not match the DopplerView mask shape: "
            f"{value.shape} != {mask_shape}."
        )
    return value


def _velocity_frame(
    moment0_frame: np.ndarray,
    moment2_frame: np.ndarray,
    background_mask: np.ndarray,
) -> np.ndarray:
    mean_m0 = np.nanmean(moment0_frame, dtype=np.float32)
    if not np.isfinite(mean_m0) or mean_m0 <= 0:
        return np.full(moment0_frame.shape, np.nan, dtype=np.float32)

    moment2_scaled = moment2_frame * MOMENT2_MATLAB_SCALE
    f_rms = np.sqrt(np.maximum(moment2_scaled / mean_m0, 0.0)).astype(np.float32)
    background = np.nanmean(f_rms[background_mask], dtype=np.float32)
    delta = (f_rms * f_rms - background * background).astype(np.float32, copy=False)
    return (VELOCITY_SCALE_MM_S * np.sign(delta) * np.sqrt(np.abs(delta))).astype(
        np.float32,
        copy=False,
    )


def _annulus_mask(shape: tuple[int, int]) -> np.ndarray:
    ny, nx = shape
    y_grid, x_grid = np.ogrid[:ny, :nx]
    cy = np.float32(ny / 2.0)
    cx = np.float32(nx / 2.0)
    y = y_grid.astype(np.float32) - cy
    x = x_grid.astype(np.float32) - cx
    outer = (x / (np.float32(nx / 2.0) * ANNULUS_OUTER_RADIUS)) ** 2 + (
        y / (np.float32(ny / 2.0) * ANNULUS_OUTER_RADIUS)
    ) ** 2
    inner = (x / (np.float32(nx / 2.0) * ANNULUS_INNER_RADIUS)) ** 2 + (
        y / (np.float32(ny / 2.0) * ANNULUS_INNER_RADIUS)
    ) ** 2
    return (outer <= 1.0) & (inner > 1.0)


def _validate_masks(
    artery_mask: np.ndarray,
    vein_mask: np.ndarray,
    background_mask: np.ndarray,
) -> None:
    if not np.any(artery_mask):
        raise ValueError("DV artery mask has no pixels inside the POC annulus.")
    if not np.any(vein_mask):
        raise ValueError("DV vein mask has no pixels inside the POC annulus.")
    if not np.any(background_mask):
        raise ValueError("No background pixels are available inside the POC annulus.")


def _heart_rate_bpm(beat_period_seconds: np.ndarray) -> np.float32:
    if beat_period_seconds.size == 0:
        return np.float32(np.nan)
    return np.float32(60.0 / np.nanmean(beat_period_seconds))


def _per_beat_waveforms(signal: np.ndarray, systole_indexes: np.ndarray):
    return per_beat_signal_analysis(
        signal,
        systole_indexes,
        BAND_LIMITED_SIGNAL_HARMONIC_COUNT,
        index_base=0,
    )


def _legacy_global_metrics(
    vessel: str,
    *,
    filtered_signal: np.ndarray,
    raw_signal: np.ndarray,
    stats: dict[str, np.float32],
) -> dict[str, object]:
    root = f"{vessel}/Velocity"
    return {
        f"{root}/VelocitySignal/value": with_attrs(
            filtered_signal,
            {"unit": "mm/s", "dimDesc": ["frame"]},
        ),
        f"{root}/VelocitySignalRaw/value": with_attrs(
            raw_signal,
            {"unit": "mm/s", "dimDesc": ["frame"]},
        ),
        f"{root}/MeanSignal/value": with_attrs(stats["mean"], {"unit": "mm/s"}),
        f"{root}/MaxSignal/value": with_attrs(stats["max"], {"unit": "mm/s"}),
        f"{root}/MinSignal/value": with_attrs(stats["min"], {"unit": "mm/s"}),
        f"{root}/ResistivityIndexSignal/value": stats["resistivity_index"],
        f"{root}/PulsatilityIndexSignal/value": stats["pulsatility_index"],
        f"{root}/MaxMinRatioSignal/value": stats["max_min_ratio"],
    }


def _legacy_arterial_waveform_metrics(
    *,
    systole_zero_based: np.ndarray,
    beat_period_idx: np.ndarray,
    beat_period_seconds: np.ndarray,
    heart_rate_bpm: np.float32,
) -> dict[str, object]:
    systole_matlab = systole_zero_based + np.int32(1)
    return {
        "Artery/Velocity/SystolicAccelerationPeakIndexes/value": (
            systole_matlab,
            {"unit": "frame", "index_base": np.int32(1), "dimDesc": ["beat"]},
        ),
        "Artery/WaveformAnalysis/SystoleIndices/value": (
            systole_matlab,
            {"unit": "frame", "index_base": np.int32(1), "dimDesc": ["beat"]},
        ),
        "Artery/WaveformAnalysis/HeartBeat/value": with_attrs(
            heart_rate_bpm,
            {"unit": "bpm"},
        ),
        "Artery/VelocityPerBeat/beatPeriodIdx/value": with_attrs(
            beat_period_idx.astype(np.float32, copy=False).reshape(1, -1),
            {"unit": "frame", "dimDesc": ["beat"]},
        ),
        "Artery/VelocityPerBeat/beatPeriodSeconds/value": with_attrs(
            beat_period_seconds.reshape(1, -1),
            {"unit": "s", "dimDesc": ["beat"]},
        ),
        "Vein/VelocityPerBeat/beatPeriodIdx/value": with_attrs(
            beat_period_idx.astype(np.float32, copy=False).reshape(1, -1),
            {"unit": "frame", "dimDesc": ["beat"]},
        ),
        "Vein/VelocityPerBeat/beatPeriodSeconds/value": with_attrs(
            beat_period_seconds.reshape(1, -1),
            {"unit": "s", "dimDesc": ["beat"]},
        ),
    }


def _legacy_per_beat_metrics(
    vessel: str,
    *,
    per_beat,
    vti_per_beat: np.ndarray,
    beat_period_idx: np.ndarray,
    beat_period_seconds: np.ndarray,
) -> dict[str, object]:
    root = f"{vessel}/VelocityPerBeat"
    signal = per_beat.velocity_signal_per_beat
    fft = per_beat.velocity_signal_per_beat_fft
    band_limited = per_beat.velocity_signal_per_beat_band_limited

    return {
        f"{root}/VelocitySignalPerBeat/value": with_attrs(
            signal,
            {"unit": "mm/s", "dimDesc": ["beat", "sample"]},
        ),
        f"{root}/VelocitySignalPerBeatFFT_abs/value": with_attrs(
            np.abs(fft).astype(np.float32, copy=False),
            {"unit": "a.u.", "dimDesc": ["beat", "frequency_bin"]},
        ),
        f"{root}/VelocitySignalPerBeatFFT_arg/value": with_attrs(
            np.angle(fft).astype(np.float32, copy=False),
            {"unit": "rad", "dimDesc": ["beat", "frequency_bin"]},
        ),
        f"{root}/VelocitySignalPerBeatBandLimited/value": with_attrs(
            band_limited,
            {"unit": "mm/s", "dimDesc": ["beat", "sample"]},
        ),
        f"{root}/VmaxPerBeatBandLimited/value": with_attrs(
            np.nanmax(band_limited, axis=1).astype(np.float32, copy=False),
            {"unit": "mm/s", "dimDesc": ["beat"]},
        ),
        f"{root}/VminPerBeatBandLimited/value": with_attrs(
            np.nanmin(band_limited, axis=1).astype(np.float32, copy=False),
            {"unit": "mm/s", "dimDesc": ["beat"]},
        ),
        f"{root}/VTIPerBeat/value": with_attrs(
            vti_per_beat,
            {"unit": "mm", "dimDesc": ["beat"]},
        ),
        f"{root}/VTIMean/value": with_attrs(_nanmean_or_nan(vti_per_beat), {"unit": "mm"}),
        f"{root}/VTIStd/value": with_attrs(_nanstd_or_nan(vti_per_beat), {"unit": "mm"}),
        f"{root}/beatPeriodIdx/value": with_attrs(
            beat_period_idx.astype(np.float32, copy=False).reshape(1, -1),
            {"unit": "frame", "dimDesc": ["beat"]},
        ),
        f"{root}/beatPeriodSeconds/value": with_attrs(
            beat_period_seconds.reshape(1, -1),
            {"unit": "s", "dimDesc": ["beat"]},
        ),
    }


def _vti_per_beat(
    signal: np.ndarray,
    dt_seconds: np.float32,
) -> np.ndarray:
    if signal.size == 0:
        return np.empty(0, dtype=np.float32)
    return (np.nansum(signal, axis=1) * dt_seconds).astype(np.float32, copy=False)


def _waveform_stats(
    signal: np.ndarray,
    systole_indexes: np.ndarray,
) -> dict[str, np.float32]:
    cycle = _mean_interpolated_cycle(signal, systole_indexes, n_samples=60)
    vmax = np.float32(np.nanmax(cycle))
    vmin = np.float32(np.nanmin(cycle))
    vmean = np.float32(np.nanmean(cycle))
    return {
        "mean": vmean,
        "max": vmax,
        "min": vmin,
        "resistivity_index": _clamped_ratio(vmax - vmin, vmax, 0.0, 1.0),
        "pulsatility_index": _clamped_ratio(vmax - vmin, vmean, 0.0, np.inf),
        "max_min_ratio": _safe_ratio(vmax, vmin),
    }


def _mean_interpolated_cycle(
    signal: np.ndarray,
    systole_indexes: np.ndarray,
    *,
    n_samples: int,
) -> np.ndarray:
    values = np.asarray(signal, dtype=np.float32).reshape(-1)
    indexes = np.asarray(systole_indexes, dtype=np.int32).reshape(-1)
    if values.size == 0:
        return np.full(n_samples, np.nan, dtype=np.float32)
    if indexes.size < 2:
        return np.interp(
            np.linspace(0, values.size - 1, n_samples, dtype=np.float32),
            np.arange(values.size, dtype=np.float32),
            values,
        ).astype(np.float32)

    cycles = []
    for start, stop in zip(indexes[:-1], indexes[1:], strict=False):
        start = int(start)
        stop = int(stop)
        if stop <= start + 1:
            continue
        x = np.arange(start, stop, dtype=np.float32)
        cycles.append(
            np.interp(
                np.linspace(start, stop - 1, n_samples, dtype=np.float32),
                x,
                values[start:stop],
            ).astype(np.float32)
        )
    if not cycles:
        return np.full(n_samples, np.nan, dtype=np.float32)
    return np.nanmean(np.stack(cycles, axis=0), axis=0).astype(np.float32)


def _clamped_ratio(
    numerator: np.float32,
    denominator: np.float32,
    lower: float,
    upper: float,
) -> np.float32:
    value = _safe_ratio(numerator, denominator)
    if not np.isfinite(value):
        return value
    return np.float32(np.clip(value, lower, upper))


def _safe_ratio(numerator: np.float32, denominator: np.float32) -> np.float32:
    if not np.isfinite(denominator) or denominator == 0:
        return np.float32(np.nan)
    return np.float32(numerator / denominator)


def _nanmean_or_nan(values: np.ndarray) -> np.float32:
    if values.size == 0:
        return np.float32(np.nan)
    return np.float32(np.nanmean(values))


def _nanstd_or_nan(values: np.ndarray) -> np.float32:
    if values.size == 0:
        return np.float32(np.nan)
    return np.float32(np.nanstd(values))
