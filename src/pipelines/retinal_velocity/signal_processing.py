"""Signal preparation for the retinal-velocity core pipeline."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from calculations.blood_flow_velocity.signal_analysis.cardiac_cycle import (
    CardiacCycleAnalysis,
)
from calculations.math import butter_lowpass_filtfilt

from .models import RetinalVelocity

DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ = 15.0


def build_retinal_velocity(
    estimation: Mapping[str, object],
    cardiac_cycle: CardiacCycleAnalysis,
    cardiac_cycle_source: str,
    *,
    dt_seconds: float,
    index_base: int = 0,
    lowpass_freq_hz: float = DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ,
) -> RetinalVelocity:
    """Build the canonical typed result from an estimator output."""

    artery_raw = _float_array(estimation, "retinal_artery_velocity_signal")
    vein_raw = _float_array(estimation, "retinal_vein_velocity_signal")
    artery_filtered = _filtered_vessel_signal(
        artery_raw,
        "artery",
        cardiac_cycle,
        cardiac_cycle_source,
        dt_seconds,
        lowpass_freq_hz,
    )
    vein_filtered = _filtered_vessel_signal(
        vein_raw,
        "vein",
        cardiac_cycle,
        cardiac_cycle_source,
        dt_seconds,
        lowpass_freq_hz,
    )
    return RetinalVelocity(
        artery_velocity_raw=artery_raw,
        vein_velocity_raw=vein_raw,
        artery_velocity_filtered=artery_filtered,
        vein_velocity_filtered=vein_filtered,
        velocity_map=estimation.get("velocity_map"),
        moment0_average=_float_array(estimation, "moment0_avg"),
        velocity_average=_float_array(estimation, "velocity_map_avg"),
        frms_average=_float_array(estimation, "fRMS_avg"),
        frms_background_average=_float_array(estimation, "fRMS_bkg_avg"),
        delta_frms_average=_float_array(estimation, "deltafRMS_avg"),
        velocity_section_mask=np.asarray(
            estimation["velocity_section_mask"],
            dtype=bool,
        ),
        artery_frms=_float_array(estimation, "retinal_artery_fRMS_signal"),
        vein_frms=_float_array(estimation, "retinal_vein_fRMS_signal"),
        artery_frms_background=_float_array(
            estimation,
            "retinal_artery_fRMS_bkg_signal",
        ),
        vein_frms_background=_float_array(
            estimation,
            "retinal_vein_fRMS_bkg_signal",
        ),
        vessel_frms_background=_float_array(
            estimation,
            "retinal_vessel_fRMS_bkg_signal",
        ),
        artery_delta_frms=_float_array(
            estimation,
            "retinal_artery_deltafRMS_signal",
        ),
        vein_delta_frms=_float_array(
            estimation,
            "retinal_vein_deltafRMS_signal",
        ),
        cardiac_cycle=cardiac_cycle,
        cardiac_cycle_source=cardiac_cycle_source,
        dt_seconds=float(dt_seconds),
        provenance=_estimation_provenance(estimation),
        index_base=int(index_base),
    )


def _filter(signal, dt_seconds: float, lowpass_freq_hz: float) -> np.ndarray:
    return butter_lowpass_filtfilt(
        signal,
        dt_seconds=np.float32(dt_seconds),
        lowpass_freq_hz=np.float32(lowpass_freq_hz),
        order=4,
    )


def _filtered_vessel_signal(
    signal: np.ndarray,
    vessel: str,
    cardiac_cycle: CardiacCycleAnalysis,
    cardiac_cycle_source: str,
    dt_seconds: float,
    lowpass_freq_hz: float,
) -> np.ndarray:
    if cardiac_cycle_source == vessel:
        return np.asarray(cardiac_cycle.systole.signal_filtered, dtype=np.float32)
    return _filter(signal, dt_seconds, lowpass_freq_hz)


def _float_array(values: Mapping[str, object], key: str) -> np.ndarray:
    return np.asarray(values[key], dtype=np.float32)


def _estimation_provenance(
    estimation: Mapping[str, object],
) -> dict[str, object]:
    keys = {
        "velocity_estimation_method",
        "velocity_quantity",
        "velocity_unit",
        "laser_wavelength_m",
        "numerical_aperture",
        "band_ratio_frequency_scale_hz",
        "band_ratio_calibration_model",
        "band_ratio_calibration_source",
        "band_ratio_calibration_version",
        "band_lf_source_path",
        "band_hf_source_path",
    }
    return {
        key: value
        for key, value in estimation.items()
        if value is not None and (key in keys or key.startswith("band_lf_"))
    }


__all__ = [
    "DEFAULT_VELOCITY_SIGNAL_LOWPASS_HZ",
    "build_retinal_velocity",
]
