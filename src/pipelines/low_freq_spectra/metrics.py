"""Peak, valley, and peak-to-peak area metrics for gap spectra."""

from __future__ import annotations

from numbers import Integral

import numpy as np
from scipy.integrate import trapezoid

from calculations.math import nanmedian
from calculations.math.frequency import frequency_bin_spacing_khz
from calculations.math.phase_windows import (
    DEFAULT_EARLY_PHASE_INTERVAL, DEFAULT_LATE_PHASE_INTERVAL,
    early_late_phase_windows, phase_window_medians, late_minus_early,
)
from input_output.schema import EyeFlowOutputPaths
from input_output.schema.holodoppler import HD_SPECTROGRAM_FREQUENCY_START_KEY, HD_SPECTROGRAM_FREQUENCY_END_KEY
from utils.logger import Logger

DEFAULT_VALLEY_NEIGHBORHOOD_FRACTION = 1 / 4
PHASE_METRIC_NAMES = ("A0", "A1", "A2", "R1", "R2", "Z1", "Z2")

SPECTRUM_METRIC_DEFINITIONS = {
    "index_max1": "zero-based index of the maximum in the first half of the bin domain",
    "max1": "maximum spectrum value in the first half of the bin domain, ignoring NaNs",
    "index_max2": "zero-based index of the maximum in the second half of the bin domain",
    "max2": "maximum spectrum value in the second half of the bin domain, ignoring NaNs",
    "index_min1": "zero-based central valley index floor(bin_count/2) - 1",
    "min1": "spectrum value at zero-based index floor(bin_count/2) - 1",
    "index_min2": "zero-based central valley index floor(bin_count/2)",
    "min2": "spectrum value at zero-based index floor(bin_count/2)",
    "A0": "NaN-aware central-neighborhood median height times peak-to-peak width in kHz",
    "A1": "trapezoidal spectrum integral from index_max1 to index_max2 with calibrated kHz center spacing; NaNs propagate",
    "A2": "arithmetic mean of max1 and max2 times peak-to-peak width in kHz",
    "R1": "A1 divided by A0; zero or missing A0 produces NaN",
    "R2": "A2 divided by A0; zero or missing A0 produces NaN",
    "Z1": "(A1 - A0) divided by A0; zero or missing A0 produces NaN",
    "Z2": "(A2 - A0) divided by A0; zero or missing A0 produces NaN",
}

PHASE_METRIC_DEFINITIONS = {
    "A0": "fixed peak-to-peak width in kHz times NaN-aware central-neighborhood median per phase",
    "A1": "trapezoidal spectrum integral between fixed peak bins per phase with calibrated kHz center spacing; interval NaNs propagate",
    "A2": "fixed peak-to-peak width in kHz times mean spectrum height at those two bins per phase",
    "R1": "A1_t divided by A0_t per phase; zero or missing A0_t produces NaN",
    "R2": "A2_t divided by A0_t per phase; zero or missing A0_t produces NaN",
    "Z1": "(A1_t - A0_t) divided by A0_t per phase; zero or missing A0_t produces NaN",
    "Z2": "(A2_t - A0_t) divided by A0_t per phase; zero or missing A0_t produces NaN",
}


def analyze_spectrum(
    spectrum, *, sampling_frequency_hz: float, valley_neighborhood_length: int | None = None,
    fft_bin_count: int | None = None, binning_factor: int = 1,
    frequency_start: float = 0.0, frequency_end: float = 0.5,
) -> dict[str, float]:
    """Analyze a spectrum in its stored frequency order.

    Select one maximum in each half, including endpoints and ignoring NaNs.
    Ties select the lowest index. An all-NaN half has a NaN index and value.
    Valley indexes are the two central bins, even if their values are NaN.
    The peak-height mean requires both inputs. Areas span the peak-index distance
    times the calibrated FFT-average spacing in kHz without wrapping or reordering. A1 integrates samples
    including both peak bins using that positive spacing and the trapezoidal
    rule; interval NaNs propagate.
    The even-length valley neighborhood is centered between the central bins;
    its default length is one quarter of the bin count, rounded to the nearest
    even number (ties upward), with at least two bins. Its median ignores NaNs, returning
    NaN if all samples are missing. A0 uses the valley median and A2 uses the
    peak-height mean, both with the same peak-to-peak width in kHz.
    Any one-dimensional domain with at least two bins is supported. Odd domains
    split at floor(bin_count/2). Unsupported shapes return NaN metrics with a warning.
    """
    values = np.asarray(spectrum)
    result = {name: float("nan") for name in SPECTRUM_METRIC_DEFINITIONS}
    if not _supported_bin_domain(values):
        Logger.log_warning(
            "Gap spectrum metrics require a one-dimensional bin domain with at least "
            f"2 bins; received shape {values.shape}. "
            "Returning NaN for all spectrum metrics."
        )
        return result
    values = np.asarray(values, dtype=np.float64)
    spacing = frequency_bin_spacing_khz(
        values.size, sampling_frequency_hz, fft_bin_count=fft_bin_count, binning_factor=binning_factor,
        frequency_start=frequency_start, frequency_end=frequency_end,
    )
    middle = values.size // 2
    result["index_max1"], result["max1"] = _peak(values, 0, middle)
    result["index_max2"], result["max2"] = _peak(values, middle, values.size)
    for name, start, stop in (("1", 0, middle), ("2", middle, values.size)):
        if np.isnan(result[f"index_max{name}"]):
            Logger.log_warning(
                f"Gap spectrum bins [{start}, {stop}) contain only NaN samples; "
                f"index_max{name} and max{name} remain NaN."
            )
    result.update(
        index_min1=float(middle - 1), min1=float(values[middle - 1]),
        index_min2=float(middle), min2=float(values[middle]),
    )
    h0 = (result["max1"] + result["max2"]) / 2
    hv = float(spectrum_valley_height(values, valley_neighborhood_length=valley_neighborhood_length))
    if np.isfinite(result["index_max1"]) and np.isfinite(result["index_max2"]):
        start = int(result["index_max1"])
        stop = int(result["index_max2"]) + 1
        width = (stop - start - 1) * spacing
        result["A0"] = hv * width
        result["A1"] = float(trapezoid(values[start:stop], dx=spacing))
        result["A2"] = h0 * width
        if result["A0"] != 0:
            result["R1"] = result["A1"] / result["A0"]
            result["R2"] = result["A2"] / result["A0"]
            result["Z1"] = (result["A1"] - result["A0"]) / result["A0"]
            result["Z2"] = (result["A2"] - result["A0"]) / result["A0"]
    return result


def _supported_bin_domain(values: np.ndarray) -> bool:
    return values.ndim == 1 and values.size >= 2


def _valley_neighborhood_indices(length: int | None, bin_count: int) -> np.ndarray:
    if length is None:
        length = max(2, 2 * int(bin_count * DEFAULT_VALLEY_NEIGHBORHOOD_FRACTION / 2 + 0.5))
    if (
        isinstance(length, bool) or not isinstance(length, Integral)
        or length < 2 or length > bin_count or length % 2
    ):
        raise ValueError(
            f"valley_neighborhood_length must be an even integer from 2 through {bin_count}."
        )
    middle = bin_count // 2
    return np.arange(middle - length // 2, middle + length // 2)


def spectrum_valley_height(spectrum, *, valley_neighborhood_length: int | None = None):
    """Return Hv over the central quarter of the last (frequency) axis, ignoring NaNs."""
    values = np.asarray(spectrum, dtype=np.float64)
    if values.ndim < 1 or values.shape[-1] < 2:
        raise ValueError("Valley height requires a frequency axis with at least two bins.")
    bins = _valley_neighborhood_indices(valley_neighborhood_length, values.shape[-1])
    return nanmedian(values[..., bins], axis=-1, dtype=np.float64)


def _peak(values: np.ndarray, start: int, stop: int) -> tuple[float, float]:
    half = values[start:stop]
    valid = np.flatnonzero(~np.isnan(half))
    if valid.size == 0:
        return float("nan"), float("nan")
    index = start + int(valid[np.argmax(half[valid])])
    return float(index), float(values[index])


def pack_spectrum_metrics(
    spectrum, *, source_attrs: dict[str, object], sampling_frequency_hz: float,
    valley_neighborhood_length: int | None = None,
    phase_spectrum=None, phase_attrs: dict[str, object] | None = None,
    early_interval=DEFAULT_EARLY_PHASE_INTERVAL, late_interval=DEFAULT_LATE_PHASE_INTERVAL,
    frequency_calibration_attrs: dict[str, object] | None = None,
) -> dict[str, object]:
    """Pack scalar metrics with units and references to the aggregate spectrum."""
    schema = EyeFlowOutputPaths.active()
    calibration_attrs = dict(frequency_calibration_attrs or {})
    calibration = dict(
        fft_bin_count=calibration_attrs.get("fft_bin_count"),
        binning_factor=calibration_attrs.get("frequency_binning_factor", 1),
        frequency_start=calibration_attrs.get(HD_SPECTROGRAM_FREQUENCY_START_KEY, 0.0),
        frequency_end=calibration_attrs.get(HD_SPECTROGRAM_FREQUENCY_END_KEY, 0.5),
    )
    result = analyze_spectrum(
        spectrum, sampling_frequency_hz=sampling_frequency_hz,
        valley_neighborhood_length=valley_neighborhood_length,
        **calibration,
    )
    values = np.asarray(spectrum)
    valley_bins = (
        _valley_neighborhood_indices(valley_neighborhood_length, values.size)
        if _supported_bin_domain(values) else np.array([], dtype=int)
    )
    attrs = {
        **calibration_attrs,
        "source_path": "/" + schema.spectrograms_masked_time_aggregated_median,
        "source_bin_count": int(np.size(spectrum)),
        "minimum_bin_count": 2,
        "bin_domain_split_index": values.size // 2 if _supported_bin_domain(values) else 0,
        "peak_half_valid_bin_counts": (
            [
                int(np.count_nonzero(~np.isnan(half)))
                for half in np.array_split(values, [values.size // 2])
            ]
            if _supported_bin_domain(values) else []
        ),
        "dimDesc": [],
        "peak_selection": "maximum_per_bin_domain_half_ignore_NaN_lowest_index_on_ties",
        "interval_endpoints": "both_included",
        "area_interval_width": "(index_max2 - index_max1) * frequency_bin_spacing_khz",
        "area_coordinate": "unwrapped_frequency_khz_in_HD_bin_order",
        "A1_integration_method": "trapezoidal_khz_spacing",
        "sampling_frequency_hz": sampling_frequency_hz,
        "sampling_frequency_source_path": "/" + schema.sampling_frequency_hz,
        "frequency_bin_spacing_khz": (
            frequency_bin_spacing_khz(values.size, sampling_frequency_hz, **calibration)
            if _supported_bin_domain(values) else np.nan
        ),
        "frequency_bin_order": calibration_attrs.get(
            "frequency_bin_order", source_attrs.get("frequency_bin_order", "HD_export_order"),
        ),
        "nonfinite_value_policy": "missing_required_values_and_interval_NaNs_propagate",
        "valley_neighborhood_length": len(valley_bins),
        "valley_neighborhood_indices": valley_bins,
        "valley_neighborhood_selection": (
            "quarter_of_bin_domain" if valley_neighborhood_length is None else "explicit_length"
        ),
        "valley_neighborhood_default_rounding": "nearest_even_ties_up_minimum_2",
        "valley_aggregation_method": "nanmedian",
        "valley_nonfinite_value_policy": "ignore_NaN_all_NaN_neighborhood_returns_NaN",
    }
    output = {}
    for name, value in result.items():
        unit = (
            "1"
            if name.startswith("index_") or name in {"R1", "R2", "Z1", "Z2"}
            else source_attrs.get("unit", "a.u.")
        )
        metric_attrs = {**attrs, "unit": unit, "definition": SPECTRUM_METRIC_DEFINITIONS[name]}
        if name in {"A0", "A1", "A2"}:
            metric_attrs["unit"] = f"{unit}*kHz"
        if name.startswith("index_"):
            metric_attrs["index_base"] = 0
            metric_attrs["index_storage"] = "integer_valued_float32_with_NaN_for_missing_index"
        output[f"{schema.spectrograms_masked_metrics_root}/{name}"] = (
            np.float32(value),
            metric_attrs,
        )
    if phase_spectrum is not None:
        if phase_attrs is None:
            raise ValueError("Phase-spectrum attributes must provide phase_samples.")
        if np.ndim(phase_spectrum) != 2 or np.shape(phase_spectrum)[1] != values.size:
            raise ValueError("Phase and pooled spectra must share the same frequency-bin domain.")
        output.update(_pack_phase_metrics(
            phase_spectrum, fixed_metrics=result, source_attrs=phase_attrs,
            metric_attrs=attrs, valley_neighborhood_length=valley_neighborhood_length,
            sampling_frequency_hz=sampling_frequency_hz,
            calibration=calibration,
            early_interval=early_interval, late_interval=late_interval,
        ))
    return output


def analyze_phase_spectrum(
    spectrum, *, index_max1: float, index_max2: float, sampling_frequency_hz: float,
    valley_neighborhood_length: int | None = None,
    fft_bin_count: int | None = None, binning_factor: int = 1,
    frequency_start: float = 0.0, frequency_end: float = 0.5,
) -> dict[str, np.ndarray]:
    """Calculate phase curves with peak bins fixed from the pooled spectrum.

    Input axes are (phase, freq). A0 uses the central neighborhood median;
    A1 is the trapezoidal integral between fixed peaks with calibrated kHz
    spacing; A2 uses the mean height at those peaks. All areas span their index
    distance times FFT-average center spacing in kHz, without sorting the signed display axis.
    Missing peak indexes leave all curves NaN. At each phase, missing peak
    heights invalidate A2, interval NaNs invalidate A1, and A0 ignores valley
    NaNs unless the whole neighborhood is missing. Zero A0 yields NaN ratios
    and offsets.
    """
    values = np.asarray(spectrum, dtype=np.float64)
    if values.ndim != 2 or any(size == 0 for size in values.shape):
        raise ValueError("Phase spectra must have nonempty (phase, freq) axes.")
    result = {name: np.full(values.shape[0], np.nan) for name in PHASE_METRIC_NAMES}
    if values.shape[1] < 2 or np.isnan(index_max1) or np.isnan(index_max2):
        return result
    if not (
        np.isfinite(index_max1) and np.isfinite(index_max2)
        and index_max1 == int(index_max1) and index_max2 == int(index_max2)
        and 0 <= index_max1 < index_max2 < values.shape[1]
    ):
        raise ValueError("Fixed peak indices must be ordered integer bin indices within the spectrum.")
    start, stop = int(index_max1), int(index_max2) + 1
    spacing = frequency_bin_spacing_khz(
        values.shape[1], sampling_frequency_hz, fft_bin_count=fft_bin_count, binning_factor=binning_factor,
        frequency_start=frequency_start, frequency_end=frequency_end,
    )
    width = (stop - start - 1) * spacing
    result["A0"] = width * spectrum_valley_height(values, valley_neighborhood_length=valley_neighborhood_length)
    result["A1"] = trapezoid(values[:, start:stop], dx=spacing, axis=1)
    result["A2"] = width * (values[:, start] + values[:, stop - 1]) / 2
    for name, numerator in (
        ("R1", result["A1"]), ("R2", result["A2"]),
        ("Z1", result["A1"] - result["A0"]), ("Z2", result["A2"] - result["A0"]),
    ):
        np.divide(
            numerator, result["A0"], out=result[name],
            where=(result["A0"] != 0) & ~np.isnan(result["A0"]),
        )
    return result


def _pack_phase_metrics(
    spectrum, *, fixed_metrics: dict[str, float], source_attrs: dict[str, object],
    metric_attrs: dict[str, object], valley_neighborhood_length: int | None,
    sampling_frequency_hz: float,
    calibration: dict[str, object],
    early_interval, late_interval,
) -> dict[str, object]:
    schema = EyeFlowOutputPaths.active()
    root = schema.spectrograms_masked_metrics_root
    phase = np.asarray(source_attrs["phase_samples"])
    if phase.dtype.kind != "f":
        phase = phase.astype(np.float64)
    if phase.shape != (np.shape(spectrum)[0],):
        raise ValueError("Phase samples must be strictly increasing coordinates in [0, 1).")
    windows = early_late_phase_windows(phase, early_interval=early_interval, late_interval=late_interval)
    curves = analyze_phase_spectrum(
        spectrum, index_max1=fixed_metrics["index_max1"], index_max2=fixed_metrics["index_max2"],
        valley_neighborhood_length=valley_neighborhood_length,
        sampling_frequency_hz=sampling_frequency_hz,
        **calibration,
    )
    attrs = {
        **{key: value for key, value in source_attrs.items() if not key.startswith("aggregation_")},
        **metric_attrs,
        "source_path": "/" + schema.spectrograms_masked_time_aggregated_median_beat,
        "peak_source_path": "/" + schema.spectrograms_masked_time_aggregated_median,
        "fixed_peak_index_source_paths": [f"/{root}/index_max1", f"/{root}/index_max2"],
        "fixed_peak_indices": [fixed_metrics["index_max1"], fixed_metrics["index_max2"]],
        "peak_selection": "fixed_pooled_spectrum_indices_at_all_phases",
        "calculation_method": "fixed_peak_bin_geometry_on_beat_median_spectrum",
        "nonfinite_value_policy": "missing_peak_heights_and_interval_NaNs_propagate_valley_ignores_NaN",
        "phase_samples": phase,
        "phase_endpoint": "excluded",
        "dimDesc": ["phase"],
    }
    output = {}
    for name, curve in curves.items():
        unit = "1" if name in {"R1", "R2", "Z1", "Z2"} else source_attrs.get("unit", "a.u.")
        if name in {"A0", "A1", "A2"}:
            unit = f"{unit}*kHz"
        curve_path = f"{root}/{name}_t"
        curve_attrs = {
            **attrs, "unit": unit,
            "definition": PHASE_METRIC_DEFINITIONS[name],
        }
        output[curve_path] = (curve.astype(np.float32), curve_attrs)
        scalar_attrs = {
            key: value for key, value in curve_attrs.items()
            if key not in {"phase_samples", "phase_endpoint"}
        }
        scalar_attrs.update(source_path="/" + curve_path, dimDesc=[])
        aggregates = phase_window_medians(curve, windows)
        for window, (mask, interval) in windows.items():
            aggregate = float(aggregates[window])
            output[f"{root}/{name}_{window}"] = (np.float32(aggregate), {
                **scalar_attrs,
                "definition": f"NaN-aware median of {name}_t over the {window} phase interval",
                "aggregation_method": "nanmedian",
                "aggregation_axes": ["phase"],
                "aggregation_weighting": "equal_weight_per_valid_phase_sample",
                "phase_interval": interval,
                "phase_interval_endpoints": "both_included",
                "source_phase_sample_count": int(np.count_nonzero(mask)),
                "valid_sample_count": int(np.count_nonzero(~np.isnan(curve[mask]))),
                "nonfinite_value_policy": "ignore_NaN_all_NaN_window_returns_NaN",
            })
        output[f"{root}/delta_{name}"] = (np.float32(late_minus_early(aggregates["early"], aggregates["late"])), {
            **scalar_attrs,
            "definition": f"{name}_late minus {name}_early",
            "source_paths": [f"/{root}/{name}_late", f"/{root}/{name}_early"],
            "calculation_method": "late_median_minus_early_median",
            "early_phase_interval": windows["early"][1],
            "late_phase_interval": windows["late"][1],
            "nonfinite_value_policy": "missing_early_or_late_median_propagates_NaN",
        })
    return output
