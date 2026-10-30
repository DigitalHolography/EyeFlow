"""Pack selected gap-spectrum shape descriptors for regular EyeFlow runs."""

from __future__ import annotations

import numpy as np
from scipy.stats import wasserstein_distance

from calculations.math.phase_windows import (
    DEFAULT_EARLY_PHASE_INTERVAL, DEFAULT_LATE_PHASE_INTERVAL,
    early_late_phase_windows, phase_window_medians, late_minus_early,
)
from calculations.math.statistics import nanmedian
from calculations.math.spectra import PERCENTILES, describe_spectra, frequency_distribution, interval_area
from input_output.schema import EyeFlowOutputPaths

VERTICAL_PERCENTAGES = (5, 25, 50)
SHAPE_PHASE_METRIC_NAMES = tuple(
    name for p in VERTICAL_PERCENTAGES for name in (f"Av_{p}_t", f"SEv_{p}")
) + ("SpectralEntropy_t", "F50_t", "F90_t")


def pack_distribution_metrics(
    phase_spectrum, *, pooled_spectrum, phase_attrs, frequencies_khz, peak_indices,
    early_interval=DEFAULT_EARLY_PHASE_INTERVAL, late_interval=DEFAULT_LATE_PHASE_INTERVAL,
):
    """Pack probabilities, entropy and magnitude-frequency transport descriptors.

    Phase spectra have (phase, freq) axes; coordinates are calibrated signed
    kHz centers. Full discrete bin masses include both fixed pooled peaks.
    Summary curves take medians after per-phase evaluation; Wasserstein uses
    separately normalized early/late median spectra instead.
    """
    schema = EyeFlowOutputPaths.active()
    ROOT = schema.spectrograms_masked_metrics_root
    PHASE_PATH = schema.spectrograms_masked_time_aggregated_median_beat
    POOLED_PATH = schema.spectrograms_masked_time_aggregated_median
    values, pooled = np.asarray(phase_spectrum, dtype=np.float64), np.asarray(pooled_spectrum)
    if values.ndim != 2 or pooled.ndim != 1 or values.shape[1] != pooled.size:
        raise ValueError("MedianBeat must have (phase, freq) axes matching MedianPhaseBeat.")
    attrs = dict(phase_attrs)
    phase = np.asarray(attrs["phase_samples"])
    if phase.shape != values.shape[:1]:
        raise ValueError("phase_samples must match MedianBeat.")
    windows = early_late_phase_windows(phase, early_interval=early_interval, late_interval=late_interval)
    frequencies = np.asarray(frequencies_khz, dtype=np.float64)
    if frequencies.shape != pooled.shape or not np.isfinite(frequencies).all():
        raise ValueError("Distribution metrics require matching finite calibrated kHz coordinates.")
    peaks = np.asarray(peak_indices, dtype=float)
    if peaks.shape != (2,):
        raise ValueError("Exactly two pooled peak indices are required.")
    valid_peaks = (
        np.isfinite(peaks).all() and np.equal(peaks, np.floor(peaks)).all()
        and 0 <= peaks[0] < pooled.size // 2 <= peaks[1] < pooled.size
    )
    # Keep a fixed schema even when a pooled peak is missing.
    indices = np.arange(int(peaks[0]), int(peaks[1]) + 1) if valid_peaks else np.arange(pooled.size)
    selected = values[:, indices] if valid_peaks else np.full((values.shape[0], indices.size), np.nan)
    descriptions = describe_spectra(selected, frequencies[indices])
    common = {
        **{key: value for key, value in attrs.items()
           if key.startswith(("cycle_boundary_", "beat_duration_", "frequency_", "fft_", "sampling_", "spectrogram_frequency_"))},
        "metric_generator": "pipelines.low_freq_spectra.shape_metrics",
        "distribution_metric_version": 1,
        "source_path": "/" + PHASE_PATH,
        "peak_source_path": "/" + POOLED_PATH,
        "fixed_peak_indices": peaks,
        "included_bin_indices": indices,
        "included_bin_count": len(indices),
        "frequency_source_path": "/" + schema.spectrograms_frequency_bin_table,
        "frequency_coordinate": "absolute_calibrated_frequency_khz",
        "support_definition": "inclusive_fixed_pooled_peak_interval",
        "probability_measure": "full_HD_band_power_at_bin_centers_including_both_endpoint_bins",
        "valley_normalization": "divide_by_true_minimum_inside_fixed_peak_interval_without_subtraction",
        "probability_normalization": "R/sum(R)=power/sum(power); valley_scale_cancels",
        "zero_valley_policy": "R_is_NaN_but_probability_descriptors_remain_defined_if_total_power_is_positive",
        "nonfinite_value_policy": "any_nonfinite_or_negative_included_power_invalidates_that_phase",
    }
    output = {}

    def put(name, data, *, unit="1", dimensions=(), definition, **metadata):
        data = np.asarray(data, dtype=np.float64)
        if np.any(np.isfinite(data) & (np.abs(data) > np.finfo(np.float32).max)):
            raise ValueError(f"{name} would overflow float32.")
        path = f"{ROOT}/{name}"
        output[path] = (data.astype(np.float32), {
            **common, "unit": unit, "dimDesc": list(dimensions),
            "definition": definition, **metadata,
        })

    put("PowerQuantilePercent", PERCENTILES, unit="%", dimensions=("power_percentile",),
        definition="Cumulative power percentage at which absolute frequency is reported")
    put("DistributionFrequencyMagnitude", np.abs(frequencies), unit="kHz", dimensions=("freq",),
        definition="Absolute calibrated bin centers in original stored frequency order")
    put("ValleyMinimum_t", descriptions["ValleyMinimum"], unit=attrs.get("unit", "a.u."),
        dimensions=("phase",), phase_samples=phase,
        definition="True minimum included spectral power per phase; distinct from the central-quarter Hv median")
    for name in ("ValleyNormalizedSpectrum", "SpectralProbability"):
        full = np.full(values.shape, np.nan)
        full[:, indices] = descriptions[name]
        definition = (
            "Power divided by its phasewise true spectral minimum; finite values are >=1"
            if name == "ValleyNormalizedSpectrum" else
            "Discrete bin probability power/sum(power); finite included probabilities sum to 1 per phase"
        )
        put(name + "_t", full, dimensions=("phase", "freq"), phase_samples=phase,
            excluded_bin_policy="NaN_outside_fixed_peak_support", definition=definition)

    def pack_curve(name, curve, *, unit, dimensions=(), **details):
        put(name + "_t", curve, unit=unit, dimensions=("phase", *dimensions),
            phase_samples=phase, definition=f"Phasewise {name}", **details)
        aggregate_details = {
            **details, "source_path": f"/{ROOT}/{name}_t",
            "aggregation_method": "nanmedian", "aggregation_axes": ["phase"],
            "aggregation_weighting": "equal_weight_per_valid_phase_sample",
        }
        put(name, nanmedian(curve, axis=0, dtype=np.float64), unit=unit, dimensions=dimensions,
            definition=f"Median of phasewise {name} over the full cardiac cycle",
            phase_interval=[0., 1.], valid_phase_sample_count=int(np.isfinite(curve.reshape(curve.shape[0], -1)).all(axis=1).sum()),
            **aggregate_details)
        medians = phase_window_medians(curve, windows)
        for window, (mask, interval) in windows.items():
            put(name + "_" + window, medians[window], unit=unit, dimensions=dimensions,
                definition=f"Median of phasewise {name} in the closed {window} phase window",
                phase_interval=interval, phase_interval_endpoints="both_included",
                source_phase_sample_count=int(mask.sum()),
                valid_phase_sample_count=int(np.isfinite(curve[mask].reshape(int(mask.sum()), -1)).all(axis=1).sum()) if mask.any() else 0,
                **aggregate_details)
        put("delta_" + name, medians["late"] - medians["early"], unit=unit, dimensions=dimensions,
            definition=f"{name}_late minus {name}_early", **details,
            source_paths=[f"/{ROOT}/{name}_late", f"/{ROOT}/{name}_early"],
            early_phase_interval=windows["early"][1], late_phase_interval=windows["late"][1])

    pack_curve("SpectralEntropy", descriptions["SpectralEntropy"], unit="bits",
               entropy_definition="-sum(p_i*log2(p_i)); zero_probability_contributes_zero",
               entropy_bin_domain="original_signed_frequency_bins_inside_fixed_peaks",
               entropy_maximum_bits=float(np.log2(indices.size)))
    quantile_details = {
        "quantile_method": "discrete_inverse_CDF_no_sub_bin_interpolation",
        "power_percentile_source_path": f"/{ROOT}/PowerQuantilePercent",
        "quantile_endpoint_definition": "0_and_100_are_smallest_and_largest_positive_mass_frequency_magnitudes",
    }
    pack_curve("PowerQuantileFrequency", descriptions["PowerQuantileFrequency"], unit="kHz",
               dimensions=("power_percentile",), **quantile_details)
    for percentile in (50, 90):
        pack_curve(f"F{percentile}", descriptions["PowerQuantileFrequency"][:, percentile], unit="kHz",
                   power_percentile=percentile, **quantile_details)

    # Use the same phase-median spectra as the existing early/late spectrum
    # outputs, then normalize each distribution separately. Do not take a
    # difference of scalar entropy or F50 values as a substitute for W1.
    window_spectra = phase_window_medians(selected, windows)
    window_distributions = {}
    x = np.unique(np.abs(frequencies[indices]))
    for window, signal in window_spectra.items():
        description = describe_spectra(signal, frequencies[indices])
        probability = description["SpectralProbability"]
        window_distributions[window] = probability
        full = np.full(pooled.shape, np.nan)
        full[indices] = probability
        put("Wasserstein" + window.title() + "Probability", full, dimensions=("freq",),
            definition=f"Normalized distribution of the {window} phase-median spectrum",
            phase_interval=windows[window][1], aggregation_order="median_spectrum_then_probability_normalization")
        _, masses = frequency_distribution(np.abs(frequencies[indices]), probability)
        cumulative = np.cumsum(masses)
        if np.isfinite(cumulative).all():
            cumulative[-1] = 1.
        put("Wasserstein" + window.title() + "CDF", cumulative, dimensions=("frequency_magnitude",),
            definition=f"CDF of the {window} normalized spectrum on sorted absolute frequencies",
            coordinate_source_path=f"/{ROOT}/WassersteinFrequencyMagnitude")
        put("Wasserstein" + window.title() + "QuantileFrequency", description["PowerQuantileFrequency"],
            unit="kHz", dimensions=("power_percentile",),
            definition=f"Quantile curve of the {window} phase-median spectrum", **quantile_details)
    put("WassersteinFrequencyMagnitude", x, unit="kHz", dimensions=("frequency_magnitude",),
        definition="Sorted unique calibrated frequency magnitudes for early/late CDF comparison")
    early, late = window_distributions["early"], window_distributions["late"]
    distance = np.nan
    if np.isfinite(early).all() and np.isfinite(late).all():
        distance = wasserstein_distance(np.abs(frequencies[indices]), np.abs(frequencies[indices]), early, late)
    put("WassersteinEarlyLate", distance, unit="kHz",
        definition="Wasserstein-1 distance between normalized early/late phase-median spectra on absolute calibrated frequency",
        comparison="early_vs_late", direction="unsigned",
        distance_interpretation="mean_frequency_displacement_required_to_transform_the_early_power_distribution_into_the_late_one",
        early_phase_interval=windows["early"][1], late_phase_interval=windows["late"][1],
        early_probability_source_path=f"/{ROOT}/WassersteinEarlyProbability",
        late_probability_source_path=f"/{ROOT}/WassersteinLateProbability")
    return output


def pack_vertical_metrics(
    pooled_spectrum, phase_spectrum, *, phase_attrs, base_metrics, percentages=VERTICAL_PERCENTAGES,
    early_interval=DEFAULT_EARLY_PHASE_INTERVAL, late_interval=DEFAULT_LATE_PHASE_INTERVAL,
):
    """Pack inward peak-interval integrals and temporal-median area fractions.

    Each width is a percentage of its own peak-to-half-domain boundary span.
    The denominator is the existing calibrated phasewise A_t, preserving the
    legacy peak selection and area convention.
    """
    from .metrics import _area_ratio

    schema = EyeFlowOutputPaths.active()
    ROOT = schema.spectrograms_masked_metrics_root
    PHASE_PATH = schema.spectrograms_masked_time_aggregated_median_beat
    POOLED_PATH = schema.spectrograms_masked_time_aggregated_median
    pooled, phases = np.asarray(pooled_spectrum, dtype=float), np.asarray(phase_spectrum, dtype=float)
    if phases.ndim != 2 or pooled.ndim != 1 or phases.shape[1] != pooled.size:
        raise ValueError("Pooled and phase spectra must share the frequency-bin domain.")
    phase_samples = np.asarray(phase_attrs["phase_samples"])
    if phase_samples.shape != phases.shape[:1]:
        raise ValueError("phase_samples must match the phase spectrum.")
    windows = early_late_phase_windows(phase_samples, early_interval=early_interval, late_interval=late_interval)
    area_attrs = dict(base_metrics[f"{ROOT}/A"][1])
    spacing = float(area_attrs["frequency_bin_spacing_khz"])
    unit = area_attrs["unit"]
    phase_area = np.asarray(base_metrics[f"{ROOT}/A_t"][0], dtype=float)
    peaks = np.asarray([base_metrics[f"{ROOT}/index_max{n}"][0] for n in (1, 2)], dtype=float)
    valid_peaks = np.isfinite(peaks).all()
    percentages = tuple(percentages)
    if any(not np.isfinite(p) or not 0 < p <= 50 for p in percentages):
        raise ValueError("Vertical percentages must lie in (0, 50].")
    common = {
        key: value for key, value in area_attrs.items()
        if key.startswith(("frequency_", "fft_", "sampling_", "spectrogram_frequency_"))
    }
    common.update({
        "metric_generator": "pipelines.low_freq_spectra.shape_metrics",
        "vertical_metric_version": 2,
        "vertical_midpoint_bin_index": pooled.size // 2 - 0.5,
        "vertical_midpoint_definition": "boundary_between_frequency_halves_at_floor_bin_count_over_2_minus_0.5",
        "source_bin_count": pooled.size,
        "frequency_bin_spacing_khz": spacing,
        "frequency_bin_order": area_attrs.get("frequency_bin_order", "HD_export_order"),
        "area_coordinate": "unwrapped_frequency_khz_in_HD_bin_order",
        "integration_method": "trapezoidal_with_linearly_interpolated_bin_boundaries",
        "interval_endpoints": "both_included",
        "index_base": 0,
        "fixed_peak_indices": peaks,
        "peak_source_path": "/" + POOLED_PATH,
        "fixed_peak_index_source_paths": [f"/{ROOT}/index_max{i}" for i in (1, 2)],
        "peak_selection": "fixed_pooled_maximum_per_bin_domain_half_ignore_NaN_lowest_index_on_ties",
    })
    output = {}

    def put(name, value, attrs):
        path = f"{ROOT}/{name}"
        value = np.asarray(value, dtype=np.float64)
        if np.any(np.isfinite(value) & (np.abs(value) > np.finfo(np.float32).max)):
            raise ValueError(f"{path} cannot be stored as float32 without overflow.")
        output[path] = (value.astype(np.float32), attrs)

    for percent in percentages:
        fraction = percent / 100
        pooled_vertical = np.nan
        curve = np.full(phases.shape[0], np.nan)
        bounds = np.full((2, 2), np.nan)
        if valid_peaks:
            left, right = peaks
            midpoint = pooled.size // 2 - 0.5
            bounds = np.array([
                [left, left + (midpoint - left) * fraction],
                [right - (right - midpoint) * fraction, right],
            ])
            pooled_vertical = sum(interval_area(pooled, *b, spacing) for b in bounds)
            curve = sum(interval_area(phases, *b, spacing) for b in bounds)
        attrs = {
            **common,
            "vertical_fraction_per_side": fraction,
            "vertical_width_definition": "percentage/100_of_each_peak_to_frequency_domain_midpoint_distance",
            "integration_intervals_bin_indices": bounds,
            "nonfinite_value_policy": "nonfinite_samples_supporting_either_interval_return_NaN",
            "quantity": "frequency_interval_power_integral",
            "unit": unit,
            "dimDesc": [],
        }
        area_name, ratio_name, summary_name = f"Av_{percent}", f"SEv_{percent}", f"SEIv_{percent}"
        definition = f"sum of the inward {percent}% integrals of each peak-to-frequency-domain-midpoint span"
        put(area_name, pooled_vertical, {
            **attrs, "source_path": "/" + POOLED_PATH, "definition": definition,
        })
        curve_attrs = {
            **attrs, "source_path": "/" + PHASE_PATH,
            "definition": definition + " per phase with fixed pooled peak bounds",
            "phase_samples": phase_samples, "phase_endpoint": "excluded", "dimDesc": ["phase"],
            **{k: v for k, v in phase_attrs.items() if k.startswith(("cycle_boundary_", "beat_duration_"))},
        }
        put(area_name + "_t", curve, curve_attrs)
        ratio = _area_ratio(curve, phase_area)
        ratio_attrs = {
            **curve_attrs, "unit": "1", "quantity": "low_frequency_power_fraction",
            "definition": f"{area_name}_t divided by A_t per phase",
            "numerator_source_path": f"/{ROOT}/{area_name}_t",
            "denominator_source_path": f"/{ROOT}/A_t",
            "denominator_validity_rule": "finite_and_strictly_positive_A_t",
            "nonfinite_value_policy": "nonfinite_numerator_denominator_or_result_returns_NaN",
        }
        put(ratio_name, ratio, ratio_attrs)
        for values, source_name, base, source_attrs in (
            (curve, area_name + "_t", area_name, curve_attrs),
            (ratio, ratio_name, summary_name, ratio_attrs),
        ):
            summary_attrs = {k: v for k, v in source_attrs.items() if k not in {"phase_samples", "phase_endpoint"}}
            summary_attrs.update({
                "source_path": f"/{ROOT}/{source_name}", "dimDesc": [],
                "aggregation_method": "nanmedian", "aggregation_axes": ["phase"],
                "aggregation_weighting": "equal_weight_per_valid_phase_sample",
                "nonfinite_value_policy": "ignore_NaN_all_invalid_or_empty_window_returns_NaN",
            })
            if base == summary_name:
                put(base, nanmedian(values, dtype=np.float64), {
                    **summary_attrs,
                    "definition": f"median of valid phasewise {source_name} over the full cardiac cycle",
                    "phase_interval": [0., 1.], "phase_interval_endpoints": "source_endpoint_excluded",
                    "source_phase_sample_count": values.size,
                    "valid_sample_count": int(np.isfinite(values).sum()),
                    "invalid_sample_count": int((~np.isfinite(values)).sum()),
                })
            medians = phase_window_medians(values, windows)
            for window, (mask, interval) in windows.items():
                put(base + "_" + window, medians[window], {
                    **summary_attrs,
                    "definition": f"median of phasewise {source_name} in the {window} phase interval",
                    "phase_interval": interval, "phase_interval_endpoints": "both_included",
                    "source_phase_sample_count": int(mask.sum()),
                    "valid_sample_count": int(np.isfinite(values[mask]).sum()),
                    "invalid_sample_count": int((~np.isfinite(values[mask])).sum()),
                })
            put("delta_" + base, late_minus_early(medians["early"], medians["late"]), {
                **{k: v for k, v in summary_attrs.items() if not k.startswith("aggregation_")},
                "definition": f"{base}_late minus {base}_early",
                "calculation_method": "late_median_minus_early_median",
                "source_paths": [f"/{ROOT}/{base}_late", f"/{ROOT}/{base}_early"],
                "early_phase_interval": windows["early"][1],
                "late_phase_interval": windows["late"][1],
                "nonfinite_value_policy": "missing_early_or_late_summary_returns_NaN",
            })
    return output

