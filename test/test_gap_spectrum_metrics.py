"""Peak selection, missing-value contracts, and scalar gap-spectrum HDF5 outputs."""

from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from input_output.writers.h5_output import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
from pipelines.low_freq_spectra.metrics import (
    FLOOD_PERCENTAGES, PHASE_METRIC_OUTPUT_NAMES, PHASE_SUMMARY_OUTPUT_NAMES,
    SPECTRUM_METRIC_OUTPUT_NAMES, analyze_spectrum,
)
from pipelines.low_freq_spectra.gap_analysis import run_gap_analysis
from test_spectrogram_frequency import _write_frequency_table


def _spectrum():
    values = np.ones(32, dtype=np.float32)
    values[2], values[28], values[15], values[16] = 10, 6, 2, 4
    return values


def test_metrics_match_hand_calculated_khz_areas_and_leave_input_unchanged():
    spectrum = _spectrum()
    original = spectrum.copy()
    result = analyze_spectrum(spectrum, sampling_frequency_hz=32000)
    assert {name: result[name] for name in ("index_max1", "max1", "index_max2", "max2", "A1", "A0")} == {
        "index_max1": 2,
        "max1": 10,
        "index_max2": 28,
        "max2": 6,
        "A1": 37,
        "A0": 26,
    }
    np.testing.assert_array_equal(spectrum, original)


def test_khz_areas_use_sampling_rate_and_peak_distance_with_trapezoidal_endpoints():
    result = analyze_spectrum(_spectrum(), sampling_frequency_hz=8000)
    spacing = 0.25
    assert result["A0"] == 26 * spacing
    assert result["A1"] == 37 * spacing  # The endpoint samples contribute half their height.
    unit_spacing = analyze_spectrum(_spectrum(), sampling_frequency_hz=32000)
    for p in FLOOD_PERCENTAGES:
        assert result[f"A_{p}"] == pytest.approx(unit_spacing[f"A_{p}"] * spacing)


def test_peak_selection_uses_one_half_each_lowest_tie_index_and_includes_endpoints():
    spectrum = _spectrum()
    spectrum[0] = spectrum[2] = 100
    spectrum[31] = 6
    result = analyze_spectrum(spectrum, sampling_frequency_hz=32000)
    # Both globally highest samples are in the first half; still use one per half.
    assert result["index_max1"] == 0
    assert result["index_max2"] == 28
    spectrum[31] = 7
    assert analyze_spectrum(spectrum, sampling_frequency_hz=32000)["index_max2"] == 31
    spectrum[15], spectrum[16] = 200, 300
    result = analyze_spectrum(spectrum, sampling_frequency_hz=32000)
    assert (result["index_max1"], result["index_max2"]) == (15, 16)
    assert result["A1"] == 250
    assert result["A_50"] == 125


def test_runner_persists_scalar_metrics_from_the_raw_aggregate_spectrum(
    tmp_path,
):
    bin_count, valley_length, expected_length = 32, None, 8
    schema = EyeFlowOutputPaths.active()
    spectrum = _spectrum() if bin_count == 32 else np.ones(bin_count)
    with h5py.File(tmp_path / "output.h5", "w") as output:
        output.create_dataset(schema.sampling_frequency_hz, data=32000)
        _write_frequency_table(output, bin_count, 32000)
        source = output.create_dataset(
            schema.spectrograms_masked,
            data=np.broadcast_to(spectrum, (9, bin_count)),
        )
        source.attrs["unit"] = "power_unit"
        ctx = SimpleNamespace(output=SimpleNamespace(h5=PipelineH5Output(output)))
        moments = SimpleNamespace(
            cycle_boundary_indexes=np.array([0, 2, 8]),
            cycle_durations_seconds=np.array([0.2, 0.6]),
            cardiac_cycle_source="artery",
            continuous=lambda vessel, raw: np.zeros(9),
        )
        with patch("pipelines.low_freq_spectra.metrics.Logger.log_warning") as warning:
            metrics = run_gap_analysis(ctx, moments, valley_neighborhood_length=valley_length)
        ctx.output.h5.write_many(metrics)
        median = output[schema.spectrograms_masked_time_aggregated_median]
        assert median.name == "/ProcessingSpectrograms/Gaps/TimeAggregated/MedianPhaseBeat"
        np.testing.assert_array_equal(median[()], spectrum)
        group = output[schema.spectrograms_masked_metrics_root]
        expected_keys = {name for name in PHASE_SUMMARY_OUTPUT_NAMES.values() if name.startswith("SEI")} | set(SPECTRUM_METRIC_OUTPUT_NAMES.values()) | set(PHASE_METRIC_OUTPUT_NAMES.values()) | {
            key for name in PHASE_SUMMARY_OUTPUT_NAMES.values()
            for key in (f"{name}_early", f"{name}_late", f"delta_{name}")
        }
        assert group.name == "/ProcessingSpectrograms/Gaps/Metrics"
        # The selected shape descriptors have a separate schema/reference test.
        assert expected_keys <= set(group)
        assert not {"index_min1", "index_min2", "min1", "min2", "A2", "R1", "A2_t", "R1_t",
                    "A2_early", "A2_late", "delta_A2", "R1_early", "R1_late", "delta_R1"} & set(group)
        assert not {"h0", "h1", "hv", "rho_h", "rho_v"} & set(group)
        assert "ProcessingSpectrograms/Gaps/TimeAggregated/Metrics" not in output
        for name in sorted(expected_keys):
            dataset = group[name]
            assert dataset.shape == ((128,) if name in PHASE_METRIC_OUTPUT_NAMES.values() else ())
            assert dataset.dtype == np.float32
            assert not dataset.name.endswith("/value")
            assert list(dataset.attrs["dimDesc"]) == (["phase"] if name in PHASE_METRIC_OUTPUT_NAMES.values() else [])
            metric_name = name.removeprefix("delta_").split("_")[0]
            assert dataset.attrs["unit"] == (
                "1" if name.startswith("index_") or metric_name in {"SE", "SEI", "SEII"}
                else "power_unit*kHz" if metric_name in {"A0", "A"} else "power_unit"
            )
            assert dataset.attrs["sampling_frequency_hz"] == 32000
            assert dataset.attrs["frequency_bin_spacing_khz"] == pytest.approx(32 / bin_count)
            assert dataset.attrs["A_integration_method"] == "trapezoidal_khz_spacing"
            if metric_name in {"SE", "SEI", "SEII"}:
                assert "A1" not in dataset.attrs["definition"]
            assert dataset.attrs["valley_neighborhood_length"] == expected_length
            middle = bin_count // 2
            np.testing.assert_array_equal(
                dataset.attrs["valley_neighborhood_indices"],
                np.arange(middle - expected_length // 2, middle + expected_length // 2),
            )
            assert dataset.attrs["valley_aggregation_method"] == "nanmedian"
            if name.startswith("index_"):
                assert dataset.attrs["index_base"] == 0
            expected = analyze_spectrum(
                spectrum, sampling_frequency_hz=32000, valley_neighborhood_length=valley_length,
                fft_bin_count=4 * bin_count, binning_factor=4,
            )
            calculation_name = next(
                (internal for internal, external in SPECTRUM_METRIC_OUTPUT_NAMES.items() if external == name),
                None,
            )
            if calculation_name is not None:
                assert dataset.attrs["source_path"] == median.name
                assert dataset[()] == np.float32(expected[calculation_name])
            elif name.startswith("delta_"):
                assert dataset[()] == 0  # A constant source has identical early/late medians.
            else:
                phase_calculation_name = next(
                    internal for internal in PHASE_METRIC_OUTPUT_NAMES
                    if name in {PHASE_METRIC_OUTPUT_NAMES[internal],
                                PHASE_SUMMARY_OUTPUT_NAMES[internal],
                                f"{PHASE_SUMMARY_OUTPUT_NAMES[internal]}_early",
                                f"{PHASE_SUMMARY_OUTPUT_NAMES[internal]}_late"}
                )
                reference = "A0" if phase_calculation_name == "SE" or phase_calculation_name.startswith("SEII_") else phase_calculation_name.replace("SE", "A", 1)
                signal = phase_calculation_name.replace("SEII", "A", 1) if phase_calculation_name.startswith("SEII_") else "A1"
                expected_value = (
                    (expected[signal] - expected[reference]) / expected[reference]
                    if phase_calculation_name.startswith("SE") else expected[phase_calculation_name]
                )
                np.testing.assert_allclose(dataset[()], expected_value)
        assert "ProcessingBandRatio/Spectrograms" not in output
        warning.assert_not_called()
