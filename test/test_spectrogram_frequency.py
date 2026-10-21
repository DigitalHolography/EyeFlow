"""Reference FFT averaging and HD spectral calibration contracts."""

import json
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from calculations.math.frequency import frequency_bin_spacing_khz, frequency_coordinates_khz
from input_output.h5_access import PipelineH5Output, RawH5SourceReader
from input_output.schema import EyeFlowOutputPaths, HolodopplerSource
from pipelines.low_freq_spectra.frequency import write_frequency_bin_table
from pipelines.low_freq_spectra.metrics import pack_spectrum_metrics


def _write_frequency_table(h5, bin_count, sampling_frequency_hz):
    """Supply calibrated coordinates to tests of consumers, including synthetic domains."""
    output = SimpleNamespace(h5=PipelineH5Output(h5))
    hd = SimpleNamespace(spectrogram_frequency_parameters=lambda count: {
        "fft_bin_count": 4 * count, "frequency_binning_factor": 4,
        "frequency_selection": "full_FFT_domain",
    })
    schema = EyeFlowOutputPaths.active()
    if schema.sampling_frequency_hz not in h5:
        h5.create_dataset(schema.sampling_frequency_hz, data=sampling_frequency_hz)
    write_frequency_bin_table(output, bin_count, hd=hd)


@pytest.mark.parametrize("bin_count, sampling_frequency, expected_hz", [
    (1, 10, [0]),
    (5, 10, [0, 2, 4, -4, -2]),
    (6, 12, [0, 2, 4, -6, -4, -2]),
])
def test_unbinned_coordinates_follow_standard_even_and_odd_fft_order(bin_count, sampling_frequency, expected_hz):
    np.testing.assert_allclose(
        frequency_coordinates_khz(bin_count, sampling_frequency),
        np.asarray(expected_hz) / 1000,
    )


@pytest.mark.parametrize("fft_size", [128, 256, 512])
def test_coordinates_match_full_fft_four_sample_averages(fft_size):
    raw = np.fft.fftfreq(fft_size, d=1 / 37037)
    expected = raw.reshape(-1, 4).mean(axis=1) / 1000
    actual = frequency_coordinates_khz(
        fft_size // 4, 37037, fft_bin_count=fft_size, binning_factor=4,
    )
    np.testing.assert_allclose(actual, expected)
    assert actual[0] == pytest.approx(1.5 * 37037 / fft_size / 1000)
    assert actual[len(actual) // 2] > -37037 / 2000
    assert actual[-1] == pytest.approx(-2.5 * 37037 / fft_size / 1000)
    for half in np.split(actual, 2):
        np.testing.assert_allclose(np.diff(half), 4 * 37037 / fft_size / 1000)


def test_first_bin_averages_dc_and_does_not_assume_dc_is_the_center():
    actual = frequency_coordinates_khz(64, 36864, fft_bin_count=256, binning_factor=4)
    assert actual[0] == pytest.approx(0.216)  # mean(0, 144, 288, 432) Hz
    assert frequency_bin_spacing_khz(64, 36864, fft_bin_count=256, binning_factor=4) == 0.576


@pytest.mark.parametrize("fft_size, bin_count", [(256, 64), (512, 128)])
def test_hd_batch_size_comes_from_hdf5_before_config(tmp_path, fft_size, bin_count):
    with h5py.File(tmp_path / "hd.h5", "w") as hd:
        hd.create_dataset("HD_parameters", data=json.dumps({
            "batch_size": fft_size, "sampling_freq": 37037, "spectrogram_frequency_binning": 4,
            "registration_ref_batch_size": 1024,
        }))
        source = HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD"), {
            "batch_size": 128,
        })
        parameters = source.spectrogram_frequency_parameters(bin_count)
        assert parameters["fft_bin_count"] == fft_size
        assert parameters["frequency_binning_factor"] == 4
        assert parameters["fft_bin_count_source_path"] == "/HD_parameters"
        assert parameters["spectrogram_frequency_start"] == 0
        assert parameters["spectrogram_frequency_end"] == 0.5
        assert parameters["frequency_selection_source"] == "full_spectrum_fallback"


def test_hd_spectral_parameters_fall_back_to_config(tmp_path):
    with h5py.File(tmp_path / "hd.h5", "w") as hd:
        source = HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD"), {"batch_size": 256})
        parameters = source.spectrogram_frequency_parameters(64)
        assert parameters["fft_bin_count"] == 256
        assert parameters["fft_bin_count_source_path"] == "parameters.json"


@pytest.mark.parametrize("start, end", [(0, 0.5), (0.25, 0.25), (0, 0.25), (0.125, 0.375)])
def test_frequency_table_matches_current_hd_selection_before_averaging(tmp_path, start, end):
    schema = EyeFlowOutputPaths.active()
    fft_size, factor, fs = 256, 4, 37037
    raw = np.fft.fftfreq(fft_size, d=1 / fs)
    selected = np.r_[raw[:int(fft_size * end)], raw[int(fft_size * start) + fft_size // 2:]]
    expected = selected.reshape(-1, factor).mean(axis=1) / 1000
    with h5py.File(tmp_path / "hd.h5", "w") as hd, h5py.File(tmp_path / "out.h5", "w") as out:
        hd.create_dataset("HD_parameters", data=json.dumps({
            "batch_size": fft_size, "spectrogram_frequency_binning": factor,
            "spectrogram_frequency_start": start, "spectrogram_frequency_end": end,
        }))
        source = HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD"))
        out.create_dataset(schema.sampling_frequency_hz, data=fs)
        write_frequency_bin_table(SimpleNamespace(h5=PipelineH5Output(out)), expected.size, hd=source)
        table = out[schema.spectrograms_frequency_bin_table]
        np.testing.assert_allclose(table[:], expected, rtol=1e-7)
        assert table.attrs["frequency_selection_source"] == "/HD_parameters"
        assert table.attrs["spectrogram_frequency_start"] == start
        assert table.attrs["spectrogram_frequency_end"] == end
        assert table.attrs["frequency_bin_spacing_khz"] == pytest.approx(0.578703125)
        assert table.attrs["bin_domain_split_index"] == np.count_nonzero(expected >= 0)
        if start == 0 and end == 0.5:
            assert table.shape == (64,)
            assert table[0] == pytest.approx(0.217013671875)
            assert table[5] == pytest.approx(3.110529296875)


@pytest.mark.parametrize("bounds", [
    {}, {"spectrogram_frequency_start": 0.25}, {"spectrogram_frequency_end": 0.25},
    {"spectrogram_frequency_start": None, "spectrogram_frequency_end": 0.25},
    {"spectrogram_frequency_start": "unknown", "spectrogram_frequency_end": 0.25},
    {"spectrogram_frequency_start": 0, "spectrogram_frequency_end": float("nan")},
])
def test_unresolved_hd_bounds_use_full_spectrum_even_if_config_has_a_crop(tmp_path, bounds):
    with h5py.File(tmp_path / "hd.h5", "w") as hd:
        hd.create_dataset("HD_parameters", data=json.dumps({"batch_size": 256, **bounds}))
        source = HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD"), {
            "spectrogram_frequency_start": 0.25, "spectrogram_frequency_end": 0.25,
        })
        parameters = source.spectrogram_frequency_parameters(64)
        assert parameters["spectrogram_frequency_start"] == 0
        assert parameters["spectrogram_frequency_end"] == 0.5
        assert parameters["frequency_selection_source"] == "full_spectrum_fallback"
        assert parameters["frequency_binning_factor"] == 4


def test_recorded_binning_must_match_the_exported_count(tmp_path):
    with h5py.File(tmp_path / "hd.h5", "w") as hd:
        hd.create_dataset("HD_parameters", data=json.dumps({
            "batch_size": 256, "spectrogram_frequency_binning": 8,
        }))
        source = HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD"))
        with pytest.raises(ValueError, match="spectrogram_frequency_binning"):
            source.spectrogram_frequency_parameters(64)


@pytest.mark.parametrize("batch_size, bin_count", [(None, 32), (0, 32), (256.5, 32),
                                                        (True, 32), (float("nan"), 32), (256, 63)])
def test_unknown_or_inconsistent_export_geometry_skips_table_without_losing_power(tmp_path, batch_size, bin_count):
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "hd.h5", "w") as hd, h5py.File(tmp_path / "out.h5", "w") as out:
        hd.create_dataset("HD_parameters", data=json.dumps({"batch_size": batch_size}))
        source = HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD"))
        out.create_dataset(schema.sampling_frequency_hz, data=37037)
        out.create_dataset(schema.spectrograms_spatial_average, data=np.ones((2, bin_count)))
        write_frequency_bin_table(SimpleNamespace(h5=PipelineH5Output(out)), bin_count, hd=source)
        assert schema.spectrograms_frequency_bin_table not in out
        assert schema.spectrograms_spatial_average in out


@pytest.mark.parametrize("fft_size, start, end, area_scale", [
    (128, 0, 0.5, 1), (256, 0.25, 0.25, 0.5),
])
def test_calibration_propagates_to_pooled_phase_and_window_areas(fft_size, start, end, area_scale):
    spectrum = np.ones(32)
    spectrum[2], spectrum[28] = 10, 6
    phase_spectrum = np.stack((spectrum, spectrum * 2, spectrum * 3))
    options = dict(sampling_frequency_hz=37037, source_attrs={"unit": "power_unit"},
                   phase_spectrum=phase_spectrum,
                   phase_attrs={"unit": "power_unit", "phase_samples": np.array([0., 1 / 3, 2 / 3])})
    old = pack_spectrum_metrics(spectrum, **options)
    corrected = pack_spectrum_metrics(spectrum, **options, frequency_calibration_attrs={
        "fft_bin_count": fft_size, "frequency_binning_factor": 4,
        "spectrogram_frequency_start": start, "spectrogram_frequency_end": end,
    })
    for path, (actual, attrs) in corrected.items():
        name = path.rsplit("/", 1)[-1].removeprefix("delta_")
        area = name.split("_")[0] in {"A0", "A1", "A2"}
        np.testing.assert_allclose(actual, old[path][0] * (area_scale if area else 1))
        assert attrs["frequency_bin_spacing_khz"] == pytest.approx(4 * 37037 / (fft_size * 1000))
        assert attrs["unit"] == ("power_unit*kHz" if area else old[path][1]["unit"])
