"""Reference FFT averaging and HD spectral calibration contracts."""

import json
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from calculations.math.frequency import frequency_coordinates_khz
from input_output.readers.h5 import RawH5SourceReader
from input_output.writers.h5_output import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths, HolodopplerSource
from pipelines.low_freq_spectra.frequency import write_frequency_bin_table


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


def test_unbinned_coordinates_follow_standard_even_and_odd_fft_order():
    bin_count, sampling_frequency = 6, 12
    expected_hz = [0, 2, 4, -6, -4, -2]
    np.testing.assert_allclose(
        frequency_coordinates_khz(bin_count, sampling_frequency),
        np.asarray(expected_hz) / 1000,
    )


def test_frequency_table_matches_current_hd_selection_before_averaging(tmp_path):
    start, end = 0, 0.5
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


def test_unknown_or_inconsistent_export_geometry_skips_table_without_losing_power(tmp_path):
    batch_size, bin_count = 256, 63
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "hd.h5", "w") as hd, h5py.File(tmp_path / "out.h5", "w") as out:
        hd.create_dataset("HD_parameters", data=json.dumps({"batch_size": batch_size}))
        source = HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD"))
        out.create_dataset(schema.sampling_frequency_hz, data=37037)
        out.create_dataset(schema.spectrograms_spatial_average, data=np.ones((2, bin_count)))
        write_frequency_bin_table(SimpleNamespace(h5=PipelineH5Output(out)), bin_count, hd=source)
        assert schema.spectrograms_frequency_bin_table not in out
        assert schema.spectrograms_spatial_average in out
