"""MATLAB reconstruction, signed DC preservation, and explicit harmonic semantics."""

import json
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.signal import butter, filtfilt, resample

from app_settings import (
    DEFAULT_BAND_LIMITED_HARMONICS_SETTING,
    AppSettingsStore,
    validate_band_limited_harmonics_setting,
)
from calculations.blood_flow_velocity.signal_analysis.per_beat.runner import (
    PerBeatAnalysisInput,
    run_per_beat_analysis,
)
from calculations.math import band_limit, butter_lowpass_filtfilt
from calculations.math.bandlimiting import DEFAULT_BAND_LIMITED_HARMONICS, validate_harmonics
from input_output.schema import EyeFlowOutputPaths
from input_output.schema.eyeflow_output import processing_path
from pipeline_engine.run_service import resolve_run_spec
from pipelines.signal_settings import band_limited_harmonics
from pipelines.velocity_analysis.outputs.per_beat import pack_velocity_per_beat_outputs


def test_application_and_calculation_defaults_agree_with_shipped_settings():
    from pathlib import Path

    defaults = json.loads((Path(__file__).resolve().parents[1] / "default_settings.json").read_text())
    assert defaults["band_limited_harmonics"] == DEFAULT_BAND_LIMITED_HARMONICS_SETTING
    assert DEFAULT_BAND_LIMITED_HARMONICS_SETTING == DEFAULT_BAND_LIMITED_HARMONICS == 12
    assert validate_harmonics(np.int64(7)) == 7
    assert validate_band_limited_harmonics_setting(np.int64(7)) == 7


@pytest.mark.parametrize("count", [63, 64, 128])
@pytest.mark.parametrize("harmonics", [0, 2, 12, 100])
def test_magnitude_matches_matlab_positive_spectrum_recipe(count, harmonics):
    rng = np.random.default_rng(24)
    values = (20 + rng.normal(size=(count, 3))).astype(np.float32)
    original = values.copy()
    spectrum = np.fft.fft(values, axis=0)
    stop = min(harmonics, count // 2) + 1
    retained = spectrum[:stop] * 2
    retained[0] = spectrum[0]
    expected = np.abs(np.fft.ifft(retained, n=count, axis=0))
    actual = band_limit(values, harmonics=harmonics, axis=0)
    np.testing.assert_allclose(actual, expected, rtol=1e-6, atol=1e-6)
    np.testing.assert_array_equal(values, original)
    assert actual.dtype == np.float32


@pytest.mark.parametrize("count", [63, 64])
def test_signed_reconstruction_preserves_sign_dc_and_nyquist(count):
    phase = np.arange(count) / count
    values = (-2 + 5 * np.cos(2 * np.pi * phase) + np.cos(2 * np.pi * 13 * phase)).astype(
        np.float32
    )
    expected = -2 + 5 * np.cos(2 * np.pi * phase)
    result = band_limit(values, reconstruction="signed")
    np.testing.assert_allclose(result, expected, atol=1e-6)
    np.testing.assert_allclose(result.mean(), values.mean(), atol=1e-6)
    assert result.min() < 0
    np.testing.assert_allclose(
        band_limit(values, harmonics=100, reconstruction="signed"), values, atol=1e-6
    )


@pytest.mark.parametrize("count", [-1, 1.5, True, np.bool_(True)])
def test_invalid_harmonic_counts_rejected(count):
    with pytest.raises(ValueError, match="harmonics"):
        band_limit(np.ones(32), harmonics=count)


def test_nonfinite_series_does_not_poison_other_batch_columns():
    values = np.ones((64, 2), dtype=np.float32)
    values[4, 0] = np.nan
    actual = band_limit(values, reconstruction="signed", axis=0)
    assert np.isnan(actual[:, 0]).all()
    np.testing.assert_allclose(actual[:, 1], 1)


def test_butterworth_uses_matlab_twelve_sample_padding_and_filters_short_valid_records():
    values = np.linspace(0, 5, 15, dtype=np.float32) ** 2
    b, a = butter(4, float(np.float32(15) / (np.float32(0.5) / np.float32(0.01))))
    expected = filtfilt(b, a, values, padlen=12).astype(np.float32)
    actual = butter_lowpass_filtfilt(
        values, dt_seconds=np.float32(0.01), lowpass_freq_hz=np.float32(15)
    )
    np.testing.assert_array_equal(actual, expected)
    assert not np.array_equal(actual, values)


def test_filtered_beats_subdivide_continuous_signal_then_follow_matlab_fft_recipe():
    phase = np.arange(257) / 128
    signal = (50 + 8 * np.sin(2 * np.pi * phase) + 3 * np.cos(2 * np.pi * 25 * phase)).astype(
        np.float32
    )
    segments = np.broadcast_to(signal, (2, 3, signal.size)).copy()
    result = run_per_beat_analysis(
        PerBeatAnalysisInput(
            arterial_velocity_signal=signal,
            venous_velocity_signal=signal,
            arterial_velocity_segments=segments,
            venous_velocity_segments=segments,
            cycle_boundary_indexes=np.array([0, 128, 256]),
            index_base=0,
            band_limited_harmonics=12,
            cardiac_cycle=None,
            dt_seconds=0.01,
        )
    )
    b, a = butter(4, float(np.float32(15) / (np.float32(0.5) / np.float32(0.01))))
    continuous = filtfilt(b, a, signal, padlen=12).astype(np.float32)
    expected = np.stack([resample(continuous[start : start + 129], 129)[:-1] for start in (0, 128)])
    actual = result.artery_filtered.signal
    np.testing.assert_allclose(actual.velocity_signal_per_beat, expected, rtol=1e-6, atol=1e-5)
    spectrum = np.fft.fft(expected, axis=1)[:, :13] * 2
    spectrum[:, 0] /= 2
    np.testing.assert_allclose(
        actual.velocity_signal_per_beat_band_limited,
        np.abs(np.fft.ifft(spectrum, n=128, axis=1)),
        rtol=1e-6,
        atol=1e-5,
    )
    np.testing.assert_allclose(
        result.artery_filtered.segments.velocity_signal_per_beat_per_segment[:, :, 0, 0],
        expected.T,
        rtol=1e-6,
        atol=1e-5,
    )
    outputs = pack_velocity_per_beat_outputs(result)
    schema = EyeFlowOutputPaths.active()
    assert schema.artery_per_beat_filtered.velocity_signal_band_limited in outputs
    assert not any("Safe" in path or path.endswith("/value") for path in outputs)
    for path in outputs:
        assert processing_path(path, "ProcessingBandRatio").startswith("ProcessingBandRatio/Velocity/PerBeat/")
    attrs = outputs[schema.artery_per_beat_filtered.velocity_signal_band_limited][1]
    assert attrs["band_limit_harmonics"] == 12
    assert attrs["band_limit_reconstruction"] == "matlab_magnitude"


def test_shared_setting_persists_and_import_rejects_invalid_count_atomically(tmp_path):
    store = AppSettingsStore(tmp_path / "settings.json")
    store.save_band_limited_harmonics(7)
    assert store.load_band_limited_harmonics() == 7
    imported = tmp_path / "invalid.json"
    imported.write_text(json.dumps({"band_limited_harmonics": -1}), encoding="utf-8")
    with pytest.raises(ValueError, match="harmonics"):
        store.import_file(imported)
    assert store.load_band_limited_harmonics() == 7
    assert (
        band_limited_harmonics(SimpleNamespace(attrs={"band_limited_harmonics": np.int64(7)})) == 7
    )
    assert (
        band_limited_harmonics(SimpleNamespace(attrs={"BandLimitedSignalHarmonicCount": 13})) == 12
    )


@pytest.mark.parametrize("count", [-1, 1.5, True])
def test_run_validation_rejects_invalid_harmonics_before_resolving_inputs(count):
    with pytest.raises(ValueError, match="harmonics"):
        resolve_run_spec(
            input_paths=[], target_names=[], pipelines=[], band_limited_harmonics=count
        )
