"""Frequency-band quantitative velocity estimation contracts."""

from __future__ import annotations

import json
import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from pipelines.retinal_velocity.estimation import (
    DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
    DEFAULT_LASER_WAVELENGTH_METERS,
    DEFAULT_NUMERICAL_APERTURE,
    _safe_band_ratio,
    estimate_retinal_velocity,
)
from input_output.schema import DopplerViewSource, HolodopplerSource
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine.context import RawH5SourceReader
from pipelines.vessel_inputs import load_retinal_source_data
from pipelines.retinal_velocity.outputs import pack_retinal_velocity_outputs
from pipelines.retinal_velocity.semantics import resolve_velocity_semantics
from pipelines.waveform_velocity.outputs.continuous import (
    pack_continuous_velocity_outputs,
)


LF_PATH = "band_0_3000_9000"
HF_PATH = "band_1_9000_18000"


def test_safe_band_ratio_maps_exact_zero_lf_to_zero_without_epsilon() -> None:
    high = np.asarray([[[2.0, 7.0, 0.0]]], dtype=np.float32)
    low = np.asarray([[[1.0, 0.0, 0.5]]], dtype=np.float32)

    ratio = _safe_band_ratio(high, low, frame_slice=slice(0, 1))

    np.testing.assert_array_equal(
        ratio,
        np.asarray([[[2.0, 0.0, 0.0]]], dtype=np.float32),
    )
    assert np.all(np.isfinite(ratio))


def test_frequency_band_estimator_converts_ratio_frequency_to_mm_per_second() -> None:
    shape = (2, 16, 16)
    low = np.ones(shape, dtype=np.float32)
    high = np.ones(shape, dtype=np.float32)
    artery = np.zeros(shape[1:], dtype=bool)
    vein = np.zeros(shape[1:], dtype=bool)
    artery[8, 10] = True
    vein[10, 8] = True
    high[:, artery | vein] = 4.0

    with tempfile.TemporaryDirectory() as tmp_dir:
        with h5py.File(Path(tmp_dir) / "scratch.h5", "w") as scratch:
            result = estimate_retinal_velocity(
                band_lf=low,
                band_hf=high,
                velocity_estimation_method="frequency_bands",
                artery_mask=artery,
                vein_mask=vein,
                optic_disc_center=(7.5, 7.5),
                local_background_dist=1,
                scratch_h5=scratch,
            )
            velocity = np.asarray(result.maps.velocity)

    # At the default 1 Hz-per-ratio calibration, the vessel and inpainted
    # neighbourhood are 4 Hz and 1 Hz. Their signed RMS difference is then
    # converted to mm/s by the same wavelength/NA law as moment mode.
    expected_delta_hz = np.float32(np.sqrt(4.0**2 - 1.0**2))
    expected = np.float32(
        1e3
        * DEFAULT_LASER_WAVELENGTH_METERS
        * expected_delta_hz
        / DEFAULT_NUMERICAL_APERTURE
    )
    np.testing.assert_allclose(velocity[:, artery], expected, rtol=1e-5)
    np.testing.assert_allclose(velocity[:, vein], expected, rtol=1e-5)
    np.testing.assert_array_equal(result.maps.moment0_average, np.ones(shape[1:]))
    np.testing.assert_allclose(result.maps.frms_average[artery | vein], 4.0)
    assert result.provenance["velocity_estimation_method"] == "frequency_bands"
    assert result.provenance["velocity_quantity"] == "physical_velocity"
    assert result.provenance["velocity_unit"] == "mm/s"
    assert (
        result.provenance["band_ratio_frequency_scale_hz"]
        == DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ
    )
    assert result.provenance["band_lf_source_path"] == f"/{LF_PATH}"
    assert result.provenance["band_hf_source_path"] == f"/{HF_PATH}"


def test_frequency_band_frequency_scale_is_explicit_and_changes_velocity() -> None:
    shape = (1, 12, 12)
    low = np.ones(shape, dtype=np.float32)
    high = np.ones(shape, dtype=np.float32)
    artery = np.zeros(shape[1:], dtype=bool)
    vein = np.zeros_like(artery)
    artery[6, 8] = True
    high[:, artery] = 4.0

    with h5py.File(
        "scratch.h5", "w", driver="core", backing_store=False
    ) as scratch:
        result = estimate_retinal_velocity(
            band_lf=low,
            band_hf=high,
            velocity_estimation_method="frequency_bands",
            band_ratio_frequency_scale_hz=2.0,
            artery_mask=artery,
            vein_mask=vein,
            optic_disc_center=(5.5, 5.5),
            local_background_dist=1,
            scratch_h5=scratch,
        )
        velocity = np.asarray(result.maps.velocity).copy()

    expected_delta_hz = np.float32(2.0 * np.sqrt(4.0**2 - 1.0**2))
    expected_velocity = np.float32(
        1e3
        * DEFAULT_LASER_WAVELENGTH_METERS
        * expected_delta_hz
        / DEFAULT_NUMERICAL_APERTURE
    )
    np.testing.assert_allclose(
        velocity[:, artery],
        expected_velocity,
        rtol=1e-5,
    )
    np.testing.assert_allclose(result.maps.frms_average[artery], 8.0)
    assert result.provenance["band_ratio_frequency_scale_hz"] == 2.0


@pytest.mark.parametrize("scale", [0.0, -1.0, np.nan, np.inf])
def test_frequency_band_estimator_rejects_invalid_frequency_scale(
    scale: float,
) -> None:
    values = np.ones((1, 4, 4), dtype=np.float32)
    mask = np.zeros((4, 4), dtype=bool)
    with h5py.File(
        "scratch.h5", "w", driver="core", backing_store=False
    ) as scratch:
        with pytest.raises(ValueError, match="band_ratio_frequency_scale_hz"):
            estimate_retinal_velocity(
                band_lf=values,
                band_hf=values,
                velocity_estimation_method="frequency_bands",
                band_ratio_frequency_scale_hz=scale,
                artery_mask=mask,
                vein_mask=mask,
                optic_disc_center=(1.5, 1.5),
                local_background_dist=1,
                scratch_h5=scratch,
            )


def test_frequency_band_estimator_reports_zero_and_near_zero_lf_counts() -> None:
    low = np.ones((1, 8, 8), dtype=np.float32)
    low[0, 0, 0] = 0.0
    low[0, 0, 1] = 1e-8
    high = np.ones_like(low)
    mask = np.zeros((8, 8), dtype=bool)

    with h5py.File(
        "scratch.h5", "w", driver="core", backing_store=False
    ) as scratch:
        result = estimate_retinal_velocity(
            band_lf=low,
            band_hf=high,
            velocity_estimation_method="frequency_bands",
            artery_mask=mask,
            vein_mask=mask,
            optic_disc_center=(3.5, 3.5),
            local_background_dist=1,
            scratch_h5=scratch,
            retain_velocity_video=False,
        )

    assert result.provenance["band_lf_zero_sample_count"] == 1
    assert result.provenance["band_lf_near_zero_sample_count"] == 1
    assert result.provenance["band_lf_vessel_zero_sample_count"] == 0
    assert result.provenance["band_lf_vessel_near_zero_sample_count"] == 0
    assert result.provenance["band_lf_neighborhood_zero_sample_count"] == 1
    assert result.provenance["band_lf_neighborhood_near_zero_sample_count"] == 1
    assert result.provenance["band_lf_low_relative_threshold"] == 1e-6


@pytest.mark.parametrize(
    ("path", "value", "message"),
    [
        (LF_PATH, -1.0, "negative values"),
        (HF_PATH, np.nan, "non-finite values"),
    ],
)
def test_frequency_band_estimator_rejects_invalid_psd_values(
    path: str,
    value: float,
    message: str,
) -> None:
    low = np.ones((1, 8, 8), dtype=np.float32)
    high = np.ones_like(low)
    (low if path == LF_PATH else high)[0, 0, 0] = value
    mask = np.zeros((8, 8), dtype=bool)

    with tempfile.TemporaryDirectory() as tmp_dir:
        with h5py.File(Path(tmp_dir) / "scratch.h5", "w") as scratch:
            with pytest.raises(ValueError, match=message):
                estimate_retinal_velocity(
                    band_lf=low,
                    band_hf=high,
                    velocity_estimation_method="frequency_bands",
                    artery_mask=mask,
                    vein_mask=mask,
                    optic_disc_center=(3.5, 3.5),
                    local_background_dist=1,
                    scratch_h5=scratch,
                    retain_velocity_video=False,
                )


def test_frequency_band_source_requires_exact_paths_and_does_not_require_moments() -> None:
    with tempfile.TemporaryDirectory() as tmp_dir:
        root = Path(tmp_dir)
        hd_path = root / "sample_HD.h5"
        dv_path = root / "sample_DV.h5"
        with h5py.File(hd_path, "w") as hd:
            hd.create_dataset(LF_PATH, data=np.ones((2, 4, 6), dtype=np.float32))
            hd.create_dataset(HF_PATH, data=np.full((2, 4, 6), 2.0, dtype=np.float32))
            hd.create_dataset("sampling_freq", data=np.float32(100.0))
            hd.create_dataset("batch_stride", data=np.float32(10.0))
            hd.create_dataset(
                "HD_parameters",
                data=json.dumps({"pixel_pitch": [20e-6, 20e-6]}),
            )
        with h5py.File(dv_path, "w") as dv:
            retina = dv.create_group("segmentation/Retina")
            retina.create_dataset("artery_mask", data=np.zeros((4, 6), dtype=bool))
            retina.create_dataset("vein_mask", data=np.zeros((4, 6), dtype=bool))

        with h5py.File(hd_path, "r") as hd, h5py.File(dv_path, "r") as dv:
            source = load_retinal_source_data(
                HolodopplerSource(RawH5SourceReader(h5file=hd, label="HD")),
                DopplerViewSource(RawH5SourceReader(h5file=dv, label="DV")),
                velocity_estimation_method="frequency_bands",
            )
            assert source.image_maps.moment0 is None
            assert source.image_maps.moment2 is None
            assert source.image_maps.band_lf.name == f"/{LF_PATH}"
            assert source.image_maps.band_hf.name == f"/{HF_PATH}"


def test_frequency_band_source_reports_missing_exact_dataset_path() -> None:
    with tempfile.TemporaryDirectory() as tmp_dir:
        hd_path = Path(tmp_dir) / "sample_HD.h5"
        with h5py.File(hd_path, "w") as hd:
            hd.create_dataset(LF_PATH, data=np.ones((1, 2, 2), dtype=np.float32))
            hd.create_dataset("band_1_other", data=np.ones((1, 2, 2), dtype=np.float32))
        with h5py.File(hd_path, "r") as hd:
            source = HolodopplerSource(
                RawH5SourceReader(h5file=hd, label="HD")
            )
            with pytest.raises(KeyError, match=f"/{HF_PATH}"):
                source.frequency_band_datasets()


@pytest.mark.parametrize(
    ("available_paths", "missing_paths"),
    [
        ((), (f"/{LF_PATH}", f"/{HF_PATH}")),
        ((LF_PATH,), (f"/{HF_PATH}",)),
        ((HF_PATH,), (f"/{LF_PATH}",)),
    ],
)
def test_frequency_band_source_lists_every_missing_exact_path(
    available_paths: tuple[str, ...],
    missing_paths: tuple[str, ...],
) -> None:
    with tempfile.TemporaryDirectory() as tmp_dir:
        hd_path = Path(tmp_dir) / "missing_bands_HD.h5"
        with h5py.File(hd_path, "w") as hd:
            for path in available_paths:
                hd.create_dataset(path, data=np.ones((1, 2, 2), dtype=np.float32))
        with h5py.File(hd_path, "r") as hd:
            source = HolodopplerSource(
                RawH5SourceReader(h5file=hd, label="HD")
            )
            with pytest.raises(KeyError) as error:
                source.frequency_band_datasets()

    message = str(error.value)
    assert "frequency_bands" in message
    assert "missing_bands_HD.h5" in message
    for path in missing_paths:
        assert path in message


@pytest.mark.parametrize(
    ("low", "high", "error_type", "message"),
    [
        (
            np.ones((2, 2), dtype=np.float32),
            np.ones((2, 2, 2), dtype=np.float32),
            ValueError,
            "3-D",
        ),
        (
            np.full((1, 2, 2), b"x", dtype="S1"),
            np.ones((1, 2, 2), dtype=np.float32),
            TypeError,
            "numeric",
        ),
        (
            np.ones((1, 2, 2), dtype=np.float32),
            np.ones((2, 2, 2), dtype=np.float32),
            ValueError,
            "identical",
        ),
    ],
)
def test_frequency_band_source_validates_rank_dtype_and_shape(
    low: np.ndarray,
    high: np.ndarray,
    error_type: type[Exception],
    message: str,
) -> None:
    with tempfile.TemporaryDirectory() as tmp_dir:
        hd_path = Path(tmp_dir) / "invalid_bands_HD.h5"
        with h5py.File(hd_path, "w") as hd:
            hd.create_dataset(LF_PATH, data=low)
            hd.create_dataset(HF_PATH, data=high)
        with h5py.File(hd_path, "r") as hd:
            source = HolodopplerSource(
                RawH5SourceReader(h5file=hd, label="HD")
            )
            with pytest.raises(error_type, match=message):
                source.frequency_band_datasets()


def test_frequency_band_estimator_rejects_spatial_mask_mismatch() -> None:
    low = np.ones((1, 4, 5), dtype=np.float32)
    high = np.ones_like(low)
    wrong_mask = np.zeros((5, 4), dtype=bool)
    with h5py.File("scratch.h5", "w", driver="core", backing_store=False) as scratch:
        with pytest.raises(ValueError, match="spatial shape"):
            estimate_retinal_velocity(
                band_lf=low,
                band_hf=high,
                velocity_estimation_method="frequency_bands",
                artery_mask=wrong_mask,
                vein_mask=wrong_mask,
                optic_disc_center=(2.0, 2.0),
                local_background_dist=1,
                scratch_h5=scratch,
            )


def test_frequency_band_estimator_is_independent_of_frame_chunk_size() -> None:
    rng = np.random.default_rng(42)
    shape = (9, 14, 12)
    low = (0.5 + rng.random(shape)).astype(np.float32)
    high = (low * (0.5 + 3.0 * rng.random(shape))).astype(np.float32)
    artery = np.zeros(shape[1:], dtype=bool)
    vein = np.zeros_like(artery)
    artery[5:7, 4:6] = True
    vein[8:10, 7:9] = True

    results: list[dict[str, np.ndarray]] = []
    for chunk_size in (1, 4, 32):
        video = np.empty(shape, dtype=np.float32)
        with (
            patch(
                "pipelines.retinal_velocity.estimation."
                "SCRATCH_FRAME_CHUNK_SIZE",
                chunk_size,
            ),
            h5py.File(
                "scratch.h5", "w", driver="core", backing_store=False
            ) as scratch,
        ):
            result = estimate_retinal_velocity(
                band_lf=low,
                band_hf=high,
                velocity_estimation_method="frequency_bands",
                artery_mask=artery,
                vein_mask=vein,
                optic_disc_center=(5.5, 6.5),
                local_background_dist=1,
                scratch_h5=scratch,
                velocity_video_output=video,
            )
            results.append(
                {
                    "velocity_map": np.asarray(result.maps.velocity).copy(),
                    "moment0_avg": result.maps.moment0_average.copy(),
                    "velocity_map_avg": result.maps.velocity_average.copy(),
                    "fRMS_avg": result.maps.frms_average.copy(),
                    "fRMS_bkg_avg": result.maps.frms_background_average.copy(),
                    "deltafRMS_avg": result.maps.delta_frms_average.copy(),
                    "retinal_artery_velocity_signal": result.artery.velocity.copy(),
                    "retinal_vein_velocity_signal": result.vein.velocity.copy(),
                }
            )

    for actual in results[1:]:
        for key, expected in results[0].items():
            np.testing.assert_array_equal(actual[key], expected)


def test_frequency_band_output_units_and_human_label_are_physical() -> None:
    schema = EyeFlowOutputPaths.active()
    values = np.asarray([1.0, 2.0], dtype=np.float32)
    provenance = {
        "velocity_estimation_method": "frequency_bands",
        "velocity_quantity": "physical_velocity",
        "velocity_unit": "mm/s",
        "band_ratio_frequency_scale_hz": 1.0,
        "band_ratio_calibration_model": "linear_origin",
        "band_ratio_calibration_source": "eyeflow_setting",
        "band_ratio_calibration_version": "1",
        "laser_wavelength_m": 8.52e-7,
        "numerical_aperture": 0.124,
    }
    analysis = SimpleNamespace(
        provenance=provenance,
        continuous=lambda _vessel, raw=False: values,
    )

    outputs = pack_continuous_velocity_outputs(analysis)
    semantics = resolve_velocity_semantics(analysis)

    assert outputs[schema.analysis.retinal_artery_velocity_signal][1]["unit"] == "mm/s"
    assert outputs[schema.analysis.retinal_vein_velocity_signal][1]["unit"] == "mm/s"
    assert (
        outputs[schema.analysis.retinal_artery_velocity_signal][1][
            "band_ratio_frequency_scale_hz"
        ]
        == 1.0
    )
    assert semantics.axis_label == "Velocity (mm/s)"


def test_frequency_maps_are_persisted_in_hz_with_calibration_provenance() -> None:
    schema = EyeFlowOutputPaths.active()
    provenance = {
        "velocity_estimation_method": "frequency_bands",
        "velocity_quantity": "physical_velocity",
        "velocity_unit": "mm/s",
        "band_ratio_frequency_scale_hz": 2.0,
        "band_ratio_calibration_model": "linear_origin",
        "band_ratio_calibration_source": "eyeflow_setting",
        "band_ratio_calibration_version": "1",
        "laser_wavelength_m": 8.52e-7,
        "numerical_aperture": 0.124,
    }
    spectral = SimpleNamespace(
        fundamental_hz=1.0,
        heart_rate_bpm=60.0,
        heart_rate_ste_bpm=0.0,
        period_seconds=1.0,
    )
    analysis = SimpleNamespace(
        provenance=provenance,
        maps=SimpleNamespace(
            velocity_average=np.ones((2, 2), dtype=np.float32),
            frms_average=np.ones((2, 2), dtype=np.float32),
            frms_background_average=np.ones((2, 2), dtype=np.float32),
            delta_frms_average=np.zeros((2, 2), dtype=np.float32),
        ),
        cardiac_cycle=SimpleNamespace(spectral=spectral),
        cycle_boundary_indexes=np.asarray([0, 1], dtype=np.int32),
        cycle_durations_seconds=np.asarray([1.0], dtype=np.float32),
    )

    outputs = pack_retinal_velocity_outputs(analysis)
    attrs = outputs[schema.analysis.fRMS_avg][1]

    assert attrs["unit"] == "Hz"
    assert attrs["quantity"] == "rms_frequency"
    assert attrs["velocity_estimation_method"] == "frequency_bands"
    assert attrs["band_ratio_frequency_scale_hz"] == 2.0


def test_dimensionless_velocity_metadata_is_rejected() -> None:
    with pytest.raises(ValueError, match="physical_velocity"):
        resolve_velocity_semantics(
            {
                "velocity_estimation_method": "frequency_bands",
                "velocity_quantity": "relative_velocity_index",
                "velocity_unit": "1",
            }
        )
