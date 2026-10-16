"""Independent timing/downstream products and workflow failure isolation."""

import json
from unittest.mock import patch

import h5py
import numpy as np
import pytest
from test_frequency_band_pipeline_integration import _write_frequency_band_run

from input_output.reports.pdf_report import _extract_parameters_from_h5
from input_output.schema import EyeFlowOutputPaths
from input_output.schema.eyeflow_output import processing_path
from pipeline_engine.run_service import execute_run, resolve_run_spec
from pipelines import load_pipeline_catalog
from pipelines.velocity.runner import (
    detect_velocity_cardiac_cycles,
    estimate_retinal_velocity,
    pack_velocity_outputs,
)
from pipelines.low_freq_spectra.spectrograms import preprocess_spectrograms
from pipelines.low_freq_spectra.per_beat import gap_spectrograms_per_beat


def _dual_input(root):
    holo = _write_frequency_band_run(root)
    hd_path = root / holo.stem / f"{holo.stem}_HD/h5/{holo.stem}_HD_output.h5"
    with h5py.File(hd_path, "a") as hd:
        low = hd["band_0_3000_9000"][:]
        high = hd["band_1_9000_18000"][:]
        # Different rhythms expose accidental timing reuse (bands remain 1 Hz).
        time = np.arange(high.shape[0], dtype=np.float32) / np.float32(20.0)
        moment_frequency = np.full_like(high, 300.0)
        moment_frequency[:, high[0] > 1.0] = (
            1200.0 + 240.0 * np.sin(2.0 * np.pi * 1.5 * time)
        )[:, None]
        hd.create_dataset("moment0", data=np.ones_like(low))
        hd.create_dataset("moment2", data=moment_frequency ** 2)
    return holo, hd_path


def _spec(holo, *, artifacts=False):
    available, _ = load_pipeline_catalog()
    targets = [
        "velocity_analysis",
        "waveform_shape_metrics",
        "absolute_waveform_metrics",
        "blood_volume_rate",
        "lowrank_waveform_decomposition",
    ]
    velocity_options = ["segments"]
    if artifacts:
        targets.append("pdf_report")
        velocity_options += [
            "segment_velocity_maps",
            "velocity_profile_analysis",
            "velocity_profile_fft",
        ]
    return resolve_run_spec(
        input_paths=[holo],
        target_names=targets,
        pipelines=available,
        pipeline_options={
            "velocity_analysis": velocity_options,
            "waveform_shape_metrics": (),
            "absolute_waveform_metrics": (),
            "blood_volume_rate": ("masked_edges",),
            "lowrank_waveform_decomposition": (),
        },
        band_ratio_frequency_scale_hz=2.0,
    )


def test_gap_spectrum_is_computed_once_and_shared_by_both_velocity_workflows(tmp_path):
    holo, hd_path = _dual_input(tmp_path)
    with h5py.File(hd_path, "a") as hd:
        frames, height, width = hd["moment0"].shape
        spectra = hd.create_dataset("spectrograms", data=np.full(
            (frames, 3, height // 2, width // 2), 7, dtype=np.float32,
        ))
        spectra[frames // 2, 1] = 700
    dv_path = tmp_path / holo.stem / f"{holo.stem}_DV/h5/{holo.stem}_DV.h5"
    with h5py.File(dv_path, "a") as dv:
        shape = dv["segmentation/Retina/artery_mask"].shape
        dv.create_dataset("segmentation/Retina/vessel_segmentation_logits", data=np.zeros(shape))
        dv.create_dataset("segmentation/Choroid/vessel_mask", data=np.zeros(shape))
    with (
        patch(
            "pipelines.low_freq_spectra.runner.preprocess_spectrograms",
            wraps=preprocess_spectrograms,
        ) as preprocess,
        patch(
            "pipelines.low_freq_spectra.gap_analysis.gap_spectrograms_per_beat",
            wraps=gap_spectrograms_per_beat,
        ) as subdivide,
    ):
        result = execute_run(_spec(holo))
    assert result.succeeded, result.failures
    preprocess.assert_called_once()
    subdivide.assert_called_once()
    with h5py.File(result.outputs[0], "r") as output:
        schema = EyeFlowOutputPaths.active()
        spectrum = output[schema.spectrograms_masked]
        assert spectrum.shape == (frames, 3)
        expected = np.full((frames, 3), 7, dtype=np.float32)
        expected[frames // 2, 1] = 700
        np.testing.assert_array_equal(spectrum[:], expected)
        np.testing.assert_array_equal(output[schema.spectrograms_spatial_average][:], expected)
        np.testing.assert_array_equal(output[schema.spectrograms_spatial_median][:], expected)
        assert spectrum.name == "/ProcessingSpectrograms/Gaps/SpatialMedian"
        assert "ProcessingSpectrograms/Gaps/Raw" not in output
        assert "velocity_estimation_method" not in spectrum.attrs
        assert "ProcessingMoments/Velocity" in output
        assert "ProcessingBandRatio/Velocity" in output
        assert "spectrogramsMasked" not in output
        assert "spectrogramsMaskedPerBeat" not in output
        assert "ProcessingBandRatio/Spectrograms" not in output
        raw = output[schema.spectrograms_masked_per_beat]
        boundaries = output[schema.cardiac_cycle.systolic_peak_frame_indices][:]
        assert raw.shape == (128, 3, len(boundaries) - 1)
        median = output[schema.spectrograms_masked_time_aggregated_median]
        assert median.name == "/ProcessingSpectrograms/Gaps/TimeAggregated/MedianPhaseBeat"
        assert median.shape == (3,)
        np.testing.assert_array_equal(median[:], [7, 7, 7])
        assert median.attrs["source_path"] == raw.name
        assert "ProcessingSpectrograms/Gaps/TimeAggregated/Median" not in output
        assert "ProcessingBandRatio/Spectrograms/Gaps/TimeAggregated/MedianPhaseBeat" not in output
        beat_median = output[schema.spectrograms_masked_time_aggregated_median_beat]
        assert beat_median.name == "/ProcessingSpectrograms/Gaps/TimeAggregated/MedianBeat"
        assert beat_median.shape == (128, 3)
        assert beat_median.attrs["source_path"] == raw.name
        assert "ProcessingBandRatio/Spectrograms/Gaps/TimeAggregated/MedianBeat" not in output
        for path in (
            schema.spectrograms_masked_time_aggregated_median_early,
            schema.spectrograms_masked_time_aggregated_median_late,
        ):
            window = output[path]
            assert window.shape == (3,)
            np.testing.assert_array_equal(window[:], [7, 7, 7])
            assert window.attrs["source_path"] == beat_median.name
        spectrum_metrics = output[schema.spectrograms_masked_metrics_root]
        assert len(spectrum_metrics) == 33
        assert spectrum_metrics["index_max1"][()] == 0
        assert spectrum_metrics["index_max2"][()] == 1
        assert spectrum_metrics["R1"][()] == 1
        assert spectrum_metrics["R2"][()] == 1
        assert "ProcessingBandRatio/Spectrograms/Gaps/Metrics" not in output
        for path, shape in (
            (schema.spectrograms_masked_per_beat_band_limited, raw.shape),
            (schema.spectrograms_masked_per_beat_svd_continuous, (128 * raw.shape[2], raw.shape[1])),
        ):
            assert output[path].shape == shape
            assert output[path].attrs["source_path"] == "/" + schema.spectrograms_masked_per_beat
            assert output[path].attrs["cycle_boundary_workflow"] == "doppler_moments"
        np.testing.assert_array_equal(raw.attrs["cycle_boundary_indices"], boundaries)
        band_boundary_path = processing_path(
            schema.cardiac_cycle.systolic_peak_frame_indices, "ProcessingBandRatio",
        )
        assert not np.array_equal(
            boundaries, output[band_boundary_path][:],
        )
        assert np.max(raw[:]) > 7
        assert "velocity_quantity" not in raw.attrs
        assert raw.name == "/ProcessingSpectrograms/Gaps/PerBeat/Raw"
        assert output[schema.spectrograms_masked_per_beat_band_limited].attrs["band_limit_harmonics"] == 12
        assert output[schema.spectrograms_masked_per_beat_band_limited].attrs["band_limit_reconstruction"] == "signed"
        for root in ("ProcessingMoments", "ProcessingBandRatio"):
            assert f"{root}/VelocityPerBeat" not in output
            assert f"{root}/VelocityPerBeatSafe" not in output
            assert f"{root}/Velocity/Continuous/Global/BandLimited" not in output
            for paths in (schema.artery_per_beat, schema.artery_per_beat_filtered):
                signal = output[processing_path(paths.velocity_signal, root)]
                limited = output[processing_path(paths.velocity_signal_band_limited, root)]
                assert signal.shape == limited.shape
                assert limited.attrs["band_limit_harmonics"] == 12
                assert limited.attrs["band_limit_reconstruction"] == "matlab_magnitude"
                assert limited.attrs["source_path"].lstrip("/") == signal.name.lstrip("/")
        svd = output[schema.spectrograms_masked_per_beat_svd_continuous]
        assert svd.name == "/ProcessingSpectrograms/Gaps/PerBeat/SVDContinuous"
        assert "ProcessingSpectrograms/Gaps/PerBeat/SVD" not in output
        assert list(svd.attrs["dimDesc"]) == ["time", "freq"]
        assert svd.attrs["source_beat_count"] == raw.shape[2]
        assert svd.attrs["phase_samples_per_beat"] == 128
        assert svd.attrs["svd_modes_removed"] == 2


@pytest.mark.parametrize("empty_gaps", [False, True])
def test_gap_uses_bandratio_without_moments_and_skips_empty_support(tmp_path, empty_gaps):
    holo, hd_path = _dual_input(tmp_path)
    with h5py.File(hd_path, "a") as hd:
        if not empty_gaps:
            del hd["moment2"]
        frames, height, width = hd["moment0"].shape
        hd.create_dataset("spectrograms", data=np.ones(
            (frames, 2, height // 2, width // 2), dtype=np.float32,
        ))
    dv_path = tmp_path / holo.stem / f"{holo.stem}_DV/h5/{holo.stem}_DV.h5"
    with h5py.File(dv_path, "a") as dv:
        shape = dv["segmentation/Retina/artery_mask"].shape
        dv.create_dataset(
            "segmentation/Retina/vessel_segmentation_logits",
            data=np.ones(shape) if empty_gaps else np.zeros(shape),
        )
        dv.create_dataset("segmentation/Choroid/vessel_mask", data=np.zeros(shape))
    with (
        patch("pipelines.low_freq_spectra.gap_analysis.gap_spectrograms_per_beat",
              wraps=gap_spectrograms_per_beat) as subdivide,
        patch(
            "pipelines.low_freq_spectra.runner.preprocess_spectrograms",
            wraps=preprocess_spectrograms,
        ) as preprocess,
    ):
        result = execute_run(_spec(holo))
    assert result.succeeded, result.failures
    if empty_gaps:
        subdivide.assert_not_called()
    else:
        subdivide.assert_called_once()
    preprocess.assert_called_once()
    schema = EyeFlowOutputPaths.active()
    with h5py.File(result.outputs[0], "r") as output:
        if empty_gaps:
            assert schema.spectrograms_masked not in output
            assert "ProcessingMoments/Velocity" in output
        else:
            assert schema.spectrograms_masked in output
        assert schema.spectrograms_spatial_average in output
        assert schema.spectrograms_spatial_median in output
        if empty_gaps:
            assert schema.spectrograms_masked_per_beat not in output
        else:
            raw = output[schema.spectrograms_masked_per_beat]
            boundary_path = processing_path(
                schema.cardiac_cycle.systolic_peak_frame_indices, "ProcessingBandRatio",
            )
            assert raw.attrs["cycle_boundary_workflow"] == "frequency_bands"
            assert raw.attrs["cycle_boundary_source_path"] == "/" + boundary_path
            np.testing.assert_array_equal(raw.attrs["cycle_boundary_indices"], output[boundary_path][:])
        assert "ProcessingBandRatio/Velocity" in output


def test_malformed_spectra_do_not_disable_velocity_or_its_downstream_analysis(tmp_path):
    holo, hd_path = _dual_input(tmp_path)
    with h5py.File(hd_path, "a") as hd:
        hd.create_dataset("spectrograms", data=np.ones((4, 4), dtype=np.float32))
    available, _ = load_pipeline_catalog()
    spec = resolve_run_spec(
        input_paths=[holo], target_names=["velocity_analysis"], pipelines=available,
        pipeline_options={"velocity_analysis": ()}, band_ratio_frequency_scale_hz=2,
    )
    result = execute_run(spec)
    assert not result.succeeded
    with h5py.File(result.outputs[0], "r") as output:
        statuses = json.loads(output.attrs["pipeline_status"])
        assert statuses["shared"]["velocity"] == "succeeded"
        assert statuses["shared"]["low_freq_spectra"] == "failed"
        for method, root in (("doppler_moments", "ProcessingMoments"),
                             ("frequency_bands", "ProcessingBandRatio")):
            assert statuses["variants"][method]["velocity_analysis"] == "succeeded"
            assert f"{root}/Velocity" in output
        assert "ProcessingSpectrograms" not in output


def test_shared_spectra_survive_moments_workflow_failure_under_processing(tmp_path):
    holo, hd_path = _dual_input(tmp_path)
    with h5py.File(hd_path, "a") as hd:
        frames = hd["moment0"].shape[0]
        hd.create_dataset("spectrograms", data=np.full((frames, 3, 2, 2), 7, np.float32))

    def pack(result):
        if result.provenance["velocity_estimation_method"] == "doppler_moments":
            raise RuntimeError("synthetic moments output failure")
        return pack_velocity_outputs(result)

    available, _ = load_pipeline_catalog()
    spec = resolve_run_spec(
        input_paths=[holo], target_names=["velocity_analysis"], pipelines=available,
        pipeline_options={"velocity_analysis": ()},
        band_ratio_frequency_scale_hz=2,
    )
    with patch("pipelines.velocity.runner.pack_velocity_outputs", side_effect=pack):
        result = execute_run(spec)
    assert result.succeeded, result.failures
    schema = EyeFlowOutputPaths.active()
    with h5py.File(result.outputs[0], "r") as output:
        assert "ProcessingMoments/Velocity" not in output
        assert "ProcessingBandRatio/Velocity" in output
        np.testing.assert_array_equal(output[schema.spectrograms_spatial_average][:], 7)
        np.testing.assert_array_equal(output[schema.spectrograms_spatial_median][:], 7)
        assert "ProcessingBandRatio/Spectrograms" not in output


def test_both_workflows_use_independent_physical_velocity_cycles_and_downstream_beats(tmp_path):
    holo, _ = _dual_input(tmp_path)
    events = []
    detected = {}

    def cycles(data, *, dt_seconds):
        method = data.provenance["velocity_estimation_method"]
        events.append(("cycles", method))
        result = detect_velocity_cardiac_cycles(data, dt_seconds=dt_seconds)
        detected[method] = result[0]
        return result

    def estimate(**kwargs):
        events.append(("estimate", kwargs["velocity_estimation_method"]))
        return estimate_retinal_velocity(**kwargs)

    with (
        patch("pipelines.velocity.runner.detect_velocity_cardiac_cycles", side_effect=cycles),
        patch("pipelines.velocity.runner.estimate_retinal_velocity", side_effect=estimate),
    ):
        result = execute_run(_spec(holo, artifacts=True))
    assert result.succeeded, result.failures
    assert events == [
        ("estimate", "doppler_moments"),
        ("estimate", "frequency_bands"),
        ("cycles", "doppler_moments"),
        ("cycles", "frequency_bands"),
    ]
    schema = EyeFlowOutputPaths.active()
    with h5py.File(result.outputs[0], "r") as output:
        assert "Segmentation" in output
        assert "velocity_estimation_method" not in output.attrs
        assert json.loads(output.attrs["execution_variant_failures"]) == {}
        assert list(output.attrs["execution_variants_completed"]) == [
            "doppler_moments",
            "frequency_bands",
        ]
        for method, root in (
            ("doppler_moments", "ProcessingMoments"),
            ("frequency_bands", "ProcessingBandRatio"),
        ):
            assert output[root].attrs["velocity_estimation_method"] == method
            assert output[root].attrs["cardiac_cycle_detection_method"] == method
            assert output[root].attrs["cardiac_cycle_detection_quantity"] == "physical_velocity"
            assert output[root].attrs["cardiac_cycle_detection_unit"] == "mm/s"
            assert "Segmentation" not in output[root]
            for family in (
                schema.absolute_waveform_metrics_root,
                schema.waveform_shape_metrics_root,
                schema.lowrank_waveform_decomposition_root,
                schema.artery_velocity_profiles.transverse_velocity_profile_masked,
            ):
                assert processing_path(family, root) in output

            cycle = detected[method]
            expected_timing = {
                schema.cardiac_cycle.systolic_peak_frame_indices: cycle.systole.systole_indexes,
                schema.cardiac_cycle.systolic_cycle_duration_seconds: (
                    np.diff(cycle.systole.systole_indexes) * 0.05
                ),
                schema.cardiac_cycle.spectral_fundamental_frequency_hz: (
                    cycle.spectral.fundamental_hz
                ),
                schema.cardiac_cycle.spectral_heart_rate_bpm: cycle.spectral.heart_rate_bpm,
                schema.cardiac_cycle.spectral_heart_rate_standard_error_bpm: (
                    cycle.spectral.heart_rate_ste_bpm
                ),
                schema.cardiac_cycle.spectral_period_seconds: cycle.spectral.period_seconds,
            }
            for path, expected in expected_timing.items():
                np.testing.assert_allclose(output[processing_path(path, root)][()], expected)
            beat_count = cycle.systole.systole_indexes.size - 1
            waveform = output[processing_path(schema.artery_per_beat.velocity_signal, root)]
            flow = output[processing_path(schema.blood_volume_rate.artery.masked_edges, root)]
            assert waveform.shape[0] == beat_count
            assert flow.shape[1] == beat_count

            def check_dataset(_name, value, expected_method=method, expected_beats=beat_count):
                if isinstance(value, h5py.Dataset):
                    assert value.attrs["velocity_estimation_method"] == expected_method
                    assert value.attrs["cardiac_cycle_detection_quantity"] == "physical_velocity"
                    assert value.attrs["cardiac_cycle_detection_unit"] == "mm/s"
                    dimensions = list(value.attrs.get("dimDesc", ()))
                    if "beat" in dimensions:
                        assert value.shape[dimensions.index("beat")] == expected_beats

            output[root].visititems(check_dataset)
        for path in (
            schema.cardiac_cycle.systolic_peak_frame_indices,
            schema.cardiac_cycle.systolic_cycle_duration_seconds,
        ):
            assert not np.array_equal(
                output[path][:], output[processing_path(path, "ProcessingBandRatio")][:]
            )
        for path in (
            schema.analysis.retinal_artery_velocity_signal,
            schema.artery_per_beat.velocity_signal,
            schema.blood_volume_rate.artery.masked_edges,
        ):
            moments = output[path][:]
            bands = output[processing_path(path, "ProcessingBandRatio")][:]
            assert moments.shape != bands.shape or not np.allclose(moments, bands, equal_nan=True)
    ef_dir = result.outputs[0].parent.parent
    assert all(path.name.startswith(f"{holo.stem}_") for path in ef_dir.rglob("*") if path.is_file())
    for folder in ("moments", "bandratio"):
        assert list((ef_dir / "png" / folder).rglob("*.png"))
        assert list((ef_dir / "avi" / folder).rglob("*.avi"))
        assert list((ef_dir / "pdf" / folder).glob("*.pdf"))
    moments = _extract_parameters_from_h5(result.outputs[0], processing_root="ProcessingMoments")
    bands = _extract_parameters_from_h5(result.outputs[0], processing_root="ProcessingBandRatio")
    assert (
        moments["Average_Arterial_Velocity"]["value"] != bands["Average_Arterial_Velocity"]["value"]
    )


@pytest.mark.parametrize("failed_method", ["doppler_moments", "frequency_bands"])
@pytest.mark.parametrize("stage", ["input", "cycles", "estimate", "downstream"])
def test_one_failed_workflow_preserves_other_outputs(tmp_path, failed_method, stage):
    holo, hd_path = _dual_input(tmp_path)
    if stage == "input":
        with h5py.File(hd_path, "a") as hd:
            del hd["moment2" if failed_method == "doppler_moments" else "band_1_9000_18000"]

    def estimate(**kwargs):
        if stage == "estimate" and kwargs["velocity_estimation_method"] == failed_method:
            raise ValueError("injected estimator failure")
        return estimate_retinal_velocity(**kwargs)

    def cycles(data, *, dt_seconds):
        if stage == "cycles" and data.provenance["velocity_estimation_method"] == failed_method:
            raise ValueError("injected cycle detection failure")
        return detect_velocity_cardiac_cycles(data, dt_seconds=dt_seconds)

    from pipelines.absolute_waveform_metrics.runner import run_absolute_waveform_metrics

    def downstream(ctx):
        if stage == "downstream" and ctx.velocity_estimation_method == failed_method:
            # Exercise cleanup of a partial write and an already exported artifact.
            ctx.output.h5.write("ProcessingMoments/Partial/value", np.arange(3))
            raise ValueError("injected downstream failure")
        return run_absolute_waveform_metrics(ctx)

    spec = _spec(holo)
    with (
        patch("pipelines.velocity.runner.detect_velocity_cardiac_cycles", side_effect=cycles),
        patch("pipelines.velocity.runner.estimate_retinal_velocity", side_effect=estimate),
        patch(
            "pipelines.absolute_waveform_metrics.run_absolute_waveform_metrics",
            side_effect=downstream,
        ),
    ):
        result = execute_run(spec)
    assert result.succeeded, result.failures
    failed_root = "ProcessingMoments" if failed_method == "doppler_moments" else "ProcessingBandRatio"
    successful_root = "ProcessingBandRatio" if failed_method == "doppler_moments" else "ProcessingMoments"
    with h5py.File(result.outputs[0], "r") as output:
        if stage == "downstream":
            schema = EyeFlowOutputPaths.active()
            assert processing_path(schema.analysis.retinal_artery_velocity_signal, failed_root) in output
            assert processing_path(schema.blood_volume_rate.artery.masked_edges, failed_root) in output
            assert processing_path(schema.absolute_waveform_metrics_root, failed_root) not in output
            assert f"{failed_root}/Partial" not in output
            statuses = json.loads(output.attrs["pipeline_status"])["variants"][failed_method]
            assert statuses["absolute_waveform_metrics"] == "failed"
            assert statuses["blood_volume_rate"] == "succeeded"
        else:
            assert failed_root not in output
        assert successful_root in output
        assert output[successful_root].attrs["cardiac_cycle_detection_method"] != failed_method
        assert failed_method in json.loads(output.attrs["execution_variant_failures"])
        assert (
            processing_path(
                EyeFlowOutputPaths.active().blood_volume_rate.artery.masked_edges, successful_root
            )
            in output
        )
        if stage == "input" and failed_method == "doppler_moments":
            assert (
                output[successful_root].attrs["cardiac_cycle_detection_method"] == "frequency_bands"
            )
    ef_dir = result.outputs[0].parent.parent
    failed_folder = "moments" if failed_method == "doppler_moments" else "bandratio"
    for kind in ("png", "avi", "pdf", "eps"):
        if stage != "downstream":
            assert not (ef_dir / kind / failed_folder).exists()
    if stage == "downstream":
        assert list((ef_dir / "png" / failed_folder).glob("*.png"))


def test_both_workflows_fail_the_run(tmp_path):
    holo, _ = _dual_input(tmp_path)
    with patch(
        "pipelines.velocity.runner.estimate_retinal_velocity",
        side_effect=ValueError("injected failure"),
    ):
        result = execute_run(_spec(holo))
    assert not result.succeeded
    assert "Neither velocity workflow completed" in result.failures[0].message


def test_pulse_export_failure_keeps_persisted_analysis_and_downstream_outputs(tmp_path):
    holo, _ = _dual_input(tmp_path)
    schema = EyeFlowOutputPaths.active()
    observed = []

    def fail_after_persistence(ctx, analysis, per_beat):
        assert ctx.output.h5.get(schema.artery_per_beat.velocity_signal) is not None
        assert ctx.output.h5.get(schema.analysis.retinal_artery_velocity_signal) is not None
        observed.append(ctx.velocity_estimation_method)
        raise ValueError("injected PNG failure")

    with patch("pipelines.velocity_analysis.runner._export_pulse_pngs", side_effect=fail_after_persistence):
        result = execute_run(_spec(holo))
    assert result.succeeded, result.failures
    assert observed == ["doppler_moments", "frequency_bands"]
    with h5py.File(result.outputs[0], "r") as output:
        assert not json.loads(output.attrs["pipeline_failures"])
        for root in ("ProcessingMoments", "ProcessingBandRatio"):
            assert processing_path(schema.artery_per_beat.velocity_signal, root) in output
            assert processing_path(schema.blood_volume_rate.artery.masked_edges, root) in output


@pytest.mark.parametrize("corruption", ["moment_rank", "band_values", "band_shape"])
def test_malformed_input_skips_only_its_workflow(tmp_path, corruption):
    holo, hd_path = _dual_input(tmp_path)
    with h5py.File(hd_path, "a") as hd:
        if corruption == "moment_rank":
            del hd["moment0"]
            hd.create_dataset("moment0", data=np.ones((64, 64), dtype=np.float32))
        elif corruption == "band_values":
            hd["band_0_3000_9000"][0, 0, 0] = np.nan
        else:
            del hd["band_1_9000_18000"]
            hd.create_dataset("band_1_9000_18000", data=np.ones((10, 64, 64)))
    result = execute_run(_spec(holo))
    assert result.succeeded, result.failures
    failed_root = "ProcessingMoments" if corruption == "moment_rank" else "ProcessingBandRatio"
    successful_root = "ProcessingBandRatio" if corruption == "moment_rank" else "ProcessingMoments"
    with h5py.File(result.outputs[0], "r") as output:
        assert failed_root not in output
        assert successful_root in output
