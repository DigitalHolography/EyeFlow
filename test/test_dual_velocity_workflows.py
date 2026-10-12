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
from pipelines.velocity.runner import detect_velocity_cardiac_cycles, estimate_retinal_velocity


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
            ("doppler_moments", "Processing"),
            ("frequency_bands", "ProcessingAlt"),
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
                output[path][:], output[processing_path(path, "ProcessingAlt")][:]
            )
        for path in (
            schema.analysis.retinal_artery_velocity_signal,
            schema.artery_per_beat.velocity_signal,
            schema.blood_volume_rate.artery.masked_edges,
        ):
            moments = output[path][:]
            bands = output[processing_path(path, "ProcessingAlt")][:]
            assert moments.shape != bands.shape or not np.allclose(moments, bands, equal_nan=True)
    ef_dir = result.outputs[0].parent.parent
    assert all(path.name.startswith(f"{holo.stem}_") for path in ef_dir.rglob("*") if path.is_file())
    for folder in ("moments", "bandratio"):
        assert list((ef_dir / "png" / folder).rglob("*.png"))
        assert list((ef_dir / "avi" / folder).rglob("*.avi"))
        assert list((ef_dir / "pdf" / folder).glob("*.pdf"))
    moments = _extract_parameters_from_h5(result.outputs[0], processing_root="Processing")
    bands = _extract_parameters_from_h5(result.outputs[0], processing_root="ProcessingAlt")
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
            ctx.output.h5.write("Processing/Partial/value", np.arange(3))
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
    failed_root = "Processing" if failed_method == "doppler_moments" else "ProcessingAlt"
    successful_root = "ProcessingAlt" if failed_method == "doppler_moments" else "Processing"
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
        for root in ("Processing", "ProcessingAlt"):
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
    failed_root = "Processing" if corruption == "moment_rank" else "ProcessingAlt"
    successful_root = "ProcessingAlt" if corruption == "moment_rank" else "Processing"
    with h5py.File(result.outputs[0], "r") as output:
        assert failed_root not in output
        assert successful_root in output
