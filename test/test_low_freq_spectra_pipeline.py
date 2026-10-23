"""Shared spectral scheduling, timing fallback, and failure isolation."""

import json
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest
from test_run_service import _write_input
from test_velocity_pipeline import _retinal_velocity

from input_output.schema import EyeFlowOutputPaths
from input_output.schema.eyeflow_output import VELOCITY_WORKFLOW_ROOTS, processing_path
from pipeline_engine import ExecutionVariant, PipelineDAG, PipelineDescriptor
from pipeline_engine.run_service import execute_run, resolve_run_spec
from pipelines import load_pipeline_catalog
from pipelines.low_freq_spectra.runner import run_low_freq_spectra
from pipelines.low_freq_spectra.frequency import write_frequency_bin_table
from input_output.h5_access import PipelineH5Output
from pipelines.topology_core.runner import CAPILLARY_ANALYSIS_ZONES_STATE
from pipelines.velocity.models import VelocityWorkflows
from pipelines.velocity.runner import VELOCITY_WORKFLOWS_STATE


def _descriptor(name, run, **metadata):
    return PipelineDescriptor(
        name=name, description=name, available=True,
        pipeline_factory=lambda: SimpleNamespace(name=name, run=run, source_path=__file__),
        **metadata,
    )


def _spec(tmp_path, *, moments=True, bands=True, spectra="valid"):
    holo = _write_input(tmp_path)
    hd_path = tmp_path / "scan/scan_HD/h5/scan_HD_output.h5"
    with h5py.File(hd_path, "a") as hd:
        hd.create_dataset("HD_parameters", data=json.dumps({
            "sampling_freq": 37037.0, "batch_size": 32, "batch_stride": 256,
        }))
        if spectra == "legacy_group":
            hd.create_group("spectrograms")
        elif spectra != "missing":
            shape = {"rank": (6, 3), "empty": (0, 3, 2, 2), "frames": (7, 3, 2, 2)}.get(
                spectra, (6, 8, 2, 2),
            )
            hd.create_dataset("spectrograms", data=np.ones(
                shape, dtype=np.complex64 if spectra == "dtype" else np.float32,
            ))
    schema = EyeFlowOutputPaths.active()
    results = {}
    for method, available in (("doppler_moments", moments), ("frequency_bands", bands)):
        if not available:
            continue
        result = _retinal_velocity()
        if method == "frequency_bands":
            result.cardiac_cycle.systole.systole_indexes = np.array([0, 2, 5])
        if available == "undetected":
            result = replace(result, cardiac_cycle_source="none")
        elif available == "one_boundary":
            result.cardiac_cycle.systole.systole_indexes = np.array([0])
        results[method] = result

    def source(ctx):
        ctx.state.set(VELOCITY_WORKFLOWS_STATE, VelocityWorkflows(
            moments=results.get("doppler_moments"), bandratio=results.get("frequency_bands"),
        ))
        ctx.state.set(CAPILLARY_ANALYSIS_ZONES_STATE, np.ones((2, 2), dtype=bool))
        ctx.output.h5.write(schema.segmentation.gaps, np.ones((2, 2), dtype=bool))
        for method, result in results.items():
            ctx.output.h5.write(
                processing_path(schema.cardiac_cycle.systolic_peak_frame_indices,
                                VELOCITY_WORKFLOW_ROOTS[method]),
                result.cycle_boundary_indexes,
            )
        ctx.publish_execution_variants([
            ExecutionVariant(
                name=method, state={"velocity": result},
                output_namespace=VELOCITY_WORKFLOW_ROOTS[method],
                artifact_namespace=method, provenance={"velocity_estimation_method": method},
            ) for method, result in results.items()
        ])

    available, _ = load_pipeline_catalog()
    low_freq = next(item for item in available if item.name == "low_freq_spectra")
    return resolve_run_spec(
        input_paths=[holo], target_names=["consumer"], pipelines=[
            _descriptor("velocity", source, produces_execution_variants=True,
                        dag_produces=("velocity_workflows", "capillary_analysis_zones")),
            low_freq,
            _descriptor("consumer", lambda ctx: {"ProcessingMoments/Consumer/value": 42},
                        dag_requires=("velocity",)),
        ],
    )


@pytest.mark.parametrize("moments,bands,expected", [
    (True, True, "doppler_moments"),
    (True, False, "doppler_moments"),
    (False, True, "frequency_bands"),
    ("undetected", True, "frequency_bands"),
    ("one_boundary", True, "frequency_bands"),
    ("undetected", "undetected", None),
])
def test_shared_pipeline_prefers_moments_then_detected_bandratio(tmp_path, moments, bands, expected):
    spec = _spec(tmp_path, moments=moments, bands=bands)
    with patch("pipelines.low_freq_spectra.runner.run_low_freq_spectra",
               wraps=run_low_freq_spectra) as run:
        result = execute_run(spec)
    assert result.succeeded, result.failures
    run.assert_called_once()
    ctx = run.call_args.args[0]
    assert ctx.execution_variant is None
    assert ctx.output.h5.processing_root is None
    assert spec.plan.targets == ("consumer",)
    assert spec.plan.names.index("velocity") < spec.plan.names.index("low_freq_spectra")
    schema = EyeFlowOutputPaths.active()
    with h5py.File(result.outputs[0], "r") as output:
        assert output[schema.sampling_frequency_hz][()] == 37037
        assert output[schema.sampling_frequency_hz].attrs["unit"] == "Hz"
        status = json.loads(output.attrs["pipeline_status"])
        assert status["shared"]["low_freq_spectra"] == "succeeded"
        assert all("low_freq_spectra" not in branch for branch in status["variants"].values())
        assert schema.spectrograms_masked in output
        assert output[schema.spectrograms_masked].name == (
            "/ProcessingSpectrograms/Gaps/SpatialMedian"
        )
        assert "ProcessingSpectrograms/Gaps/Raw" not in output
        assert output[schema.spectrograms_masked].attrs["quantity"] == "spatial_median_band_power"
        assert output[schema.spectrograms_masked].attrs["spatial_aggregation_method"] == "median"
        assert "velocity_estimation_method" not in output[schema.spectrograms_masked].attrs
        if expected is None:
            assert schema.spectrograms_masked_per_beat not in output
            assert schema.spectrograms_masked_time_aggregated_median not in output
            assert schema.spectrograms_masked_time_aggregated_median_beat not in output
            assert schema.spectrograms_masked_time_aggregated_median_early not in output
            assert schema.spectrograms_masked_time_aggregated_median_late not in output
        else:
            root = VELOCITY_WORKFLOW_ROOTS[expected]
            boundary_path = processing_path(schema.cardiac_cycle.systolic_peak_frame_indices, root)
            for path in (schema.spectrograms_masked_per_beat,
                         schema.spectrograms_masked_per_beat_band_limited,
                         schema.spectrograms_masked_per_beat_svd_continuous,
                         schema.spectrograms_masked_time_aggregated_median,
                         schema.spectrograms_masked_time_aggregated_median_beat,
                         schema.spectrograms_masked_time_aggregated_median_early,
                         schema.spectrograms_masked_time_aggregated_median_late):
                product = output[path]
                assert product.attrs["quantity"] == "spatial_median_band_power"
                assert product.attrs["spatial_aggregation_method"] == "median"
                assert product.attrs["cycle_boundary_workflow"] == expected
                assert product.attrs["cycle_boundary_source_path"] == "/" + boundary_path
                np.testing.assert_array_equal(product.attrs["cycle_boundary_indices"],
                                              output[boundary_path][:])
        assert "ProcessingBandRatio/Spectrograms" not in output
        assert "ProcessingSpectrograms/Gaps/TimeAggregated/Median" not in output

    names = {"spectrograms_gaps", "spectrograms_gaps_monochrome", "gaps_mask"}
    if expected is not None:
        names.add("spectrogram_gaps_median")
        names.update(f"gaps_{name}" for name in ("A0", "A", "A2", "R1", "SE"))
        names.update(f"gaps_positive_freq_median_phase_{suffix}" for suffix in (
            "early_late", "beat", "all",
        ))
    ef_dir = result.outputs[0].parent.parent
    for extension in ("png", "eps"):
        expected_names = names | ({"gaps_qc"} if extension == "png" else set())
        assert {path.name for path in (ef_dir / extension / "spectrograms").iterdir()} == {
            f"scan_{name}.{extension}" for name in expected_names
        }
        assert not list((ef_dir / extension / "moments").glob("*spectrogram*"))
        assert not list((ef_dir / extension / "bandratio").glob("*spectrogram*"))


def test_positive_figures_export_for_32_bin_input_with_the_stored_peak(tmp_path):
    spec = _spec(tmp_path)
    with h5py.File(tmp_path / "scan/scan_HD/h5/scan_HD_output.h5", "a") as hd:
        del hd["spectrograms"]
        values = np.ones((6, 32, 2, 2))
        values[:, 5] = 10
        hd.create_dataset("spectrograms", data=values)
        hd["HD_parameters"][()] = json.dumps({"sampling_freq": 37037, "batch_size": 128})
    result = execute_run(spec)
    assert result.succeeded, result.failures
    schema = EyeFlowOutputPaths.active()
    with h5py.File(result.outputs[0], "r") as output:
        assert output[f"{schema.spectrograms_masked_metrics_root}/index_max1"][()] == 5
    ef_dir = result.outputs[0].parent.parent
    for extension in ("png", "eps"):
        for suffix in ("early_late", "beat", "all"):
            assert (ef_dir / extension / "spectrograms" / f"scan_gaps_positive_freq_median_phase_{suffix}.{extension}").is_file()


@pytest.mark.parametrize("bin_count", [2, 32, 64])
def test_frequency_bin_table_is_written_without_detected_timing_or_figure_export(tmp_path, bin_count):
    spec = _spec(tmp_path, moments="undetected", bands="undetected")
    with h5py.File(tmp_path / "scan/scan_HD/h5/scan_HD_output.h5", "a") as hd:
        del hd["spectrograms"]
        hd.create_dataset("spectrograms", data=np.ones((6, bin_count, 2, 2)))
        hd["HD_parameters"][()] = json.dumps({
            "sampling_freq": 64000, "batch_size": 4 * bin_count, "batch_stride": 256,
        })
    with patch("pipelines.low_freq_spectra.runner.export_spectrogram_figures", return_value=[]) as export:
        result = execute_run(spec)
    assert export.call_args.kwargs["dt_seconds"] == 256 / 64000
    assert result.succeeded, result.failures
    schema = EyeFlowOutputPaths.active()
    with h5py.File(result.outputs[0], "r") as output:
        table = output[schema.spectrograms_frequency_bin_table]
        assert table.name == "/ProcessingSpectrograms/FrequencyBinTable"
        assert table.shape == (bin_count,)
        assert table.dtype == np.float32
        raw = np.fft.fftfreq(4 * bin_count, d=1 / 64000)
        np.testing.assert_allclose(table[:], raw.reshape(bin_count, 4).mean(axis=1) / 1000)
        assert "column_names" not in table.attrs
        assert table.attrs["unit"] == "kHz"
        assert list(table.attrs["dimDesc"]) == ["freq"]
        assert table.attrs["index_base"] == 0
        assert table.attrs["sampling_frequency_hz"] == 64000
        assert table.attrs["sampling_frequency_source_path"] == "/" + schema.sampling_frequency_hz
        assert table.attrs["source_path"] == "/spectrograms"
        assert table.attrs["fft_bin_count"] == 4 * bin_count
        assert table.attrs["frequency_binning_factor"] == 4
        assert schema.spectrograms_spatial_average in output
        assert schema.spectrograms_masked_per_beat not in output


@pytest.mark.parametrize("sampling_frequency", [None, 0, np.nan])
def test_frequency_bin_table_requires_valid_calibration(tmp_path, sampling_frequency):
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as h5:
        if sampling_frequency is not None:
            h5.create_dataset(schema.sampling_frequency_hz, data=sampling_frequency)
        with patch("pipelines.low_freq_spectra.frequency.Logger.log_warning") as warning:
            write_frequency_bin_table(SimpleNamespace(h5=PipelineH5Output(h5)), 64)
        assert schema.spectrograms_frequency_bin_table not in h5
        warning.assert_called_once()


def test_masked_spatial_median_feeds_persisted_downstream_spectra(tmp_path):
    spec = _spec(tmp_path)
    with h5py.File(tmp_path / "scan/scan_HD/h5/scan_HD_output.h5", "a") as hd:
        # The median is 4 despite the included 1000-power outlier; mean is 252.25.
        hd["spectrograms"][:] = (
            np.array([[1, 3], [5, 1000]])[None, None]
            + 10 * np.arange(8)[None, :, None, None]
        )
    result = execute_run(spec)
    assert result.succeeded, result.failures
    schema = EyeFlowOutputPaths.active()
    expected = 4 + 10 * np.arange(8)
    with h5py.File(result.outputs[0], "r") as output:
        assert "ProcessingSpectrograms/Gaps/SpatialMedian" in output
        assert "ProcessingSpectrograms/Gaps/Raw" not in output
        np.testing.assert_array_equal(
            output[schema.spectrograms_masked][:], np.broadcast_to(expected, (6, 8)),
        )
        np.testing.assert_array_equal(
            output[schema.spectrograms_spatial_average][:],
            np.broadcast_to(expected + 248.25, (6, 8)),
        )
        np.testing.assert_array_equal(
            output[schema.spectrograms_masked_time_aggregated_median][:], expected,
        )
        np.testing.assert_array_equal(
            output[schema.spectrograms_masked_time_aggregated_median_beat][:],
            np.broadcast_to(expected, (128, 8)),
        )
        for path in (schema.spectrograms_masked_per_beat,
                     schema.spectrograms_masked_per_beat_band_limited):
            data = output[path][:]
            np.testing.assert_allclose(data, np.broadcast_to(expected[None, :, None], data.shape))
        assert output[schema.spectrograms_masked_per_beat].attrs["source_path"] == (
            "/ProcessingSpectrograms/Gaps/SpatialMedian"
        )
        data = output[schema.spectrograms_masked_per_beat_svd_continuous][:]
        np.testing.assert_allclose(data, np.broadcast_to(expected, data.shape))


@pytest.mark.parametrize("spectra", ["missing", "legacy_group"])
def test_optional_source_skips_all_spectral_processing(tmp_path, spectra):
    spec = _spec(tmp_path, spectra=spectra)
    with patch("pipelines.low_freq_spectra.runner.preprocess_spectrograms") as preprocess:
        result = execute_run(spec)
    assert result.succeeded, result.failures
    preprocess.assert_not_called()
    ef_dir = result.outputs[0].parent.parent
    assert not (ef_dir / "png/spectrograms").exists()
    assert not (ef_dir / "eps/spectrograms").exists()
    with h5py.File(result.outputs[0], "r") as output:
        assert "ProcessingSpectrograms" not in output
        assert output["ProcessingMoments/Consumer/value"][()] == 42
        assert output["ProcessingBandRatio/Consumer/value"][()] == 42


@pytest.mark.parametrize("spectra", ["rank", "dtype", "empty", "frames"])
def test_bad_spectra_fail_only_the_shared_pipeline(tmp_path, spectra):
    result = execute_run(_spec(tmp_path, spectra=spectra))
    assert not result.succeeded
    with h5py.File(result.outputs[0], "r") as output:
        status = json.loads(output.attrs["pipeline_status"])
        assert status["shared"]["low_freq_spectra"] == "failed"
        assert status["shared"]["velocity"] == "succeeded"
        assert all(branch["consumer"] == "succeeded" for branch in status["variants"].values())
        assert "ProcessingSpectrograms" not in output  # Rolled back, including raw reductions.


def test_automatic_companions_close_dependencies_and_do_not_become_targets():
    descriptors = [
        _descriptor("second", lambda ctx: None, auto_run_after="first"),
        _descriptor("first", lambda ctx: None, auto_run_after="source", dag_requires=("extra",)),
        _descriptor("extra", lambda ctx: None),
        _descriptor("source", lambda ctx: None),
        _descriptor("unrelated", lambda ctx: None),
    ]
    dag = PipelineDAG(descriptors)
    plan = dag.resolve_targets(["source"])
    assert plan.targets == ("source",)
    assert set(plan.names) == {"source", "extra", "first", "second"}
    assert plan.names.index("source") < plan.names.index("first") < plan.names.index("second")
    assert dag.resolve_targets(["unrelated"]).names == ("unrelated",)
    assert dag.resolve_targets([]).names == ()


def test_low_frequency_pipeline_cannot_be_selected_directly():
    available, _ = load_pipeline_catalog()
    with pytest.raises(ValueError, match="hidden.*low_freq_spectra"):
        resolve_run_spec(input_paths=[], target_names=["low_freq_spectra"], pipelines=available)


@pytest.mark.parametrize("class_pipeline", [False, True])
def test_decorators_preserve_companion_and_shared_execution_metadata(monkeypatch, class_pipeline):
    from pipeline_engine import base

    monkeypatch.setattr(base, "PIPELINE_REGISTRY", {})
    if class_pipeline:
        @base.registerPipeline(name="companion", auto_run_after="source", run_once=True)
        class Companion(base.ProcessPipeline):
            def run(self, ctx):
                return None
    else:
        @base.pipeline(name="companion", auto_run_after="source", run_once=True)
        def companion(ctx):
            return None
    descriptor = base.PIPELINE_REGISTRY["companion"]
    instance = descriptor.instantiate()
    assert descriptor.auto_run_after == instance.auto_run_after == "source"
    assert descriptor.run_once and instance.run_once
