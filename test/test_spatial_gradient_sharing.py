"""Continuous gradients are shared across methods; beat products stay local."""

import json
from dataclasses import replace
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from calculations.topology import OpticDisc
from calculations.vessel_segments.measurement import SegmentMeasurementSettings
from calculations.vessel_segments.sampling.transforms import (
    interpolate_segment_masks, rotate_segment_masks,
)
from input_output.output_manager import OutputManager
from input_output.runs.layout import HoloRunLayout
from input_output.schema import PixelPitch
from pipeline_engine import ExecutionVariant, PipelineContext
from pipeline_engine.runtime import run_pipelines_to_output
from pipelines.spatial_gradient_moment0 import profiles, runner
from pipelines.spatial_gradient_moment0.cache import continuous_gradient_key
from test_pipeline_failure_isolation import _descriptor
from test_topology_chunks import _prepared


@pytest.fixture
def gradient_input(tmp_path):
    plan = _prepared()
    masks = interpolate_segment_masks(plan.native.segment_masks, 16)
    plan = replace(
        plan, interpolated_masks=masks,
        rotated_masks=rotate_segment_masks(masks, plan.rotation_degrees),
    )
    vessel = plan.native.labels > 0
    source = SimpleNamespace(
        segmentation=SimpleNamespace(
            vessels=SimpleNamespace(artery=vessel, vein=vessel.copy()),
            optic_disc=OpticDisc(None, (4.0, 8.0), 2.0, 2.0),
        ),
        holodoppler=SimpleNamespace(
            pixel_pitch=PixelPitch(20e-6, 20e-6),
            timing=SimpleNamespace(dt_seconds=0.1),
        ),
    )
    path = tmp_path / "hd.h5"
    values = np.arange(9 * 9, dtype=np.float32).reshape(9, 9) ** 2
    with h5py.File(path, "w") as hd:
        hd.create_dataset("moment0ff", data=np.stack([
            values * (1 + frame / 8) for frame in range(8)
        ]))
    return path, source, {"artery": plan, "vein": plan}


def _extract(ctx, source, plans, settings=SegmentMeasurementSettings(0.02)):
    return profiles.extract_spatial_gradient_segments(
        ctx, source, plans, profile_settings=settings,
    )


def test_reuses_identical_continuous_arrays_and_preserves_method_specific_beats(gradient_input):
    path, source, plans = gradient_input
    with h5py.File(path) as hd, h5py.File(
        "gradient-sharing.h5", "w", driver="core", backing_store=False,
    ) as output:
        cache = {}
        contexts = [PipelineContext(
            work_h5=output, holodoppler_h5=hd, doppler_vision_h5=None,
            computation_cache=cache, velocity_estimation_method=method,
        ) for method in ("doppler_moments", "frequency_bands")]
        with patch.object(
            profiles, "analyze_segment_profiles", wraps=profiles.analyze_segment_profiles,
        ) as analyze:
            first = _extract(contexts[0], source, plans)
            second = _extract(contexts[1], source, plans)
        analyze.assert_called_once()
        assert first is second
        for segment in first:
            assert np.any(np.isfinite(segment.transverse.masked))
            for pair in (segment.transverse, segment.longitudinal, segment.mean_images):
                assert not pair.unmasked.flags.writeable
                assert not pair.masked.flags.writeable

        # Compare with independent extraction, not just another cache lookup.
        independent = PipelineContext(
            work_h5=output, holodoppler_h5=hd, doppler_vision_h5=None,
        )
        baseline = _extract(independent, source, plans)
        for shared, separate in zip(first, baseline):
            np.testing.assert_array_equal(shared.segment_signal, separate.segment_signal)
            for name in ("transverse", "longitudinal", "mean_images"):
                for masked in (False, True):
                    np.testing.assert_array_equal(
                        getattr(shared, name).select(masked=masked),
                        getattr(separate, name).select(masked=masked),
                    )
        for timing in ([0, 4, 7], [0, 7]):
            actual = profiles.pack_spatial_gradient_profile_outputs(*second, timing)
            expected = profiles.pack_spatial_gradient_profile_outputs(*baseline, timing)
            for name in expected:
                np.testing.assert_array_equal(actual[name].data, expected[name].data)
                assert actual[name].attrs == expected[name].attrs
        moments = profiles.pack_spatial_gradient_profile_outputs(*first, [0, 4, 7])
        bands = profiles.pack_spatial_gradient_profile_outputs(*second, [0, 7])
        profile_path = next(name for name in moments if "Masked/SpatialGradientProfile/value" in name)
        assert moments[profile_path].data.shape[2] == 2
        assert bands[profile_path].data.shape[2] == 1


@pytest.mark.parametrize("changed", ["dataset", "mask", "rotation", "settings", "filter"])
def test_cache_requires_compatible_source_masks_topology_and_settings(gradient_input, changed):
    path, source, plans = gradient_input
    masks = {name: getattr(source.segmentation.vessels, name) for name in plans}
    settings = SegmentMeasurementSettings(0.02)
    filters = (7, 6, 8, 0.6, 5)
    with h5py.File(path, "a") as hd:
        original = continuous_gradient_key(hd["moment0ff"], masks, plans, settings, filters)
        dataset = hd["moment0ff"]
        # Distinct Python wrappers for the same HDF5 object must still hit.
        assert original == continuous_gradient_key(hd["moment0ff"], masks, plans, settings, filters)
        if changed == "dataset":
            dataset = hd.create_dataset("other", data=dataset[:])
        elif changed == "mask":
            masks = {**masks, "artery": ~masks["artery"]}
        elif changed == "rotation":
            plans = {**plans, "artery": replace(
                plans["artery"], rotation_degrees=plans["artery"].rotation_degrees + 1,
            )}
        elif changed == "settings":
            settings = replace(settings, pixel_size_mm=0.03)
        else:
            filters = (5, 6, 8, 0.6, 5)
        assert original != continuous_gradient_key(dataset, masks, plans, settings, filters)


@pytest.mark.parametrize("failure", [None, "extraction", "packing", "persistence"])
def test_runtime_sharing_keeps_variant_failures_and_retries_independent(
    gradient_input, tmp_path, failure,
):
    path, source, plans = gradient_input
    products = {}
    calls = []
    caches = []
    real_analyze = profiles.analyze_segment_profiles
    real_pack = runner.pack_spatial_gradient_profile_outputs

    def produce(ctx):
        ctx.publish_execution_variants([
            ExecutionVariant(
                name=method, state={}, output_namespace=root, artifact_namespace=method,
                provenance={"velocity_estimation_method": method},
            ) for method, root in (
                ("doppler_moments", "ProcessingMoments"),
                ("frequency_bands", "ProcessingBandRatio"),
            )
        ])

    def analyze(*args, **kwargs):
        calls.append(True)
        if failure == "extraction" and len(calls) == 1:
            raise ValueError("first extraction failed")
        return real_analyze(*args, **kwargs)

    def pack(*args, **kwargs):
        if failure == "packing" and len(caches) == 1:
            raise ValueError("first packing failed")
        return real_pack(*args, **kwargs)

    def gradient(ctx):
        caches.append(ctx.computation_cache)
        outputs = runner.run_spatial_gradient_moment0(ctx)
        if failure == "persistence" and len(caches) == 1:
            return {**outputs, "Invalid/value": object()}
        products[ctx.execution_variant.name] = ctx.state[runner.SPATIAL_GRADIENT_PRODUCTS_STATE]
        return outputs

    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    with (
        patch.object(runner, "load_velocity_inputs", return_value=source),
        patch.object(runner, "prepared_topologies", return_value=plans),
        patch.object(runner, "cardiac_cycle_indexes", side_effect=lambda ctx: (
            np.array([0, 4, 7]) if ctx.execution_variant.name == "doppler_moments"
            else np.array([0, 7])
        )),
        patch.object(runner, "export_lumen_size_pngs"),
        patch.object(profiles, "analyze_segment_profiles", side_effect=analyze),
        patch.object(runner, "pack_spatial_gradient_profile_outputs", side_effect=pack),
    ):
        def run():
            return run_pipelines_to_output(
                output_manager=manager, holodoppler_h5=path, doppler_vision_h5=None,
                target_names=("spatial_gradient_moment0",), copy_hd_passthrough=False,
                pipelines=[
                    _descriptor("source", produce, variants=True),
                    _descriptor("spatial_gradient_moment0", gradient, ("source",)),
                ],
            )
        output_path = run()
        if failure is None:
            assert products["doppler_moments"].artery_segments is products["frequency_bands"].artery_segments
        else:
            assert set(products) == {"frequency_bands"}
        assert caches[0] is caches[1]
        assert len(calls) == (2 if failure == "extraction" else 1)
        with h5py.File(output_path) as output:
            root = "SpatialGradientProfiles/Artery/Transverse/Masked/SpatialGradientProfile/value"
            assert output[f"ProcessingBandRatio/{root}"].shape[2] == 1
            if failure is None:
                assert output[f"ProcessingMoments/{root}"].shape[2] == 2
            else:
                assert f"ProcessingMoments/{root}" not in output
                statuses = json.loads(output.attrs["pipeline_status"])["variants"]
                assert statuses["doppler_moments"]["spatial_gradient_moment0"] == "failed"
                assert statuses["frequency_bands"]["spatial_gradient_moment0"] == "succeeded"
        # A second input run must compute again, even with the same source path.
        if failure is None:
            run()
            assert len(calls) == 2
            assert caches[2] is caches[3] and caches[2] is not caches[0]
