"""DV alignment, topology ownership, and persisted capillary analysis zones."""

import json
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from input_output.h5_access import PipelineH5Output, PipelineInputSource, RawH5SourceReader
from input_output.schema import EyeFlowOutputPaths
from input_output.schema.doppler_view import (
    DV_CHOROID_VESSEL_MASK_PATH,
    DV_RETINAL_VESSEL_LOGITS_PATH,
)
from pipeline_engine import PipelineDAG
from pipeline_engine.context import PipelineInputs, PipelineState
from pipelines import load_pipeline_catalog
from pipelines.topology_core.runner import capillary_zones, run_topology_core
from pipelines.low_freq_spectra.spectrograms import preprocess_spectrograms


@pytest.fixture
def source_pair(tmp_path, request):
    shape = getattr(request, "param", (41, 61))
    with (
        h5py.File(tmp_path / "hd.h5", "w") as hd,
        h5py.File(tmp_path / "dv.h5", "w") as dv,
    ):
        hd.create_dataset("moment0", data=np.ones((2, *shape), dtype=np.float32))
        hd.create_dataset("moment2", data=np.ones((2, *shape), dtype=np.float32))
        hd.create_dataset(
            "HD_parameters",
            data=json.dumps({
                "pixel_pitch": [0.0002, 0.0002],
                "sampling_freq": 20,
                "batch_stride": 1,
            }),
        )
        dv.create_dataset("segmentation/Retina/artery_mask", data=np.zeros(shape, bool))
        dv.create_dataset("segmentation/Retina/vein_mask", data=np.zeros(shape, bool))
        yield hd, dv


def _ctx(hd, dv):
    return SimpleNamespace(
        inputs=PipelineInputs(
            hd=PipelineInputSource(RawH5SourceReader(h5file=hd, label="HD"), {}),
            dv=PipelineInputSource(RawH5SourceReader(h5file=dv, label="DV"), {}),
        ),
        state=PipelineState(),
        attrs={},
        require_inputs=lambda *slots: None,
    )


def _write_capillary_sources(dv, shape):
    logits = np.full(shape, -2.0, dtype=np.float32)
    logits[20, 12] = 0.25  # Above threshold, separated from the other samples.
    logits[20, 16] = 0.0
    logits[20, 20] = 0.075
    logits[20, 28] = 0.1
    logits[20, 36] = np.nextafter(np.float32(0.1), np.float32(np.inf))
    choroid = np.zeros(shape, dtype=np.int32)
    choroid[30, 30] = 7
    disc = np.zeros(shape, dtype=np.uint8)
    disc[12, 34] = 255
    dv.create_dataset(DV_RETINAL_VESSEL_LOGITS_PATH, data=logits)
    dv.create_dataset(DV_CHOROID_VESSEL_MASK_PATH, data=choroid)
    dv.create_dataset("segmentation/OpticDisc/mask", data=disc)
    return logits, choroid, disc


def _reference_zones(retinal, choroid, optic_disc):
    # Build the expected exclusion by stamping each vessel's 7x7 neighborhood,
    # independently of the production morphology implementation.
    ny, nx = retinal.shape
    excluded = optic_disc.copy()
    for y, x in zip(*np.nonzero(retinal | choroid)):
        excluded[max(0, y - 3):min(ny, y + 4), max(0, x - 3):min(nx, x + 4)] = True
    y, x = np.ogrid[:ny, :nx]
    cutoff = (x - nx / 2.0) ** 2 + (y - ny / 2.0) ** 2 <= (0.4 * nx) ** 2
    return cutoff & ~excluded


@pytest.mark.parametrize("source_pair", [(41, 61), (61, 61)], indirect=True)
def test_dv_logits_read_yx_and_publish_gap_exclusion_with_one_y_flip(source_pair, tmp_path):
    hd, dv = source_pair
    height, width = hd["moment0"].shape[-2:]
    # This asymmetric vessel position distinguishes transpose and Y inversion,
    # including the square case where dimensions cannot reveal axis mistakes.
    vessel_y, vessel_x = 15, 35
    logits = np.full((height, width), -2, dtype=np.float32)
    logits[vessel_y, vessel_x] = 0.25
    dv.create_dataset(DV_RETINAL_VESSEL_LOGITS_PATH, data=logits)
    dv.create_dataset(DV_CHOROID_VESSEL_MASK_PATH, data=np.zeros_like(logits))
    ctx = _ctx(hd, dv)
    native = ctx.inputs.dv.as_dopplerview().retinal_vessel_mask()
    np.testing.assert_array_equal(native, logits > 0.1)
    assert native[vessel_y, vessel_x]
    assert not native[vessel_x, vessel_y]
    assert not native[height - 1 - vessel_y, vessel_x]

    metrics = run_topology_core(ctx)
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as output:
        PipelineH5Output(output).write_many(metrics)
        gaps = output[schema.segmentation.gaps]
        assert gaps.shape == (width, height)
        # Native (y=15, x=35) becomes EF (x=35, y=height-1-15).
        assert not gaps[vessel_x, height - 1 - vessel_y]
        assert gaps[vessel_x, vessel_y]  # A second/missing flip would exclude here.
        assert gaps[vessel_y, height - 1 - vessel_x]  # A transpose would exclude here.


@pytest.mark.parametrize("swapped", [False, True])
def test_publishes_cutoff_gaps_with_dilated_vessels_and_aligns_all_sources(source_pair, tmp_path, swapped):
    hd, dv = source_pair
    native_shape = hd["moment0"].shape[-2:]
    dv_shape = native_shape[::-1] if swapped else native_shape
    if swapped:
        for name in ("artery_mask", "vein_mask"):
            path = f"segmentation/Retina/{name}"
            del dv[path]
            dv.create_dataset(path, data=np.zeros(dv_shape, bool))
    logits, choroid, disc = _write_capillary_sources(dv, dv_shape)
    retinal, choroid, disc = logits > 0.1, choroid != 0, disc != 0
    if swapped:
        retinal, choroid, disc = retinal.T, choroid.T, disc.T
    expected = _reference_zones(retinal, choroid, disc)
    below, equal, above = [(x, 20) if swapped else (20, x) for x in (20, 28, 36)]
    assert expected[below]  # Positive logits below the threshold are included.
    assert expected[equal]  # Exactly 0.1 is included (strict greater-than).
    assert not expected[above]  # The next representable value above 0.1 is excluded.

    ctx = _ctx(hd, dv)
    metrics = run_topology_core(ctx)
    np.testing.assert_array_equal(capillary_zones(ctx), expected)
    assert not capillary_zones(ctx).flags.writeable
    path = EyeFlowOutputPaths.active().segmentation.gaps
    assert path == "Segmentation/Gaps/Mask"
    with h5py.File(tmp_path / "out.h5", "w") as output:
        PipelineH5Output(output).write_many(metrics)
        assert "Segmentation/Capillary" not in output
        dataset = output[path]
        np.testing.assert_array_equal(dataset[:], np.flip(expected, axis=0).T)
        assert dataset.dtype == np.uint8
        assert set(np.unique(dataset[:])) == {0, 1}
        assert dataset.attrs["original_class"] == "bool"
        assert list(dataset.attrs["dimDesc"]) == ["x", "y"]
        assert dataset.attrs["image_origin"] == "lower_left"
        assert dataset.attrs["legacy_fallback"]
        np.testing.assert_array_equal(
            dataset[:], output[EyeFlowOutputPaths.active().segmentation.gaps_legacy][:],
        )
        assert "Segmentation/Gaps/value" not in output
        assert EyeFlowOutputPaths.active().segmentation.gaps_from_m0ff not in output
        dataset = output[EyeFlowOutputPaths.active().segmentation.gaps_legacy]
        assert dataset.attrs["retinal_vessel_binarization"] == "logits > 0.1"
        assert dataset.attrs["cutoff_disc_diameter_width_fraction"] == 0.8
        assert dataset.attrs["cutoff_disc_diameter_pixels"] == pytest.approx(0.8 * native_shape[1])
        assert dataset.attrs["vessel_dilation_radius_pixels"] == 3
        assert dataset.attrs["vessel_dilation_iterations"] == 1
        assert dataset.attrs["vessel_dilation_footprint"] == "7x7_square"
        assert dataset.attrs["vessel_dilation_distance_metric"] == "chebyshev"
        assert dataset.attrs["vessel_dilation_applies_to"] == "retinal_and_choroid_vessels_only"
        assert dataset.attrs["optic_disc_mask_source"] == "dopplerview_segmentation"


@pytest.mark.parametrize("geometry_only", [False, True])
def test_missing_disc_mask_excludes_centered_physical_circle(source_pair, geometry_only):
    hd, dv = source_pair
    shape = hd["moment0"].shape[-2:]
    dv.create_dataset(DV_RETINAL_VESSEL_LOGITS_PATH, data=np.zeros(shape))
    dv.create_dataset(DV_CHOROID_VESSEL_MASK_PATH, data=np.zeros(shape))
    if geometry_only:
        dv.create_dataset("segmentation/OpticDisc/center", data=[1.0, 1.0])
        dv.create_dataset("segmentation/OpticDisc/width", data=2.0)
        dv.create_dataset("segmentation/OpticDisc/height", data=2.0)
    y, x = np.ogrid[: shape[0], : shape[1]]
    distance_squared = (x - shape[1] / 2.0) ** 2 + (y - shape[0] / 2.0) ** 2
    fallback_disc = distance_squared <= (1.91 / 0.2 / 2) ** 2
    expected = (distance_squared <= (0.4 * shape[1]) ** 2) & ~fallback_disc
    ctx = _ctx(hd, dv)

    metrics = run_topology_core(ctx)

    np.testing.assert_array_equal(capillary_zones(ctx), expected)
    path = EyeFlowOutputPaths.active().segmentation.gaps_legacy
    assert metrics[path][1]["optic_disc_fallback_diameter_mm"] == 1.91
    if not geometry_only:
        disc_mask = metrics[EyeFlowOutputPaths.active().segmentation.optic_disc.mask][0]
        np.testing.assert_array_equal(disc_mask, np.flip(fallback_disc, axis=0).T)


@pytest.mark.parametrize("missing", ["retinal", "choroid", "both"])
def test_missing_vessel_sources_omit_mask_and_warn_without_losing_topology(source_pair, missing):
    hd, dv = source_pair
    shape = hd["moment0"].shape[-2:]
    if missing != "retinal" and missing != "both":
        dv.create_dataset(DV_RETINAL_VESSEL_LOGITS_PATH, data=np.zeros(shape))
    if missing != "choroid" and missing != "both":
        dv.create_dataset(DV_CHOROID_VESSEL_MASK_PATH, data=np.zeros(shape))
    ctx = _ctx(hd, dv)
    with patch("pipelines.topology_core.runner.Logger.log_warning") as warning:
        metrics = run_topology_core(ctx)
    path = EyeFlowOutputPaths.active().segmentation.gaps
    assert path not in metrics
    assert EyeFlowOutputPaths.active().segmentation.artery.mask in metrics
    assert any("Omitting legacy capillary analysis zones" in call.args[0] for call in warning.call_args_list)
    assert any("moment0ff/M0FF is missing" in call.args[0] for call in warning.call_args_list)
    with pytest.raises(RuntimeError, match="unavailable"):
        capillary_zones(ctx)


@pytest.mark.parametrize("path", [DV_RETINAL_VESSEL_LOGITS_PATH, DV_CHOROID_VESSEL_MASK_PATH])
@pytest.mark.parametrize("defect", ["shape", "rank", "nan", "string"])
def test_invalid_capillary_inputs_fail_at_the_source_boundary(source_pair, path, defect):
    hd, dv = source_pair
    shape = hd["moment0"].shape[-2:]
    values = {
        "shape": np.zeros((4, 5)),
        "rank": np.zeros((1, *shape)),
        "nan": np.full(shape, np.nan),
        "string": np.full(shape, b"mask"),
    }[defect]
    dv.create_dataset(path, data=values)
    with pytest.raises((TypeError, ValueError), match="shape|2-D|finite|numeric"):
        run_topology_core(_ctx(hd, dv))


def test_capillary_zones_are_owned_by_existing_hidden_topology_pipeline():
    available, missing = load_pipeline_catalog()
    assert "topology_core" not in {item.name for item in missing}
    descriptors = {item.name: item for item in available}
    assert "capillary_analysis_zones" in descriptors["topology_core"].dag_produces
    assert not any(item.name.startswith("capillary") for item in [*available, *missing])
    assert "topology_core" in PipelineDAG(available).dependencies_of("velocity_analysis")


@pytest.mark.parametrize("source_name", ["moment0ff", "M0FF"])
@pytest.mark.parametrize("with_legacy", [False, True])
def test_m0ff_vessel_union_gap_intersection_persistence_and_spectral_consumer(
    source_pair, tmp_path, source_name, with_legacy,
):
    hd, dv = source_pair
    shape = hd["moment0"].shape[-2:]
    if with_legacy:
        _write_capillary_sources(dv, shape)
    image = np.full(shape, 2, dtype=np.float32)
    image[-12:, -15:] = 100
    hd.create_dataset(source_name, data=np.stack([image, image, image]))
    ctx = _ctx(hd, dv)
    metrics = run_topology_core(ctx)
    paths = EyeFlowOutputPaths.active()
    segmentation = paths.segmentation
    combined = capillary_zones(ctx)
    assert not combined.flags.writeable
    assert not combined[0, 0]  # The OR must not reintroduce outer-disc pixels.
    y, x = np.indices(shape)
    if with_legacy:
        disc = np.zeros(shape, dtype=bool)
        disc[12, 34] = True
        disc_y, disc_x = 12, 34
    else:
        disc = ((x - shape[1] / 2) ** 2 + (y - shape[0] / 2) ** 2) <= (1.91 / 0.2 / 2) ** 2
        disc_y, disc_x = shape[0] // 2, shape[1] // 2
    assert not combined[disc_y, disc_x]
    allowed = _reference_zones(np.zeros(shape, bool), np.zeros(shape, bool), disc)
    saved_allowed = np.flip(allowed, axis=0).T
    with h5py.File(tmp_path / "out.h5", "w") as output:
        PipelineH5Output(output).write_many(metrics)
        new = output[segmentation.gaps_from_m0ff][:].astype(bool)
        saved = output[segmentation.gaps][:].astype(bool)
        assert new[0, shape[0] - 1]  # The intensity mask itself has no disc constraints.
        assert new[disc_x, shape[0] - 1 - disc_y]
        if with_legacy:
            legacy = output[segmentation.gaps_legacy][:].astype(bool)
            assert np.any(new & ~legacy)
            assert np.any(legacy & ~new)
            np.testing.assert_array_equal(saved, (legacy & new) & saved_allowed)
            assert np.any(saved_allowed & new & ~legacy)
            assert np.any(saved_allowed & legacy & ~new)
            assert saved.sum() < legacy.sum()
            assert saved.sum() < (new & saved_allowed).sum()
            assert not np.any(saved & ~legacy)
            assert not np.any(saved & ~new)
        else:
            assert segmentation.gaps_legacy not in output
            np.testing.assert_array_equal(saved, new & saved_allowed)
        np.testing.assert_array_equal(saved, np.flip(combined, axis=0).T)
        assert "Segmentation/Gaps/value" not in output
        attrs = output[segmentation.gaps_from_m0ff].attrs
        assert attrs["source_path"] == f"/{source_name}"
        assert attrs["image_origin"] == "lower_left"
        assert attrs["dilation_radius_pixels"] == 3
        assert attrs["dilation_applies_to"] == "bright_vessels"
        assert attrs["dilation_footprint"] == "7x7_square"
        assert not output[segmentation.gaps].attrs["legacy_fallback"]
        assert output[segmentation.gaps].attrs["combination_operation"] == "OR"
        assert output[segmentation.gaps].attrs["equivalent_gap_operation"] == "AND"
        assert output[segmentation.gaps].attrs["cutoff_disc_diameter_width_fraction"] == 0.8

    frame = np.arange(np.prod(shape), dtype=np.float32).reshape(shape)
    spectra = frame[None, None]
    power, attrs = preprocess_spectrograms(spectra, combined)[paths.spectrograms_masked]
    assert power[0, 0] == np.median(frame[combined])
    assert attrs["mask_source_path"] == "/Segmentation/Gaps/Mask"


def test_m0ff_shape_must_match_native_hd_axes(source_pair):
    hd, dv = source_pair
    shape = hd["moment0"].shape[-2:]
    hd.create_dataset("moment0ff", data=np.ones((2, *shape[::-1])))
    with pytest.raises(ValueError, match="moment0ff spatial shape"):
        run_topology_core(_ctx(hd, dv))


@pytest.mark.parametrize("empty", [False, True])
def test_shared_topology_publishes_gap_support_without_running_spectral_analysis(
    source_pair, tmp_path, empty,
):
    hd, dv = source_pair
    shape = hd["moment0"].shape[-2:]
    logits = np.ones(shape) if empty else np.zeros(shape)
    dv.create_dataset(DV_RETINAL_VESSEL_LOGITS_PATH, data=logits)
    dv.create_dataset(DV_CHOROID_VESSEL_MASK_PATH, data=np.zeros(shape))
    yy, xx = np.indices(shape)
    frame = (100 * yy + xx).astype(np.float32)
    hd.create_dataset("spectrograms", data=np.stack([
        np.stack([frame, frame + 5]), np.stack([frame + 10, frame + 15]),
    ]))
    metrics = run_topology_core(_ctx(hd, dv))
    schema = EyeFlowOutputPaths.active()
    with h5py.File(tmp_path / "out.h5", "w") as output:
        PipelineH5Output(output).write_many(metrics)
        assert schema.segmentation.artery.mask in output
        assert schema.segmentation.gaps in output
        assert schema.spectrograms_masked not in output
        assert bool(np.any(output[schema.segmentation.gaps][:])) is not empty
