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


@pytest.fixture
def source_pair(tmp_path):
    shape = (9, 14)
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
    logits[0, 0] = 0.25  # Positive logits below 0.5 must still be excluded.
    logits[0, 1] = 0.0
    choroid = np.zeros(shape, dtype=np.int32)
    choroid[-1, -2] = 7
    disc = np.zeros(shape, dtype=np.uint8)
    disc[2, 4] = 255
    dv.create_dataset(DV_RETINAL_VESSEL_LOGITS_PATH, data=logits)
    dv.create_dataset(DV_CHOROID_VESSEL_MASK_PATH, data=choroid)
    dv.create_dataset("segmentation/OpticDisc/mask", data=disc)
    return logits, choroid, disc


@pytest.mark.parametrize("swapped", [False, True])
def test_publishes_inverse_union_and_aligns_all_sources(source_pair, tmp_path, swapped):
    hd, dv = source_pair
    native_shape = hd["moment0"].shape[-2:]
    dv_shape = native_shape[::-1] if swapped else native_shape
    if swapped:
        for name in ("artery_mask", "vein_mask"):
            path = f"segmentation/Retina/{name}"
            del dv[path]
            dv.create_dataset(path, data=np.zeros(dv_shape, bool))
    logits, choroid, disc = _write_capillary_sources(dv, dv_shape)
    expected = ~((logits > 0) | (choroid != 0) | (disc != 0))
    if swapped:
        expected = expected.T

    ctx = _ctx(hd, dv)
    metrics = run_topology_core(ctx)
    np.testing.assert_array_equal(capillary_zones(ctx), expected)
    assert not capillary_zones(ctx).flags.writeable
    path = EyeFlowOutputPaths.active().segmentation.capillary.analysis_zones
    assert path == "Segmentation/Capillary/analysis_zones"
    with h5py.File(tmp_path / "out.h5", "w") as output:
        PipelineH5Output(output).write_many(metrics)
        dataset = output[path]
        np.testing.assert_array_equal(dataset[:], np.flip(expected, axis=0).T)
        assert dataset.dtype == np.uint8
        assert set(np.unique(dataset[:])) == {0, 1}
        assert dataset.attrs["original_class"] == "bool"
        assert list(dataset.attrs["dimDesc"]) == ["x", "y"]
        assert dataset.attrs["image_origin"] == "lower_left"
        assert dataset.attrs["retinal_vessel_binarization"] == "logits > 0"
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
    expected = (x - shape[1] / 2.0) ** 2 + (y - shape[0] / 2.0) ** 2 > (1.91 / 0.2 / 2) ** 2
    ctx = _ctx(hd, dv)

    metrics = run_topology_core(ctx)

    np.testing.assert_array_equal(capillary_zones(ctx), expected)
    path = EyeFlowOutputPaths.active().segmentation.capillary.analysis_zones
    assert metrics[path][1]["optic_disc_fallback_diameter_mm"] == 1.91
    if not geometry_only:
        disc_mask = metrics[EyeFlowOutputPaths.active().segmentation.optic_disc.mask][0]
        np.testing.assert_array_equal(disc_mask, np.flip(~expected, axis=0).T)


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
    path = EyeFlowOutputPaths.active().segmentation.capillary.analysis_zones
    assert path not in metrics
    assert EyeFlowOutputPaths.active().segmentation.artery.mask in metrics
    warning.assert_called_once()
    assert "Omitting capillary analysis zones" in warning.call_args.args[0]
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
