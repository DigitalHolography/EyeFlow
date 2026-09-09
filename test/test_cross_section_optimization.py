"""Regression checks for geometry, sparse storage and bounded cross sections."""

import importlib
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pytest
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
cs = importlib.import_module(
    "calculations.blood_flow_velocity.cross_section.generate_cross_section_signals"
)
from calculations.blood_flow_velocity.cross_section.gpu_cross_section import (
    measure_cross_section_gpu,
)
from calculations.blood_flow_velocity.cross_section.profile_processing import (
    process_velocity_profiles,
)
from calculations.blood_flow_velocity.cross_section.segment_array import SegmentArray
from calculations.blood_flow_velocity.cross_section.segment_geometry import (
    SegmentRingSettings,
    annulus_mask,
    image_half_diagonal,
    normalized_radius_squared,
)
from pipelines.waveform_velocity.segment_maps import interpolate_velocity_maps_per_beat


@pytest.fixture(autouse=True)
def cpu(monkeypatch):
    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "cpu")
    monkeypatch.setattr(cs, "_GPU_FAILED", False)


@pytest.mark.parametrize("shape,center", [((101, 101), (50, 50)), ((81, 151), (30, 45))])
def test_radii_match_annulus_units_and_mask_tilt(shape, center):
    r = cs._normalized_radius_grid(shape, center)
    squared = normalized_radius_squared(shape, center)
    np.testing.assert_allclose(r * r, squared, rtol=2e-7)
    np.testing.assert_array_equal(
        (squared > 0.2**2) & (squared <= 0.7**2), annulus_mask(shape, center, 0.2, 0.7)
    )
    section = annulus_mask(shape, center, 0.3, 0.6)
    vessel = np.zeros(shape, bool)
    vessel[center[1] - 1 : center[1] + 2, center[0] + 1 :] = True
    geometry = cs._circle_tilt_geometry(vessel & section, section, center)
    assert geometry is not None
    assert abs(geometry.tilt_angle) < 1e-6
    branches = SimpleNamespace(labels=vessel.astype(np.int32), branch_ids=np.array([1]))
    prepared = cs._prepare_segments(np.array([section]), branches, center)
    assert abs(prepared[0].tilt_angle_mask) < 1e-6
    pixel_radius = image_half_diagonal(*shape) * r[center[1], center[0] + 10]
    assert pixel_radius == pytest.approx(10, rel=1e-6)


def test_largest_component_defines_centroid():
    mask = np.zeros((20, 20), bool)
    mask[0, 0] = True
    mask[10:13, 10:13] = True
    assert cs._centroid_xy(mask) == (11, 11)


def test_missing_intersections_allocate_no_payload():
    buffers = cs._CrossSectionBuffers.allocate(frame_count=1000, ring_count=10, branch_count=20)
    maps = buffers.velocity_maps_per_segment
    assert maps.shape == (10, 20, 1000, 181, 181)
    assert maps.nbytes == 0
    # Reading an absent movie broadcasts a single NaN; it does not reserve a movie.
    missing = maps[2, 4]
    assert missing.strides == (0, 0, 0)
    assert maps.nbytes == 0
    for name in ("velocity", "velocity_profiles", "segment_masks", "rotated_mean_images"):
        assert getattr(buffers, name).nbytes == 0


def test_sparse_profiles_fit_only_present_segments_without_densifying():
    values = SegmentArray((3, 4, 2, 9))
    values[1, 2] = np.maximum(0, 16 - (np.arange(9) - 4) ** 2)
    with patch.object(SegmentArray, "__array__", side_effect=AssertionError("dense conversion")):
        result = process_velocity_profiles(
            values, pixel_size_mm=0.01, velocity_profile_threshold=0.5
        )
    assert result.centered_velocity.segment_indexes == ((1, 2),)
    assert np.all(np.isnan(result.centered_velocity[0, 0]))


def test_sparse_maps_export_matches_legacy_dense_input():
    maps = SegmentArray((2, 3, 5, 3, 3))
    maps[0, 2] = np.arange(45, dtype=np.float32).reshape(5, 3, 3)
    dense = np.asarray(maps)
    expected = interpolate_velocity_maps_per_beat(dense, [0, 2, 4])
    with patch.object(SegmentArray, "__array__", side_effect=AssertionError("dense input")):
        actual = interpolate_velocity_maps_per_beat(maps, [0, 2, 4])
    np.testing.assert_allclose(actual, expected, equal_nan=True)


@pytest.mark.parametrize(
    "field,value",
    [
        ("pixel_size_mm", 0),
        ("pixel_size_mm", np.nan),
        ("velocity_profile_threshold", 1),
        ("velocity_profile_threshold", np.inf),
        ("working_memory_mb", -1),
        ("submask_size_percentile_kept", np.nan),
    ],
)
def test_invalid_settings_are_rejected(field, value):
    with pytest.raises(ValueError):
        replace(cs.CrossSectionSignalSettings(True, 0.5, False, 0.01), **{field: value})


def test_missing_frames_remain_missing_and_zero_flow_remains_zero():
    cube = np.zeros((2, 7, 7), np.float32)
    cube[0] = np.nan
    raw, safe, *_ = cs._frame_velocities(cube, 0, 0, 6)
    np.testing.assert_array_equal(raw, [np.nan, 0])
    np.testing.assert_array_equal(safe, [np.nan, 0])


def test_empty_result_has_zero_branches_everywhere():
    branches = SimpleNamespace(branch_ids=np.array([], np.int32), labels=np.zeros((3, 3), np.int32))
    result = cs._empty_result(
        np.zeros((2, 3, 3)),
        np.zeros((3, 3), bool),
        SegmentRingSettings(0.1, 0.9, 0.2, 2),
        branches,
        substack_side_pixels=0,
        profile_pixel_size_mm=0,
    )
    assert result.velocity.shape == (2, 0, 2)
    assert result.safe_velocity.shape == (2, 0, 2)
    assert result.velocity_profiles.shape[1] == 0


def test_flat_and_upward_fits_rejected_and_parabolic_roots_scale():
    settings = cs.CrossSectionSignalSettings(True, 0.5, False, 0.01)
    x = np.arange(11) - 5
    assert cs._hydrodynamic_limits(np.ones(11), settings, 0.01) is None
    assert cs._half_height_roots(x * x + 10, np.arange(11), 5, 10, 0.01) is None
    profile = 25.0 - x * x
    for scale in (1e-12, 0.01, 1e3):
        assert cs._hydrodynamic_limits(profile, settings, scale) == (2, 8)


def _window_fixture():
    rng = np.random.default_rng(24)
    velocity = rng.uniform(1, 5, (7, 13, 13)).astype(np.float32)
    velocity[2, 4:6, 4:7] = np.nan
    velocity[5] = np.nan
    mask = np.zeros((13, 13), bool)
    mask[3:10, 5:8] = True
    yy, xx = np.nonzero(mask)
    seg = cs._SegmentGeometry(0, 0, (6, 6), yy, xx, 0.0)
    return velocity, cs._CrossSectionWork(seg, (2, 11, 2, 11))


def test_temporal_batches_match_full_window_with_global_fit():
    velocity, work = _window_fixture()
    settings = cs.CrossSectionSignalSettings(True, 0.5, True, 0.01)
    full = cs._CrossSectionBuffers.allocate(frame_count=7, ring_count=1, branch_count=1)
    batched = cs._CrossSectionBuffers.allocate(frame_count=7, ring_count=1, branch_count=1)
    cs._measure_windowed_work(full, velocity, work, (6, 6), settings, 9)
    budget = (cs._estimated_work_bytes(2, 9) + 1) / 1024**2
    with (
        patch.object(cs, "_extract_work", wraps=cs._extract_work) as extract,
        patch.object(cs, "_cross_section_limits", wraps=cs._cross_section_limits) as fit,
    ):
        cs._measure_windowed_work(
            batched, velocity, work, (6, 6), replace(settings, working_memory_mb=budget), 9
        )
    assert fit.call_count == 1
    assert all(call.args[4] - call.args[3] <= 2 for call in extract.call_args_list)
    for name in (
        "velocity",
        "safe_velocity",
        "velocity_maps_per_segment",
        "velocity_profiles",
        "transverse_velocity_profiles_masked",
        "rotated_mean_images_masked",
    ):
        np.testing.assert_allclose(
            getattr(batched, name)[0, 0],
            getattr(full, name)[0, 0],
            rtol=3e-6,
            atol=2e-6,
            equal_nan=True,
        )
    np.testing.assert_array_equal(
        batched.profile_integration_limits_pixels, full.profile_integration_limits_pixels
    )


def test_worker_count_obeys_memory_and_runtime_caps():
    per_work = cs._estimated_work_bytes(4, 9)
    with patch.object(cs, "cap_parallel_jobs", return_value=8):
        assert (
            cs._cross_section_worker_count(
                10, frame_count=4, side_pixels=9, memory_mb=2.1 * per_work / 1024**2
            )
            == 2
        )
        assert (
            cs._cross_section_worker_count(
                10, frame_count=400, side_pixels=9, memory_mb=2.1 * per_work / 1024**2
            )
            == 1
        )
    with pytest.raises(MemoryError):
        cs._cross_section_worker_count(1, frame_count=4, side_pixels=9, memory_mb=0.01)


def test_windows_are_extracted_only_when_measurement_starts():
    velocity, work = _window_fixture()
    segments = (work.segment, replace(work.segment, circle_index=1))
    buffers = cs._CrossSectionBuffers.allocate(frame_count=7, ring_count=2, branch_count=1)
    events = []
    original_extract = cs._extract_work
    original_measure = cs._cross_section_velocity_from_substack

    def extract(*args, **kwargs):
        events.append("extract")
        return original_extract(*args, **kwargs)

    def measure(*args, **kwargs):
        events.append("measure")
        return original_measure(*args, **kwargs)

    with (
        patch.object(cs, "_extract_work", side_effect=extract),
        patch.object(cs, "_cross_section_velocity_from_substack", side_effect=measure),
        patch.object(cs, "cap_parallel_jobs", return_value=1),
    ):
        cs._fill_cross_section_buffers(
            buffers,
            velocity,
            None,
            SimpleNamespace(labels=np.zeros((13, 13))),
            (6, 6),
            cs.CrossSectionSignalSettings(False, 0.5, True, 0.01),
            9,
            segments=segments,
        )
    assert events == ["extract", "measure", "extract", "measure"]


class NumpyDevice:
    """Exercise device-path math and transfer boundaries without CUDA hardware."""

    def __init__(self):
        self.download_shapes = []

    def __getattr__(self, name):
        return getattr(np, name)

    def asnumpy(self, value):
        self.download_shapes.append(value.shape)
        return value.copy()


def test_device_path_matches_cpu_and_downloads_only_one_movie():
    velocity, work = _window_fixture()
    stack, mask = cs._extract_work(velocity, work, 9)
    settings = cs.CrossSectionSignalSettings(True, 0.5, True, 0.01)
    cpu_result = cs._cross_section_velocity_from_substack(
        stack,
        mask,
        (6, 6),
        (6, 6),
        0,
        settings,
        9,
        angle_override=31,
    )
    device = NumpyDevice()
    gpu_result = measure_cross_section_gpu(
        SimpleNamespace(cupy=device, ndimage=ndimage),
        stack,
        mask,
        (6, 6),
        (6, 6),
        0,
        settings,
        9,
        angle_override=31,
    )
    for kind in ("masked", "unmasked"):
        for field in ("raw", "safe_velocity", "transverse_profiles", "longitudinal_profiles"):
            np.testing.assert_allclose(
                getattr(getattr(gpu_result, kind), field),
                getattr(getattr(cpu_result, kind), field),
                rtol=3e-6,
                atol=2e-6,
                equal_nan=True,
            )
    assert gpu_result.limits == cpu_result.limits
    assert sum(len(shape) == 3 for shape in device.download_shapes) == 1
    assert gpu_result.masked.rotated_stack is None


def test_gpu_failure_is_reported_once_and_explicit_mode_raises(monkeypatch):
    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "auto")
    with pytest.warns(RuntimeWarning, match="using CPU"):
        cs._disable_cross_section_gpu(RuntimeError("device failed"))
    assert cs._cross_section_backend() is None
    cs._disable_cross_section_gpu(RuntimeError("device failed again"))
    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "cupy")
    with pytest.raises(RuntimeError, match="Explicit CuPy"):
        cs._disable_cross_section_gpu(RuntimeError("device failed"))


def test_failed_device_dispatch_retries_on_cpu(monkeypatch):
    velocity, work = _window_fixture()
    stack, mask = cs._extract_work(velocity, work, 9)
    settings = cs.CrossSectionSignalSettings(False, 0.5, True, 0.01)
    expected = cs._cross_section_velocity_from_substack(
        stack,
        mask,
        (6, 6),
        (6, 6),
        0,
        settings,
        9,
    )

    def fail(*args, **kwargs):
        raise RuntimeError("device allocation failed")

    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "auto")
    monkeypatch.setattr(
        cs,
        "optional_cupy_backend",
        lambda: SimpleNamespace(
            cupy=NumpyDevice(),
            ndimage=SimpleNamespace(zoom=fail),
        ),
    )
    with pytest.warns(RuntimeWarning, match="using CPU"):
        actual = cs._cross_section_velocity_from_substack(
            stack,
            mask,
            (6, 6),
            (6, 6),
            0,
            settings,
            9,
        )
    np.testing.assert_allclose(actual.masked.raw, expected.masked.raw, equal_nan=True)
    assert cs._cross_section_backend() is None


def test_cuda_parity_when_available(monkeypatch):
    cupy = pytest.importorskip("cupy")
    try:
        if cupy.cuda.runtime.getDeviceCount() == 0:
            pytest.skip("No CUDA device")
    except cupy.cuda.runtime.CUDARuntimeError:
        pytest.skip("CUDA runtime unavailable")
    from cupyx.scipy import ndimage as gpu_ndi

    velocity, work = _window_fixture()
    stack, mask = cs._extract_work(velocity, work, 9)
    settings = cs.CrossSectionSignalSettings(False, 0.5, True, 0.01)
    expected = cs._cross_section_velocity_from_substack(
        stack, mask, (6, 6), (6, 6), 0, settings, 9, angle_override=31
    )
    actual = measure_cross_section_gpu(
        SimpleNamespace(cupy=cupy, ndimage=gpu_ndi),
        stack,
        mask,
        (6, 6),
        (6, 6),
        0,
        settings,
        9,
        angle_override=31,
    )
    np.testing.assert_allclose(
        actual.masked.raw, expected.masked.raw, rtol=1e-5, atol=1e-5, equal_nan=True
    )
    # Exercise the GPU helpers used in the global-mean pass as well as the
    # device-resident path, with explicit mode preventing silent CPU fallback.
    budget = (cs._estimated_work_bytes(2, 9) + 1) / 1024**2
    small = replace(settings, working_memory_mb=budget)
    cpu_buffers = cs._CrossSectionBuffers.allocate(frame_count=7, ring_count=1, branch_count=1)
    gpu_buffers = cs._CrossSectionBuffers.allocate(frame_count=7, ring_count=1, branch_count=1)
    cs._measure_windowed_work(cpu_buffers, velocity, work, (6, 6), small, 9, angle_override=31)
    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "cupy")
    monkeypatch.setattr(
        cs, "optional_cupy_backend", lambda: SimpleNamespace(cupy=cupy, ndimage=gpu_ndi)
    )
    cs._measure_windowed_work(gpu_buffers, velocity, work, (6, 6), small, 9, angle_override=31)
    np.testing.assert_allclose(
        gpu_buffers.velocity[0, 0], cpu_buffers.velocity[0, 0], rtol=1e-5, atol=1e-5, equal_nan=True
    )
    assert not cs._GPU_FAILED
