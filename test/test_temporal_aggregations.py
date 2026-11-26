"""Frame-axis temporal statistics and HoloDoppler artifact routing."""

import h5py
import numpy as np
import pytest
from PIL import Image

from app_settings import UserSettings
from calculations import temporal_statistics
from calculations.compute_backend import optional_cupy_backend
from input_output.runs.layout import HoloRunLayout
from pipeline_engine import PipelineDAG, execute_run, resolve_run_spec
from pipelines import load_pipeline_catalog


def test_frame_statistics_use_finite_samples_and_normalized_magnitudes(monkeypatch):
    video = np.array([
        [[0., -1., np.nan], [5., 1e9, 1.]],
        [[0., 1., np.inf], [5., 1e9 + 1, 3.]],
        [[0., 2., -np.inf], [5., 1e9 + 2, np.nan]],
    ])
    monkeypatch.setattr(temporal_statistics, "_BLOCK_BYTES", 8 * 6 * 2)
    stats = temporal_statistics.temporal_image_statistics(video)

    np.testing.assert_allclose(stats["avg"], [[0., 2 / 3, np.nan], [5., 1e9 + 1, 2.]])
    np.testing.assert_allclose(stats["std"], [[0., np.sqrt(7 / 3), np.nan], [0., 1., np.sqrt(2)]])
    probabilities = np.array([0.25, 0.75])
    np.testing.assert_allclose(stats["entropy"], [
        [0., 1.5 * np.log(2), np.nan],
        [np.log(3), np.log(3), -np.sum(probabilities * np.log(probabilities))],
    ])


@pytest.mark.parametrize("scale,display", [
    ("per_frame", [[[0., 1 / 3], [2 / 3, 1.]], [[0., 1.], [1 / 3, 2 / 3]]]),
    ("global", [[[1 / 3, 4 / 9], [5 / 9, 2 / 3]], [[0., 1.], [1 / 3, 2 / 3]]]),
])
def test_matlab_summaries_reduce_flash_subtracted_display_frames(scale, display):
    video = np.array([[[100., 102.], [104., 106.]], [[1000., 1018.], [1006., 1012.]]])
    stats = temporal_statistics.temporal_image_statistics(video, display=scale)
    first, second = np.asarray(display)
    np.testing.assert_allclose(stats["avg"], (first + second) / 2)
    np.testing.assert_allclose(stats["std"], np.abs(first - second) / np.sqrt(2))
    np.testing.assert_array_equal(stats["entropy"], temporal_statistics.temporal_image_statistics(video)["entropy"])


def test_matlab_contrast_stretch_clips_outliers_and_blacks_out_invalid_pixels():
    frame = np.concatenate(([0, np.nan], np.ones(10_000), [2], np.full(10_000, 3), [4]))
    video = frame[None, None, :]
    stats = temporal_statistics.temporal_image_statistics(video, display="per_frame")

    expected = np.concatenate((np.zeros(10_002), [0.5], np.ones(10_001)))
    np.testing.assert_allclose(stats["avg"][0], expected)
    np.testing.assert_array_equal(stats["std"], np.zeros_like(stats["std"]))
    assert np.isnan(stats["entropy"][0, 1])


def test_shared_display_scale_uses_one_chunk_aligned_pass_with_a_final_partial_block(monkeypatch):
    data = np.array([
        [[0., 2., 4., 6.]],
        [[2., np.nan, 6., 10.]],
        [[2., 4., 8., 10.]],
        [[4., 8., np.nan, 16.]],
        [[np.nan, np.nan, np.nan, np.nan]],
    ], dtype=np.float32)
    reads = []

    class Video:
        shape = data.shape
        dtype = data.dtype
        chunks = (2, *data.shape[1:])

        def read_direct(self, buffer, *, source_sel):
            reads.append(source_sel)
            assert buffer.dtype == self.dtype
            np.copyto(buffer, data[source_sel])

    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "cpu")
    monkeypatch.setattr(temporal_statistics, "_BLOCK_BYTES", 3 * 4 * 8)
    stats = temporal_statistics.temporal_image_statistics(Video(), display="global")

    finite = np.isfinite(data)
    centered = data[:4].astype(np.float64) - np.nanmean(data[:4], axis=(1, 2), keepdims=True, dtype=np.float64)
    low, high = np.nanmin(centered), np.nanmax(centered)
    normalized = np.zeros(data.shape)
    normalized[:4] = np.where(finite[:4], (centered - low) / (high - low), 0.0)
    np.testing.assert_allclose(stats["avg"], normalized.mean(axis=0))
    np.testing.assert_allclose(stats["std"], normalized.std(axis=0, ddof=1))
    assert reads == [slice(0, 2), slice(2, 4), slice(4, 5)]


@pytest.mark.parametrize("scale", ["raw", "global", "per_frame"])
def test_cuda_matches_cpu_across_blocks_and_rgb_channels(monkeypatch, scale):
    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "auto")
    if optional_cupy_backend() is None:
        pytest.skip("CuPy/CUDA is unavailable")
    video = np.random.default_rng(23).standard_normal((5, 65, 65, 3)).astype(np.float32)
    video[0] = np.nan
    video[:, 0, 0, 0] = np.inf
    video[:, 0, 0, 1] = 0
    monkeypatch.setattr(temporal_statistics, "_BLOCK_BYTES", 2 * video[0].size * 8)
    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "cpu")
    expected = temporal_statistics.temporal_image_statistics(video, display=scale)
    monkeypatch.setenv("EYEFLOW_COMPUTE_BACKEND", "cupy")
    actual = temporal_statistics.temporal_image_statistics(video, display=scale)
    for name in expected:
        np.testing.assert_allclose(actual[name], expected[name], rtol=1e-12, atol=1e-13)


@pytest.mark.usefixtures("inline_output_service")
def test_temporal_pipeline_exports_one_png_per_video_and_statistic(tmp_path):
    layout = HoloRunLayout.from_holo(tmp_path / "scan.holo")
    layout.holo_path.touch()
    hd_path = layout.root_dir / "scan_HD/h5/scan_HD_output.h5"
    hd_path.parent.mkdir(parents=True)
    dv_path = layout.root_dir / "scan_DV/h5/scan_DV.h5"
    dv_path.parent.mkdir(parents=True)
    with h5py.File(dv_path, "w"):
        pass
    frame = np.arange(6, dtype=np.float32).reshape(2, 3)
    video = np.stack([frame, 2 * frame])
    rgb = np.stack([video, 0.5 * video[:, :, ::-1], 0.25 * video[:, ::-1]], axis=-1)
    with h5py.File(hd_path, "w") as hd:
        hd["moment0"] = video
        hd["M0"] = hd["moment0"]
        hd["M0FF"] = video
        hd["moment1"] = video
        hd["moment2"] = video
        hd["band_0_3000_6000"] = video
        hd["band_1_9000_18000"] = video
        hd["color_Doppler_image"] = rgb
        hd["directional_Doppler_image"] = rgb
        hd["averaged_spectrogram"] = video
        hd["spectrograms"] = rgb
        hd["registration"] = np.ones((2, 3, 3))

    available, _ = load_pipeline_catalog()
    velocity_plan = PipelineDAG(available).resolve_targets(("velocity_analysis",))
    assert "temporal_aggregations" not in velocity_plan.names
    spec = resolve_run_spec(
        input_paths=[layout.holo_path], target_names=["temporal_aggregations"],
        pipelines=available, settings=UserSettings(values={}),
    )
    result = execute_run(spec)
    assert result.succeeded, result.failures

    directory = layout.ef_dir / "png/temporal_aggregations"
    grayscale_datasets = {"moment0", "moment0ff", "moment1", "moment2", "band_0_3000_6000", "band_1_9000_18000"}
    rgb_datasets = {"color_Doppler_image", "directional_Doppler_image"}
    statistics = ("avg", "std", "entropy")
    assert {path.relative_to(directory).as_posix() for path in directory.rglob("*.png")} == {
        f"{folder}/scan_temporal_{statistic}.png"
        for folder in grayscale_datasets | rgb_datasets
        for statistic in statistics
    }
    expected_average = [
        [[0, 115, 134], [51, 89, 147], [102, 64, 159]],
        [[153, 191, 96], [204, 166, 108], [255, 140, 121]],
    ]
    for name in rgb_datasets:
        for statistic in statistics:
            with Image.open(directory / name / f"scan_temporal_{statistic}.png") as image:
                assert image.mode == "RGB"
                assert image.size == (3, 2)
                if statistic == "avg":
                    np.testing.assert_array_equal(image, expected_average)
    with Image.open(directory / "moment0/scan_temporal_avg.png") as image:
        assert image.mode == "L"
