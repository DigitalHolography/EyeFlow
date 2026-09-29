"""Numerical contracts for blood-volume-rate calculations."""

from __future__ import annotations

import tempfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from PIL import Image

from calculations.blood_volume_rate import (
    circular_lumen_flow,
    mask_derived_lumen_geometry,
    total_masked_edges_flow,
)
from calculations.math import nanmean_float32
from calculations.topology import AnnulusGeometry
from input_output.output_manager import OutputManager
from input_output.profile_datasets import _profile_dataset
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import DatasetValue
from pipelines.blood_volume_rate.outputs import (
    BLOOD_VOLUME_RATE_ENVELOPE_GRAY,
    LUMEN_DIAMETER_FIGURE_ASPECT_RATIO,
    _blood_volume_rate_figure,
    _lumen_diameter_distribution_figure,
    export_blood_volume_rate_signals,
    export_lumen_diameter_distributions,
    pack_gradient_edge_outputs,
    pack_mask_derived_outputs,
)


def test_circular_lumen_flow_supports_dynamic_and_static_geometry() -> None:
    velocity = np.asarray(
        [[[[2.0]]], [[[3.0]]]],
        dtype=np.float32,
    )
    static_diameter = np.asarray([[0.4]], dtype=np.float32)
    dynamic_diameter = np.asarray(
        [[[[0.2]]], [[[0.4]]]],
        dtype=np.float32,
    )

    np.testing.assert_allclose(
        circular_lumen_flow(velocity, static_diameter),
        velocity * np.pi * static_diameter**2 / 4.0,
    )
    np.testing.assert_allclose(
        circular_lumen_flow(velocity, dynamic_diameter),
        velocity * np.pi * dynamic_diameter**2 / 4.0,
    )


def test_mask_geometry_and_signed_flow_keep_established_model() -> None:
    labels = np.asarray(
        [
            [0, 1, 1, 0],
            [0, 1, 1, 0],
            [0, 1, 1, 0],
            [0, 1, 1, 0],
        ],
        dtype=np.int32,
    )
    topology = SimpleNamespace(
        labels=labels,
        branch_ids=np.asarray([1], dtype=np.int32),
        annulus_masks=np.ones((1, 4, 4), dtype=bool),
        ring_settings=AnnulusGeometry(0.0, 0.5, 0.5, 1, 0.5),
    )
    diameters, areas, widths = mask_derived_lumen_geometry(
        (topology,),
        pixel_size_mm=0.1,
    )
    expected_diameter = 8 * 0.1 / widths[0]
    np.testing.assert_array_equal(areas[0], [[8]])
    np.testing.assert_allclose(diameters[0], [[expected_diameter]])

    velocity = np.full((8, 1, 1, 1), -2.0, dtype=np.float32)
    rate = circular_lumen_flow(velocity, diameters[0])
    expected_rate = -2.0 * np.pi / 4.0 * expected_diameter**2
    np.testing.assert_allclose(rate, expected_rate)
    np.testing.assert_allclose(total_masked_edges_flow(rate), expected_rate)


def test_total_masked_edges_flow_uses_centered_periodic_nine_point_window() -> None:
    rate = np.zeros((10, 1, 1, 1), dtype=np.float32)
    rate[0] = 9.0

    expected = np.ones((10, 1), dtype=np.float32)
    expected[5, 0] = 0.0
    np.testing.assert_array_equal(total_masked_edges_flow(rate), expected)


def test_output_packers_keep_paths_units_and_valid_provenance() -> None:
    prepared = object()
    velocity_topology = SimpleNamespace(
        valid_segments=np.ones((1, 1), dtype=bool),
        segment_centers_xy=np.asarray([[[2.0, 3.0]]], dtype=np.float32),
        profile_rotation_degrees=np.asarray([[17.0]], dtype=np.float32),
        prepared_topology=prepared,
    )
    profiles = np.broadcast_to(
        np.asarray([0.0, 1.0, 4.0, 9.0, 16.0, 25.0], dtype=np.float32),
        (1, 1, 3, 6),
    ).copy()
    velocity_segments = SimpleNamespace(
        labels=np.asarray([[1]], dtype=np.int32),
        branch_ids=np.asarray([1], dtype=np.int32),
        topology=velocity_topology,
        transverse_profiles_masked=profiles,
        profile_pixel_size_mm=0.02,
    )
    gradient_segments = SimpleNamespace(
        labels=velocity_segments.labels,
        branch_ids=velocity_segments.branch_ids,
        topology=SimpleNamespace(
            segment_centers_xy=np.asarray([[[2.0, 3.0]]], dtype=np.float32),
            profile_rotation_degrees=np.asarray([[17.0]], dtype=np.float32),
            prepared_topology=prepared,
        ),
    )
    cycle_boundaries = np.asarray([0, 2], dtype=np.int32)
    profile_dataset = _profile_dataset(
        profiles,
        cycle_boundaries,
        index_base=0,
        valid_segments=np.ones((1, 1), dtype=bool),
    )
    edge_shape = profile_dataset.data.shape[1:]
    edges = {}
    for vessel in ("Artery", "Vein"):
        root = (
            f"Processing/SpatialGradientMetrics/{vessel}/"
            "Transverse/Masked/tbkr"
        )
        edges[f"{root}/left_edge_index"] = DatasetValue(
            np.full(edge_shape, 0.5, dtype=np.float32)
        )
        edges[f"{root}/right_edge_index"] = DatasetValue(
            np.full(edge_shape, 4.5, dtype=np.float32)
        )
    gradient_outputs = pack_gradient_edge_outputs(
        velocity_segments,
        velocity_segments,
        SimpleNamespace(
            artery_segments=gradient_segments,
            vein_segments=gradient_segments,
            outputs=edges,
        ),
        cycle_boundaries,
        index_base=0,
    )

    topology = SimpleNamespace(
        labels=np.ones((4, 4), dtype=np.int32),
        branch_ids=np.asarray([1], dtype=np.int32),
        annulus_masks=np.ones((1, 4, 4), dtype=bool),
        ring_settings=AnnulusGeometry(0.0, 0.5, 0.5, 1, 0.5),
    )
    schema = EyeFlowOutputPaths.active()
    safe_velocity = np.full((8, 1, 1, 1), -2.0, dtype=np.float32)
    mask_outputs = pack_mask_derived_outputs(
        {"artery": topology, "vein": topology},
        {
            schema.artery_per_beat_safe.velocity_signal: safe_velocity,
            schema.vein_per_beat_safe.velocity_signal: safe_velocity,
        },
        pixel_size_mm=0.1,
    )
    outputs = {**gradient_outputs, **mask_outputs}

    expected_gradient_rate = circular_lumen_flow(
        nanmean_float32(profile_dataset.data, axis=0),
        np.full(edge_shape, (4.5 - 0.5) * 0.02, dtype=np.float32),
    )

    for vessel_paths in (
        schema.blood_volume_rate.artery,
        schema.blood_volume_rate.vein,
    ):
        for path in (
            vessel_paths.dynamic_edges,
            vessel_paths.static_edges,
            vessel_paths.masked_edges,
            vessel_paths.total_masked_edges,
        ):
            assert path in outputs
            assert outputs[path].attrs["unit"] == "mm^3/s"
        assert not gradient_outputs[vessel_paths.dynamic_edges].attrs[
            "source_velocity"
        ].startswith("/")
        np.testing.assert_allclose(
            gradient_outputs[vessel_paths.dynamic_edges].data,
            expected_gradient_rate,
        )
        assert (
            gradient_outputs[vessel_paths.dynamic_edges].attrs[
                "cross_section_model"
            ]
            == "circular_pi_diameter_squared_over_4"
        )
        assert (
            mask_outputs[vessel_paths.total_masked_edges].attrs["source"]
            == f"/{vessel_paths.masked_edges}"
        )


def test_lumen_diameter_distribution_uses_five_micron_bins_and_scaled_gaussian() -> None:
    values = np.asarray([10.0, 20.0, 30.0, 40.0, 50.0])

    fig = _lumen_diameter_distribution_figure(values)
    ax = fig.axes[0]

    assert ax.get_title() == ""
    assert ax.get_xlabel() == "Lumen diameter (µm)"
    assert ax.get_ylabel() == "Count"
    assert all(np.isclose(patch.get_width(), 5.0) for patch in ax.patches)
    assert all(np.allclose(patch.get_facecolor()[:3], (0.8, 0.8, 0.8)) for patch in ax.patches)
    assert all(patch.get_edgecolor()[:3] == (0.25, 0.25, 0.25) for patch in ax.patches)
    assert all(np.isclose(patch.get_linewidth(), 0.25) for patch in ax.patches)
    assert len(ax.lines) == 1
    gaussian = ax.lines[0]
    assert gaussian.get_color() == "black"
    assert gaussian.get_linestyle() == "--"
    assert np.isclose(gaussian.get_linewidth(), 0.625)
    peak_index = int(np.argmax(gaussian.get_ydata()))
    assert np.isclose(gaussian.get_xdata()[peak_index], np.median(values), atol=0.1)
    expected_peak = values.size * 5.0 / (np.std(values) * np.sqrt(2.0 * np.pi))
    assert np.isclose(gaussian.get_ydata()[peak_index], expected_peak, rtol=1e-4)
    legend_labels = [text.get_text() for text in ax.get_legend().get_texts()]
    assert legend_labels == [
        "Median: 30.0 µm",
        "SD: 14.1 µm",
    ]
    assert all(np.isclose(text.get_fontsize(), 8.0) for text in ax.get_legend().get_texts())
    assert not ax.get_legend().get_frame_on()
    assert np.isclose(ax.xaxis.label.get_fontsize(), 8.0)
    assert np.isclose(ax.yaxis.label.get_fontsize(), 8.0)
    figure_text = [*ax.get_xticklabels(), *ax.get_yticklabels(), *ax.get_legend().get_texts()]
    figure_text.extend((ax.xaxis.label, ax.yaxis.label))
    assert all(text.get_fontfamily() == ["Times New Roman"] for text in figure_text)
    assert all(spine.get_visible() for spine in ax.spines.values())
    assert all(spine.get_edgecolor()[:3] == (0.0, 0.0, 0.0) for spine in ax.spines.values())
    assert all(np.isclose(spine.get_linewidth(), 0.4) for spine in ax.spines.values())
    width, height = fig.get_size_inches()
    assert np.isclose(width / height, LUMEN_DIAMETER_FIGURE_ASPECT_RATIO)


def test_lumen_diameter_distributions_export_png_and_eps_for_both_vessels() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        output = OutputManager.from_holo(
            Path(temp_dir) / "sample.holo",
            output_root=Path(temp_dir),
        )
        paths = export_lumen_diameter_distributions(
            output,
            np.asarray([[1.0, 2.0], [np.nan, 3.0]], dtype=np.float32),
            np.asarray([[2.0, 3.0], [4.0, 5.0]], dtype=np.float32),
            pixel_pitch_m=10e-6,
        )

        assert {path.relative_to(output.layout.ef_dir).as_posix() for path in paths} == {
            "png/lumen_diameter/artery_lumen_diameter_distribution.png",
            "eps/lumen_diameter/artery_lumen_diameter_distribution.eps",
            "png/lumen_diameter/vein_lumen_diameter_distribution.png",
            "eps/lumen_diameter/vein_lumen_diameter_distribution.eps",
        }
        for path in paths:
            assert path.is_file()
            if path.suffix == ".png":
                with Image.open(path) as image:
                    assert image.format == "PNG"
                    assert image.width >= 1000
                    assert np.isclose(
                        image.width / image.height,
                        LUMEN_DIAMETER_FIGURE_ASPECT_RATIO,
                        rtol=0.01,
                    )
            else:
                assert path.read_bytes().startswith(b"%!PS-Adobe")


def test_blood_volume_rate_figure_plots_median_and_one_sd_over_beats() -> None:
    values = np.asarray(
        [
            [1.0, 2.0, 3.0],
            [4.0, 4.0, 4.0],
            [1.0, np.nan, 5.0],
        ]
    )

    fig = _blood_volume_rate_figure(values)
    ax = fig.axes[0]

    assert ax.get_title() == ""
    assert ax.get_xlabel() == r"Cardiac Phase $t/T$"
    assert ax.get_ylabel() == r"$Q(t)$ (mm³/s)"
    assert ax.get_xlim() == (0.0, 1.0)
    assert len(ax.lines) == 2
    zero_line, median_line = ax.lines
    np.testing.assert_allclose(zero_line.get_ydata(), [0.0, 0.0])
    assert zero_line.get_color() == "black"
    assert zero_line.get_linestyle() == ":"
    np.testing.assert_allclose(median_line.get_xdata(), [0.0, 0.5, 1.0])
    np.testing.assert_allclose(median_line.get_ydata(), [2.0, 4.0, 3.0])
    assert median_line.get_color() == "black"
    assert median_line.get_linestyle() == "-"

    assert len(ax.collections) == 1
    envelope = ax.collections[0]
    np.testing.assert_allclose(
        envelope.get_facecolor()[0, :3],
        np.full(3, float(BLOOD_VOLUME_RATE_ENVELOPE_GRAY)),
    )
    vertices = envelope.get_paths()[0].vertices
    expected_sd = [np.std([1.0, 2.0, 3.0]), 0.0, np.std([1.0, 5.0])]
    for time, center, sd in zip([0.0, 0.5, 1.0], [2.0, 4.0, 3.0], expected_sd):
        y_values = vertices[np.isclose(vertices[:, 0], time), 1]
        assert np.isclose(np.min(y_values), center - sd)
        assert np.isclose(np.max(y_values), center + sd)

    assert all(spine.get_visible() for spine in ax.spines.values())
    assert all(np.isclose(spine.get_linewidth(), 0.4) for spine in ax.spines.values())
    width, height = fig.get_size_inches()
    assert np.isclose(width / height, LUMEN_DIAMETER_FIGURE_ASPECT_RATIO)


def test_blood_volume_rate_signals_export_png_and_eps_for_both_vessels() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        output = OutputManager.from_holo(
            Path(temp_dir) / "sample.holo",
            output_root=Path(temp_dir),
        )
        total = DatasetValue(
            np.asarray([[1.0, 2.0], [3.0, 5.0]], dtype=np.float32)
        )

        paths = export_blood_volume_rate_signals(output, total, total)

        assert {path.relative_to(output.layout.ef_dir).as_posix() for path in paths} == {
            "png/blood_volume_rate/artery_blood_volume_rate.png",
            "eps/blood_volume_rate/artery_blood_volume_rate.eps",
            "png/blood_volume_rate/vein_blood_volume_rate.png",
            "eps/blood_volume_rate/vein_blood_volume_rate.eps",
        }
        for path in paths:
            assert path.is_file()
            if path.suffix == ".png":
                with Image.open(path) as image:
                    assert image.format == "PNG"
                    assert np.isclose(
                        image.width / image.height,
                        LUMEN_DIAMETER_FIGURE_ASPECT_RATIO,
                        rtol=0.01,
                    )
                    image.verify()
            else:
                assert path.read_bytes().startswith(b"%!PS-Adobe")

