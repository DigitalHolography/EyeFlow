"""Tests for artery and vein velocity signal figures."""

from __future__ import annotations

import numpy as np

from pipelines.velocity_analysis.artifacts.velocity_signals import (
    VELOCITY_ENVELOPE_GRAY,
    VELOCITY_FIGURE_ASPECT_RATIO,
    _velocity_figure,
)


def _segment_values(spatial_medians: np.ndarray) -> np.ndarray:
    values = np.empty((*spatial_medians.shape, 2, 2), dtype=np.float64)
    for index in np.ndindex(spatial_medians.shape):
        center = spatial_medians[index]
        values[index] = [[center - 3.0, center - 1.0], [center + 1.0, center + 100.0]]
    return values


def test_velocity_figure_reduces_space_then_plots_beat_median_and_sd() -> None:
    spatial_medians = np.asarray(
        [
            [1.0, 2.0, 9.0],
            [4.0, 4.0, 4.0],
            [-5.0, 1.0, 7.0],
        ]
    )

    fig = _velocity_figure(_segment_values(spatial_medians))
    ax = fig.axes[0]

    assert ax.get_title() == ""
    assert ax.get_xlabel() == r"Cardiac Phase $t/T$"
    assert ax.get_ylabel() == r"$v(t)$ (mm/s)"
    assert ax.get_xlim() == (0.0, 1.0)
    assert len(ax.lines) == 2
    zero_line, median_line = ax.lines
    np.testing.assert_allclose(zero_line.get_ydata(), [0.0, 0.0])
    assert zero_line.get_color() == "black"
    assert zero_line.get_linestyle() == ":"

    expected_signal = np.median(spatial_medians, axis=1)
    expected_sd = np.std(spatial_medians, axis=1)
    np.testing.assert_allclose(median_line.get_xdata(), [0.0, 0.5, 1.0])
    np.testing.assert_allclose(median_line.get_ydata(), expected_signal)
    assert median_line.get_color() == "black"
    assert median_line.get_linestyle() == "-"

    assert len(ax.collections) == 1
    envelope = ax.collections[0]
    np.testing.assert_allclose(
        envelope.get_facecolor()[0, :3],
        np.full(3, float(VELOCITY_ENVELOPE_GRAY)),
    )
    vertices = envelope.get_paths()[0].vertices
    for time, center, sd in zip([0.0, 0.5, 1.0], expected_signal, expected_sd):
        y_values = vertices[np.isclose(vertices[:, 0], time), 1]
        assert np.isclose(np.min(y_values), center - sd)
        assert np.isclose(np.max(y_values), center + sd)

    assert all(spine.get_visible() for spine in ax.spines.values())
    assert all(np.isclose(spine.get_linewidth(), 0.4) for spine in ax.spines.values())
    width, height = fig.get_size_inches()
    assert np.isclose(width / height, VELOCITY_FIGURE_ASPECT_RATIO)
