"""Serialize blood-volume-rate calculations to the EyeFlow schema."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from calculations.blood_volume_rate import (
    TOTAL_MASKED_EDGES_WINDOW_SIZE,
    TOTAL_MASKED_EDGES_WINDOW_STRIDE,
    circular_lumen_flow,
    mask_derived_lumen_geometry,
    total_masked_edges_flow,
)
from calculations.math import nanmean_float32, nanmedian
from input_output.profile_datasets import _profile_dataset, _profile_h5_options
from input_output.schema import EyeFlowOutputPaths
from input_output.writers.eps import EpsArtifactWriter, write_eps_file
from input_output.writers.png import PngArtifactWriter
from pipeline_engine.base import DatasetValue

LUMEN_DIAMETER_BIN_WIDTH_MICRONS = 5.0
LUMEN_DIAMETER_FIGURE_DPI = 320
LUMEN_DIAMETER_FIGURE_WIDTH_INCHES = 8.0 / 2.54
LUMEN_DIAMETER_FIGURE_ASPECT_RATIO = 1.618
LUMEN_DIAMETER_FONT_FAMILY = "Times New Roman"
LUMEN_DIAMETER_FONT_SIZE = 8
LUMEN_DIAMETER_LEGEND_FONT_SIZE = 8
LUMEN_DIAMETER_HISTOGRAM_GRAY = "0.8"
LUMEN_DIAMETER_HISTOGRAM_EDGE_GRAY = "0.4"
LUMEN_DIAMETER_BOX_LINE_WIDTH = 0.4
LUMEN_DIAMETER_GAUSSIAN_LINE_WIDTH = 0.8
BLOOD_VOLUME_RATE_ENVELOPE_GRAY = "0.80"
BLOOD_VOLUME_RATE_LINE_WIDTH = 0.8  


def pack_gradient_edge_outputs(
    artery_velocity_segments,
    vein_velocity_segments,
    gradient_products,
    cycle_boundary_indexes,
    *,
    index_base: int,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, DatasetValue]:
    """Calculate dynamic- and static-edge flow for both vessel classes."""

    schema = _resolve_output_paths(output_paths)
    outputs: dict[str, DatasetValue] = {}
    vessels = (
        (
            "Artery",
            artery_velocity_segments,
            gradient_products.artery_segments,
            schema.blood_volume_rate.artery,
        ),
        (
            "Vein",
            vein_velocity_segments,
            gradient_products.vein_segments,
            schema.blood_volume_rate.vein,
        ),
    )
    for vessel_name, velocity, gradient, paths in vessels:
        _validate_profile_segment_alignment(vessel_name, velocity, gradient)
        velocity_profile = velocity.profile
        profile = _profile_dataset(
            np.asarray(velocity_profile.transverse.masked, dtype=np.float32),
            cycle_boundary_indexes,
            index_base=index_base,
            spatial_axis="x",
            valid_segments=np.asarray(
                velocity_profile.topology.valid_segments, dtype=bool
            ),
        )
        metrics_root = (
            f"Processing/SpatialGradientMetrics/{vessel_name}/"
            "Transverse/Masked/tbkr"
        )
        left_path = f"{metrics_root}/left_edge_index"
        right_path = f"{metrics_root}/right_edge_index"
        left_value = gradient_products.outputs[left_path]
        right_value = gradient_products.outputs[right_path]
        pixel_size_mm = float(velocity_profile.sample_spacing_mm)
        outputs[paths.dynamic_edges] = _gradient_edge_dataset(
            profile,
            left_value,
            right_value,
            profile_pixel_size_mm=pixel_size_mm,
            left_edge_path=left_path,
            right_edge_path=right_path,
            static_edges=False,
        )
        outputs[paths.static_edges] = _gradient_edge_dataset(
            profile,
            left_value,
            right_value,
            profile_pixel_size_mm=pixel_size_mm,
            left_edge_path=left_path,
            right_edge_path=right_path,
            static_edges=True,
        )
    return outputs


def _gradient_edge_dataset(
    profile: DatasetValue,
    left_edge: DatasetValue,
    right_edge: DatasetValue,
    *,
    profile_pixel_size_mm: float,
    left_edge_path: str,
    right_edge_path: str,
    static_edges: bool,
) -> DatasetValue:
    left = np.asarray(left_edge.data, dtype=np.float32)
    right = np.asarray(right_edge.data, dtype=np.float32)
    if static_edges:
        left = np.broadcast_to(nanmean_float32(left, axis=(0, 1)), left.shape)
        right = np.broadcast_to(nanmean_float32(right, axis=(0, 1)), right.shape)
    if not np.isfinite(profile_pixel_size_mm) or profile_pixel_size_mm <= 0:
        raise ValueError("profile_pixel_size_mm must be finite and positive.")
    velocity = nanmean_float32(profile.data, axis=0)
    diameter_mm = np.where(
        right > left,
        (right - left) * np.float32(profile_pixel_size_mm),
        np.float32(np.nan),
    )
    rate = circular_lumen_flow(velocity, diameter_mm)
    return DatasetValue(
        rate,
        {
            "unit": "mm^3/s",
            "dimDesc": ["time", "beat", "branch", "radius"],
            "definition": (
                "mean masked transverse velocity multiplied by the circular "
                "lumen area implied by the spatial-gradient edges"
            ),
            "source_velocity": "waveform_velocity.masked_transverse_profile",
            "source_left_edge_index": f"/{left_edge_path.lstrip('/')}",
            "source_right_edge_index": f"/{right_edge_path.lstrip('/')}",
            "diameter_model": "gradient_edge_separation_times_profile_pixel_size",
            "cross_section_model": "circular_pi_diameter_squared_over_4",
            "velocity_reduction": "mean_over_finite_transverse_profile_samples",
            "profile_pixel_size_mm": np.float32(profile_pixel_size_mm),
            "edge_temporal_reduction": (
                "mean_over_time_and_beats" if static_edges else "none"
            ),
        },
        h5_options=_profile_h5_options(rate.shape),
    )


def pack_mask_derived_outputs(
    prepared_topologies,
    velocity_per_beat_outputs: dict[str, object],
    *,
    pixel_size_mm: float,
    output_paths: EyeFlowOutputPaths | str | None = None,
) -> dict[str, DatasetValue]:
    """Calculate masked-edge and total masked-edge flow for both vessels."""

    schema = _resolve_output_paths(output_paths)
    diameters, _, radial_widths = mask_derived_lumen_geometry(
        (prepared_topologies["artery"], prepared_topologies["vein"]),
        pixel_size_mm=pixel_size_mm,
    )
    vessel_sources = (
        (
            "Artery",
            schema.artery_per_beat_safe.velocity_signal,
            schema.blood_volume_rate.artery,
            diameters[0],
        ),
        (
            "Vein",
            schema.vein_per_beat_safe.velocity_signal,
            schema.blood_volume_rate.vein,
            diameters[1],
        ),
    )
    outputs: dict[str, DatasetValue] = {}
    for vessel_name, velocity_path, paths, diameter_mm in vessel_sources:
        if velocity_path is None or velocity_path not in velocity_per_beat_outputs:
            raise KeyError(
                f"Required safe per-beat velocity is unavailable for {vessel_name}."
            )
        velocity = _metric_data(velocity_per_beat_outputs[velocity_path])
        rate = circular_lumen_flow(velocity, diameter_mm)
        masked = DatasetValue(
            rate,
            {
                "unit": "mm^3/s",
                "dimDesc": ["time", "beat", "branch", "radius"],
                "definition": (
                    "safe per-beat segment velocity multiplied by an equivalent "
                    "circular lumen area derived from native vessel-mask pixels"
                ),
                "source_velocity": "waveform_velocity.safe_per_beat_segment_velocity",
                "diameter_model": "masked_pixel_count_over_radial_width",
                "diameter_model_assumption": "locally_radial_vessel",
                "cross_section_model": "circular_pi_diameter_squared_over_4",
                "annulus_geometry": "native_pixel_center_section_mask",
                "annulus_edge_handling": "outer_radius_clipped_to_configured_limit",
                "native_pixel_size_mm": np.float32(pixel_size_mm),
                "radial_width_pixels": radial_widths,
            },
            h5_options=_profile_h5_options(rate.shape),
        )
        outputs[paths.masked_edges] = masked
        total = total_masked_edges_flow(rate)
        outputs[paths.total_masked_edges] = DatasetValue(
            total,
            {
                "unit": "mm^3/s",
                "dimDesc": ["time", "beat"],
                "definition": (
                    "median over radius of the sum over branches of masked-edge "
                    "blood-volume rate after a circular sliding average over tau"
                ),
                "source": f"/{paths.masked_edges.lstrip('/')}",
                "aggregation": "median_over_radius_of_sum_over_branches",
                "temporal_filter": "circular_sliding_average_over_tau",
                "temporal_window_size": np.int32(TOTAL_MASKED_EDGES_WINDOW_SIZE),
                "temporal_window_stride": np.int32(TOTAL_MASKED_EDGES_WINDOW_STRIDE),
                "temporal_boundary_mode": "circular",
                "temporal_window_alignment": "centered",
                "temporal_nan_policy": "propagate",
                "branch_reduction": "sum_over_finite_values",
                "radius_reduction": "median_over_finite_values",
            },
            h5_options=_profile_h5_options(total.shape),
        )
    return outputs


def export_lumen_diameter_distributions(
    output,
    artery_lumen_diameter_pixels,
    vein_lumen_diameter_pixels,
    *,
    pixel_pitch_m: float,
) -> list[Path]:
    """Export artery and vein lumen-diameter histograms as PNG and EPS."""

    pitch = np.asarray(pixel_pitch_m, dtype=np.float64)
    if pitch.size != 1:
        raise ValueError("pixel_pitch_m must be a scalar.")
    pitch_value = float(pitch.reshape(()))
    if not np.isfinite(pitch_value) or pitch_value <= 0.0:
        raise ValueError("pixel_pitch_m must be finite and positive.")

    paths: list[Path] = []
    for vessel, diameter_pixels in (
        ("artery", artery_lumen_diameter_pixels),
        ("vein", vein_lumen_diameter_pixels),
    ):
        diameter_microns = (
            np.asarray(diameter_pixels, dtype=np.float64).reshape(-1)
            * pitch_value
            * 1e6
        )
        diameter_microns = diameter_microns[np.isfinite(diameter_microns)]
        fig = _lumen_diameter_distribution_figure(diameter_microns)
        filename = f"lumen_diameter/{vessel}_lumen_diameter_distribution"
        png_path = output.path_for(_png_output_type(), f"{filename}.png")
        eps_path = output.path_for(_eps_output_type(), f"{filename}.eps")
        png_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(png_path, dpi=LUMEN_DIAMETER_FIGURE_DPI)
        write_eps_file(eps_path, fig, dpi=LUMEN_DIAMETER_FIGURE_DPI)
        _close_figure(fig)
        paths.extend((png_path, eps_path))
    return paths


def export_blood_volume_rate_signals(
    output,
    artery_total_masked_edges,
    vein_total_masked_edges,
) -> list[Path]:
    """Export median total masked-edge flow and its beatwise SD."""

    paths: list[Path] = []
    for vessel, values in (
        ("artery", artery_total_masked_edges),
        ("vein", vein_total_masked_edges),
    ):
        stem = f"blood_volume_rate/{vessel}"
        paths.append(
            PngArtifactWriter(output, stem).save_figure(
                _blood_volume_rate_figure(_metric_data(values)),
                "blood_volume_rate.png",
                dpi=LUMEN_DIAMETER_FIGURE_DPI,
                bbox_inches=None,
            )
        )
        paths.append(
            EpsArtifactWriter(output, stem).save_figure(
                _blood_volume_rate_figure(_metric_data(values)),
                "blood_volume_rate.eps",
                dpi=LUMEN_DIAMETER_FIGURE_DPI,
            )
        )
    return paths


def _blood_volume_rate_figure(total_masked_edges: np.ndarray):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.text import Text

    values = np.asarray(total_masked_edges, dtype=np.float64)
    if values.ndim != 2:
        raise ValueError(
            "total_masked_edges must have exactly two dimensions (time, beat)."
        )

    figure_height = (
        LUMEN_DIAMETER_FIGURE_WIDTH_INCHES
        / LUMEN_DIAMETER_FIGURE_ASPECT_RATIO
    )
    fig = Figure(
        figsize=(LUMEN_DIAMETER_FIGURE_WIDTH_INCHES, figure_height),
        constrained_layout=True,
    )
    FigureCanvasAgg(fig)
    ax = fig.subplots()

    median = nanmedian(values, axis=1, dtype=np.float64)
    finite = np.isfinite(values)
    finite_count = np.sum(finite, axis=1)
    beat_mean = np.divide(
        np.sum(np.where(finite, values, 0.0), axis=1),
        finite_count,
        out=np.full(values.shape[0], np.nan, dtype=np.float64),
        where=finite_count > 0,
    )
    squared_deviation = np.where(finite, (values - beat_mean[:, None]) ** 2, 0.0)
    variance = np.divide(
        np.sum(squared_deviation, axis=1),
        finite_count,
        out=np.full(values.shape[0], np.nan, dtype=np.float64),
        where=finite_count > 0,
    )
    standard_deviation = np.sqrt(variance)
    cardiac_cycle = np.linspace(0.0, 1.0, values.shape[0])
    ax.fill_between(
        cardiac_cycle,
        median - standard_deviation,
        median + standard_deviation,
        color=BLOOD_VOLUME_RATE_ENVELOPE_GRAY,
        linewidth=0.0,
    )
    ax.axhline(
        0.0,
        color="black",
        linestyle=":",
        linewidth=LUMEN_DIAMETER_BOX_LINE_WIDTH,
    )
    ax.plot(
        cardiac_cycle,
        median,
        color="black",
        linestyle="-",
        linewidth=BLOOD_VOLUME_RATE_LINE_WIDTH,
    )

    ax.set_xlim(0.0, 1.0)
    ax.set_xlabel(r"Cardiac Phase $t/T$", fontsize=LUMEN_DIAMETER_FONT_SIZE)
    ax.set_ylabel(r"$Q(t)$ (mm³/s)", fontsize=LUMEN_DIAMETER_FONT_SIZE)
    ax.tick_params(axis="both", labelsize=LUMEN_DIAMETER_FONT_SIZE)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(LUMEN_DIAMETER_BOX_LINE_WIDTH)
    for text in fig.findobj(match=Text):
        text.set_fontfamily(LUMEN_DIAMETER_FONT_FAMILY)
    return fig


def _lumen_diameter_distribution_figure(diameter_microns: np.ndarray):
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.lines import Line2D
    from matplotlib.text import Text

    figure_height = (
        LUMEN_DIAMETER_FIGURE_WIDTH_INCHES
        / LUMEN_DIAMETER_FIGURE_ASPECT_RATIO
    )
    fig = Figure(
        figsize=(LUMEN_DIAMETER_FIGURE_WIDTH_INCHES, figure_height),
        constrained_layout=True,
    )
    FigureCanvasAgg(fig)
    ax = fig.subplots()
    values = np.asarray(diameter_microns, dtype=np.float64).reshape(-1)
    values = values[np.isfinite(values)]

    median = float(np.median(values)) if values.size else np.nan
    standard_deviation = float(np.std(values)) if values.size else np.nan
    edges = _lumen_diameter_histogram_edges(values)
    ax.hist(
        values,
        bins=edges,
        color=LUMEN_DIAMETER_HISTOGRAM_GRAY,
        edgecolor=LUMEN_DIAMETER_HISTOGRAM_EDGE_GRAY,
        linewidth=0.25,
    )

    if values.size and standard_deviation > 0.0:
        x = np.linspace(edges[0], edges[-1], 512)
        probability_density = np.exp(
            -0.5 * ((x - median) / standard_deviation) ** 2
        ) / (standard_deviation * np.sqrt(2.0 * np.pi))
        ax.plot(
            x,
            probability_density
            * values.size
            * LUMEN_DIAMETER_BIN_WIDTH_MICRONS,
            color="black",
            linestyle="--",
            linewidth=LUMEN_DIAMETER_GAUSSIAN_LINE_WIDTH,
        )

    statistic_handles = [
        Line2D([], [], color="none", label=_statistic_label("Median", median)),
        Line2D(
            [],
            [],
            color="none",
            label=_statistic_label("SD", standard_deviation),
        ),
    ]
    ax.legend(
        handles=statistic_handles,
        loc="upper right",
        fontsize=LUMEN_DIAMETER_LEGEND_FONT_SIZE,
        handlelength=2.0,
        frameon=False,
    )
    ax.set_xlabel("Lumen diameter (µm)", fontsize=LUMEN_DIAMETER_FONT_SIZE)
    ax.set_ylabel("Count", fontsize=LUMEN_DIAMETER_FONT_SIZE)
    ax.tick_params(axis="both", labelsize=LUMEN_DIAMETER_FONT_SIZE)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(LUMEN_DIAMETER_BOX_LINE_WIDTH)
    for text in fig.findobj(match=Text):
        text.set_fontfamily(LUMEN_DIAMETER_FONT_FAMILY)
    return fig


def _lumen_diameter_histogram_edges(values: np.ndarray) -> np.ndarray:
    width = LUMEN_DIAMETER_BIN_WIDTH_MICRONS
    if values.size == 0:
        return np.asarray([0.0, width])
    lower = np.floor(float(np.min(values)) / width) * width
    upper = np.ceil(float(np.max(values)) / width) * width
    if upper <= lower:
        upper = lower + width
    bin_count = max(1, round((upper - lower) / width))
    return lower + np.arange(bin_count + 1, dtype=np.float64) * width


def _statistic_label(name: str, value: float) -> str:
    if not np.isfinite(value):
        return f"{name}: n/a"
    return f"{name}: {value:.1f} µm"


def _close_figure(fig) -> None:
    import matplotlib.pyplot as plt

    plt.close(fig)


def _png_output_type():
    from input_output.output_manager import OutputType

    return OutputType.PNG


def _eps_output_type():
    from input_output.output_manager import OutputType

    return OutputType.EPS


def _validate_profile_segment_alignment(vessel_name, velocity, gradient) -> None:
    velocity_topology = velocity.profile.topology
    gradient_topology = gradient.topology
    if (
        velocity_topology is not gradient_topology
        and not velocity_topology.is_aligned_with(gradient_topology)
    ):
        raise RuntimeError(f"{vessel_name} profile segment topologies do not align.")


def _metric_data(value) -> np.ndarray:
    if isinstance(value, DatasetValue):
        value = value.data
    elif isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], dict):
        value = value[0]
    return np.asarray(value)


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)


__all__ = [
    "export_blood_volume_rate_signals",
    "export_lumen_diameter_distributions",
    "pack_gradient_edge_outputs",
    "pack_mask_derived_outputs",
]
