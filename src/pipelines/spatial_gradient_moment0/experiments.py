"""Reproducible exports, diagnostics, metrics, and parameter sweeps."""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.signal import peak_prominences, peak_widths

from .config import (
    ConfigurationBundle,
    PreprocessingConfig,
    config_with_parameter,
)
from .preprocessing import (
    PipelineResult,
    processing_metadata,
    run_preprocessing_pipeline,
)


@dataclass
class ExperimentResult:
    experiment_id: str
    output_directory: Path
    preprocessing: PipelineResult
    profiles: np.ndarray
    metrics: dict[str, float]


@dataclass
class LumenSizeSweepResult:
    experiment_id: str
    output_directory: Path
    artery_time_branch: np.ndarray | None
    vein_time_branch: np.ndarray | None


def iter_sweep_configs(bundle: ConfigurationBundle) -> Iterable[PreprocessingConfig]:
    """Yield the Cartesian product requested by the central sweep section."""

    base_configs = (
        [bundle.resolve_preset(name) for name in bundle.sweep.presets]
        if bundle.sweep.presets
        else [bundle.preprocessing]
    )
    parameter_paths = list(bundle.sweep.parameters)
    candidate_lists = [bundle.sweep.parameters[path] for path in parameter_paths]
    combinations = itertools.product(*candidate_lists) if candidate_lists else [()]
    combinations = list(combinations)
    for base in base_configs:
        for values in combinations:
            config = base
            suffixes = []
            for path, value in zip(parameter_paths, values):
                config = config_with_parameter(config, path, value)
                suffixes.append(f"{path.split('.')[-1]}-{_slug_value(value)}")
            if suffixes:
                payload = config.to_dict()
                payload["name"] = f"{base.name}__{'__'.join(suffixes)}"
                config = PreprocessingConfig.from_dict(payload)
            yield config


def run_lumen_size_sweep(
    moment0ff,
    *,
    artery_geometry,
    vein_geometry,
    cycle_boundary_indexes,
    index_base: int,
    bundle: ConfigurationBundle,
    output_root: str | Path,
    source_dataset: str = "/moment0ff",
) -> list[LumenSizeSweepResult]:
    """Evaluate each config from masked tbkr lumen size across all branches."""

    if not bundle.sweep.enabled:
        return []
    from .runner import (
        SPATIAL_GRADIENT_METRICS_ROOT,
        _gradient_segments_from_velocity_geometry,
        _pack_vessel_spatial_gradient_profiles,
    )

    root = Path(output_root)
    root.mkdir(parents=True, exist_ok=True)
    results: list[LumenSizeSweepResult] = []
    for config in iter_sweep_configs(bundle):
        experiment_id = _experiment_id(config)
        output_directory = _unique_directory(root, experiment_id)
        output_directory.mkdir(parents=False)
        experiment_id = output_directory.name
        config_payload = config.to_dict()
        (output_directory / "configuration.json").write_text(
            json.dumps(config_payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )

        vessel_series: list[tuple[str, np.ndarray, np.ndarray]] = []
        row: dict[str, Any] = {
            "experiment_id": experiment_id,
            "pipeline_order": " -> ".join(config.pipeline),
            "gradient_method": config.gradient.method,
            "interpolation_order": processing_metadata(config)["interpolation_order"],
            **_flatten_dict(config_payload, prefix="config"),
            "output_directory": str(output_directory.resolve()),
        }
        for vessel_name, geometry in (
            ("Artery", artery_geometry),
            ("Vein", vein_geometry),
        ):
            if geometry is None:
                _add_missing_vessel_summary(row, vessel_name)
                continue
            processed_segments = _gradient_segments_from_velocity_geometry(
                moment0ff,
                geometry,
                source_dataset=source_dataset,
                preprocessing_config=config,
            )
            packed = _pack_vessel_spatial_gradient_profiles(
                processed_segments,
                vessel_name,
                cycle_boundary_indexes,
                index_base=index_base,
            )
            lumen_path = (
                f"{SPATIAL_GRADIENT_METRICS_ROOT}/{vessel_name}/Transverse/"
                "Masked/tbkr/lumen/size"
            )
            lumen_dataset = packed[lumen_path]
            time_branch = average_masked_lumen_size_tbkr(lumen_dataset)
            branch_ids = _resolved_branch_ids(processed_segments, time_branch.shape[1])
            vessel_series.append((vessel_name, time_branch, branch_ids))
            np.savez_compressed(
                output_directory
                / f"{vessel_name.lower()}_masked_lumen_size_time_branch.npz",
                lumen_size=time_branch,
                normalized_beat_time=_normalized_beat_time(len(time_branch)),
                branch_ids=branch_ids,
            )
            _add_vessel_summary(row, vessel_name, time_branch)
            del packed, processed_segments

        _write_lumen_size_branch_plot(
            output_directory / "masked_lumen_size_all_branches.png",
            vessel_series,
            title=config.name,
        )
        _append_results_row(root / "lumen_size_results.csv", row)
        results.append(
            LumenSizeSweepResult(
                experiment_id=experiment_id,
                output_directory=output_directory,
                artery_time_branch=_series_for_vessel(vessel_series, "Artery"),
                vein_time_branch=_series_for_vessel(vessel_series, "Vein"),
            )
        )
    return results


def average_masked_lumen_size_tbkr(lumen_dataset) -> np.ndarray:
    """Average masked lumen size over beat and radius, retaining time/branch."""

    from calculations.math import nanmean

    values = np.asarray(lumen_dataset.data, dtype=np.float32)
    dimensions = list((lumen_dataset.attrs or {}).get("dimDesc", ()))
    if values.ndim != 4 or dimensions != ["time", "beat", "branch", "radius"]:
        raise ValueError(
            "Masked lumen size must have dimensions (time, beat, branch, radius)."
        )
    return np.asarray(nanmean(values, axis=(1, 3)), dtype=np.float32)


def run_experiment(
    cropped_stack,
    *,
    rotation_angle: float,
    config: PreprocessingConfig,
    output_root: str | Path,
    mask: np.ndarray | None = None,
) -> ExperimentResult:
    """Run, evaluate with the existing profile/peaks, and export one experiment."""

    original = np.asarray(cropped_stack, dtype=np.float32)
    preprocessing = run_preprocessing_pipeline(
        original,
        rotation_angle=rotation_angle,
        config=config,
        capture_intermediates=config.output.save_intermediates,
    )
    profiles = transverse_profiles(preprocessing.data, mask=mask)
    metrics, peak_positions = profile_quality_metrics(profiles)
    root = Path(output_root)
    root.mkdir(parents=True, exist_ok=True)
    experiment_id = _experiment_id(config)
    output_directory = _unique_directory(root, experiment_id)
    output_directory.mkdir(parents=False)
    experiment_id = output_directory.name

    config_payload = config.to_dict()
    (output_directory / "configuration.json").write_text(
        json.dumps(config_payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    np.save(output_directory / "transverse_profiles.npy", profiles, allow_pickle=False)
    if config.output.save_final_float_tiff:
        _write_float_tiff(output_directory / "final_processed.tif", preprocessing.data)
    if config.output.save_intermediates:
        for index, (stage_name, stack) in enumerate(preprocessing.intermediates, start=1):
            _write_float_tiff(
                output_directory / f"{index:02d}_{_slug(stage_name)}.tif",
                stack,
            )
    if config.output.save_previews:
        for frame_index in _resolved_frame_indexes(
            config.output.diagnostic_frames, len(preprocessing.data)
        ):
            _write_diagnostic(
                output_directory / f"diagnostic_frame_{frame_index:05d}.png",
                original,
                preprocessing.data,
                profiles,
                peak_positions,
                frame_index,
                config.name,
            )
    row = {
        "experiment_id": experiment_id,
        "pipeline_order": " -> ".join(config.pipeline),
        "gradient_method": config.gradient.method,
        "interpolation_order": processing_metadata(config)["interpolation_order"],
        **_flatten_dict(config_payload, prefix="config"),
        **metrics,
        "output_directory": str(output_directory.resolve()),
    }
    _append_results_row(root / "results.csv", row)
    return ExperimentResult(
        experiment_id=experiment_id,
        output_directory=output_directory,
        preprocessing=preprocessing,
        profiles=profiles,
        metrics=metrics,
    )


def transverse_profiles(stack: np.ndarray, *, mask: np.ndarray | None = None) -> np.ndarray:
    """Call the existing transverse-profile reducer for every processed frame."""

    from .runner import _mean_transverse

    values = np.asarray(stack, dtype=np.float32)
    if mask is not None:
        mask_values = np.asarray(mask, dtype=bool)
        if mask_values.shape != values.shape[1:]:
            raise ValueError(
                f"Profile mask shape {mask_values.shape} does not match frames {values.shape[1:]}."
            )
        values = np.where(mask_values[None, :, :], values, np.nan)
    return np.stack([_mean_transverse(frame) for frame in values]).astype(
        np.float32, copy=False
    )


def profile_quality_metrics(
    profiles: np.ndarray,
    *,
    minimum_gap: int | None = None,
) -> tuple[dict[str, float], np.ndarray]:
    """Diagnostic metrics around the unchanged two-highest-peaks detector."""

    from .runner import (
        SPATIAL_GRADIENT_PEAK_MIN_GAP_SAMPLES,
        _fractional_peak_index,
        _two_highest_separated_indexes,
    )

    gap = (
        SPATIAL_GRADIENT_PEAK_MIN_GAP_SAMPLES
        if minimum_gap is None
        else int(minimum_gap)
    )
    curves = np.asarray(profiles, dtype=np.float32)
    if curves.ndim != 2:
        raise ValueError("Diagnostic profiles must have shape (time, x).")
    peak_positions = np.full((len(curves), 2), np.nan, dtype=np.float32)
    peak_values = np.full((len(curves), 2), np.nan, dtype=np.float32)
    prominences: list[float] = []
    widths: list[float] = []
    background_ratios: list[float] = []
    for time, curve in enumerate(curves):
        indexes = _two_highest_separated_indexes(curve, gap)
        for side, index in enumerate(indexes):
            peak_positions[time, side] = _fractional_peak_index(curve, index)
            peak_values[time, side] = curve[index]
        if len(indexes) != 2:
            continue
        finite = np.isfinite(curve)
        if np.count_nonzero(finite) < 3:
            continue
        finite_curve = np.where(finite, curve, np.nanmin(curve[finite]))
        for index in indexes:
            if 0 < index < len(curve) - 1 and (
                finite_curve[index] >= finite_curve[index - 1]
                and finite_curve[index] >= finite_curve[index + 1]
            ):
                prominence = float(peak_prominences(finite_curve, [index])[0][0])
                width = float(peak_widths(finite_curve, [index], rel_height=0.5)[0][0])
                if math.isfinite(prominence):
                    prominences.append(prominence)
                if math.isfinite(width):
                    widths.append(width)
        background = curve[finite].astype(np.float64)
        for index in indexes:
            background = background[background != float(curve[index])]
        if background.size:
            median = float(np.median(background))
            mad = float(np.median(np.abs(background - median)))
            noise = 1.4826 * mad
            if noise > 0.0:
                background_ratios.extend(
                    (float(curve[index]) - median) / noise for index in indexes
                )
    detected = np.all(np.isfinite(peak_positions), axis=1)
    separations = peak_positions[:, 1] - peak_positions[:, 0]
    return (
        {
            "frames": float(len(curves)),
            "missing_or_ambiguous_peak_fraction": float(1.0 - np.mean(detected))
            if len(curves)
            else float("nan"),
            "peak_prominence_mean": _finite_mean(prominences),
            "peak_width_half_prominence_mean": _finite_mean(widths),
            "peak_to_background_robust_mean": _finite_mean(background_ratios),
            "left_peak_x_jitter_std": _finite_std(peak_positions[:, 0]),
            "right_peak_x_jitter_std": _finite_std(peak_positions[:, 1]),
            "peak_separation_mean": _finite_mean(separations),
            "peak_separation_std": _finite_std(separations),
            "left_peak_value_mean": _finite_mean(peak_values[:, 0]),
            "right_peak_value_mean": _finite_mean(peak_values[:, 1]),
        },
        peak_positions,
    )


def _write_diagnostic(
    path: Path,
    original: np.ndarray,
    processed: np.ndarray,
    profiles: np.ndarray,
    peak_positions: np.ndarray,
    frame_index: int,
    title: str,
) -> None:
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    original_index = min(frame_index, len(original) - 1)
    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    axes[0].imshow(original[original_index], cmap="gray")
    axes[0].set_title("Existing geometry crop (float)")
    axes[1].imshow(processed[frame_index], cmap="gray")
    axes[1].set_title("Processed (display-normalized only)")
    x = np.arange(profiles.shape[1])
    axes[2].plot(x, profiles[frame_index], color="black", linewidth=1.2)
    for side, color in ((0, "tab:blue"), (1, "tab:red")):
        position = peak_positions[frame_index, side]
        if np.isfinite(position):
            axes[2].axvline(float(position), color=color, linestyle="--", linewidth=1)
    axes[2].set_title("Existing transverse profile + peaks")
    axes[2].set_xlabel("X sample")
    axes[2].set_ylabel("Processed intensity / gradient")
    for axis in axes[:2]:
        axis.axis("off")
    fig.suptitle(f"{title} — frame {frame_index}")
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _write_lumen_size_branch_plot(
    path: Path,
    vessel_series: list[tuple[str, np.ndarray, np.ndarray]],
    *,
    title: str,
) -> None:
    """Plot every artery and vein branch on one normalized-beat-time graph."""

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(figsize=(12, 6))
    line_count = sum(values.shape[1] for _, values, _ in vessel_series)
    color_map = plt.get_cmap("turbo", max(line_count, 1))
    color_index = 0
    plotted = 0
    for vessel_name, values, branch_ids in vessel_series:
        time = _normalized_beat_time(len(values))
        line_style = "-" if vessel_name == "Artery" else "--"
        for branch_index in range(values.shape[1]):
            curve = values[:, branch_index]
            color = color_map(color_index)
            color_index += 1
            if not np.any(np.isfinite(curve)):
                continue
            branch_id = branch_ids[branch_index]
            axis.plot(
                time,
                curve,
                color=color,
                linestyle=line_style,
                linewidth=1.2,
                label=f"{vessel_name} branch {branch_id}",
            )
            plotted += 1
    if not plotted:
        axis.text(
            0.5,
            0.5,
            "No finite masked lumen-size branch curves",
            ha="center",
            va="center",
            transform=axis.transAxes,
        )
    axis.set_title(f"{title}: masked lumen size, mean over beats and radii")
    axis.set_xlabel("Normalized time within beat")
    axis.set_ylabel("Lumen size (pixels)")
    axis.grid(True, alpha=0.25)
    axis.set_xlim(0.0, 1.0)
    if plotted:
        axis.legend(
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            fontsize=7,
            frameon=False,
        )
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _normalized_beat_time(time_count: int) -> np.ndarray:
    if time_count < 1:
        return np.empty(0, dtype=np.float32)
    return np.arange(time_count, dtype=np.float32) / np.float32(time_count)


def _resolved_branch_ids(segments, branch_count: int) -> np.ndarray:
    branch_ids = np.asarray(
        getattr(segments, "branch_ids", np.arange(branch_count)),
    ).reshape(-1)
    if branch_ids.size != branch_count:
        return np.arange(branch_count, dtype=np.int32)
    return branch_ids


def _add_vessel_summary(
    row: dict[str, Any],
    vessel_name: str,
    time_branch: np.ndarray,
) -> None:
    prefix = vessel_name.lower()
    finite = np.isfinite(time_branch)
    row[f"{prefix}_branch_count"] = int(time_branch.shape[1])
    row[f"{prefix}_time_branch_valid_fraction"] = (
        float(np.mean(finite)) if finite.size else float("nan")
    )
    row[f"{prefix}_lumen_size_mean_pixels"] = _finite_mean(time_branch)
    row[f"{prefix}_lumen_size_std_pixels"] = _finite_std(time_branch)


def _add_missing_vessel_summary(row: dict[str, Any], vessel_name: str) -> None:
    prefix = vessel_name.lower()
    row[f"{prefix}_branch_count"] = 0
    row[f"{prefix}_time_branch_valid_fraction"] = float("nan")
    row[f"{prefix}_lumen_size_mean_pixels"] = float("nan")
    row[f"{prefix}_lumen_size_std_pixels"] = float("nan")


def _series_for_vessel(
    vessel_series: list[tuple[str, np.ndarray, np.ndarray]],
    vessel_name: str,
) -> np.ndarray | None:
    return next(
        (values for name, values, _ in vessel_series if name == vessel_name),
        None,
    )


def _write_float_tiff(path: Path, stack: np.ndarray) -> None:
    try:
        import tifffile
    except ModuleNotFoundError as exc:
        raise ImportError("Float TIFF experiment exports require tifffile.") from exc
    tifffile.imwrite(
        path,
        np.asarray(stack, dtype=np.float32),
        photometric="minisblack",
        metadata={"axes": "TYX"},
    )


def _append_results_row(path: Path, row: dict[str, Any]) -> None:
    normalized = {
        key: json.dumps(value, sort_keys=True) if isinstance(value, (list, dict)) else value
        for key, value in row.items()
    }
    exists = path.is_file() and path.stat().st_size > 0
    if exists:
        with path.open("r", newline="", encoding="utf-8") as source:
            header = next(csv.reader(source))
        missing = [key for key in normalized if key not in header]
        if missing:
            raise ValueError(
                f"Existing results table {path} lacks columns: {', '.join(missing)}."
            )
        normalized = {key: normalized.get(key, "") for key in header}
    else:
        header = list(normalized)
    with path.open("a", newline="", encoding="utf-8") as destination:
        writer = csv.DictWriter(destination, fieldnames=header)
        if not exists:
            writer.writeheader()
        writer.writerow(normalized)


def _flatten_dict(values: dict[str, Any], *, prefix: str) -> dict[str, Any]:
    flattened: dict[str, Any] = {}
    for key, value in values.items():
        path = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict):
            flattened.update(_flatten_dict(value, prefix=path))
        else:
            flattened[path] = value
    return flattened


def _experiment_id(config: PreprocessingConfig) -> str:
    encoded = json.dumps(config.to_dict(), sort_keys=True, separators=(",", ":")).encode()
    digest = hashlib.sha256(encoded).hexdigest()[:10]
    return f"{_slug(config.name)[:80]}__{digest}"


def _unique_directory(root: Path, name: str) -> Path:
    candidate = root / name
    sequence = 2
    while candidate.exists():
        candidate = root / f"{name}__run{sequence:03d}"
        sequence += 1
    return candidate


def _slug(value: str) -> str:
    cleaned = re.sub(r"[^a-zA-Z0-9._-]+", "-", str(value)).strip("-_")
    return cleaned or "experiment"


def _slug_value(value: Any) -> str:
    if isinstance(value, list):
        return _slug("-".join(str(item) for item in value))
    return _slug(str(value))


def _resolved_frame_indexes(indexes: list[int], frame_count: int) -> list[int]:
    resolved = []
    for index in indexes:
        actual = index if index >= 0 else frame_count + index
        if 0 <= actual < frame_count and actual not in resolved:
            resolved.append(actual)
    return resolved


def _finite_mean(values) -> float:
    array = np.asarray(values, dtype=np.float64)
    finite = array[np.isfinite(array)]
    return float(np.mean(finite)) if finite.size else float("nan")


def _finite_std(values) -> float:
    array = np.asarray(values, dtype=np.float64)
    finite = array[np.isfinite(array)]
    return float(np.std(finite, ddof=0)) if finite.size else float("nan")


__all__ = [
    "ExperimentResult",
    "LumenSizeSweepResult",
    "average_masked_lumen_size_tbkr",
    "iter_sweep_configs",
    "profile_quality_metrics",
    "run_experiment",
    "run_lumen_size_sweep",
    "transverse_profiles",
]
