"""AVI mosaics for rotated per-segment velocity maps."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
import math
from pathlib import Path
from time import perf_counter, process_time
import warnings

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from input_output.writers.avi import AviArtifactWriter
from utils.logger import Logger


SEGMENT_VELOCITY_MAP_AVI_FOLDER = "segment_velocity_map"
SEGMENT_VELOCITY_MAP_FPS = 60.0
SEGMENT_VELOCITY_MAP_COLORMAP = "turbo"
SEGMENT_VELOCITY_MAP_MAX_SIDE = 1024
SEGMENT_VELOCITY_MAP_EXPORTED_BEATS = 1
SEGMENT_VELOCITY_MAP_FRAME_CHUNK_SIZE = 4


@dataclass(frozen=True)
class _MosaicTile:
    branch_index: int
    radius_index: int
    branch_id: int
    row: int
    column: int
    label: str


@dataclass(frozen=True)
class _MosaicRenderPlan:
    x_indexes: np.ndarray
    y_indexes: np.ndarray
    branch_indexes: np.ndarray
    radius_indexes: np.ndarray
    x_offset: int
    y_offset: int
    display_x_count: int
    display_y_count: int


@dataclass(frozen=True)
class _MosaicSource:
    vessel: str
    maps: np.ndarray
    tiles: tuple[_MosaicTile, ...]
    grid_side: int
    tile_side: int
    exported_beat_indexes: tuple[int, ...]
    velocity_minimum: float
    velocity_maximum: float
    render_plan: _MosaicRenderPlan

    @property
    def time_count(self) -> int:
        return int(self.maps.shape[2])

    @property
    def beat_count(self) -> int:
        return int(self.maps.shape[3])

    @property
    def video_side(self) -> int:
        mosaic_side = self.grid_side * self.tile_side
        return mosaic_side + mosaic_side % 2


@dataclass(frozen=True)
class _MosaicPreparationStats:
    candidate_count: int = 0
    candidate_discovery_seconds: float = 0.0
    range_scan_seconds: float = 0.0
    range_values_inspected: int = 0


def export_segment_velocity_map_avis(
    output,
    artery_segments,
    vein_segments,
    packed_segment_maps: dict[str, object],
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    fps: float = SEGMENT_VELOCITY_MAP_FPS,
) -> list[str]:
    """Export artery and vein per-beat segment mosaics with one color scale."""
    if not getattr(output, "available", True):
        return []
    export_wall_started = perf_counter()
    export_cpu_started = process_time()
    schema = _resolve_output_paths(output_paths)
    sources: list[_MosaicSource | None] = []
    for vessel, segments, path in (
        (
            "artery",
            artery_segments,
            schema.artery_segments.velocity_map_per_segment,
        ),
        (
            "vein",
            vein_segments,
            schema.vein_segments.velocity_map_per_segment,
        ),
    ):
        source, stats = _prepare_mosaic_source(
            vessel,
            segments,
            packed_segment_maps.get(path),
        )
        sources.append(source)
        if source is not None:
            Logger.log(
                f"Prepared {vessel} segment AVI source: "
                f"shape={source.maps.shape}, "
                f"input={source.maps.nbytes / (1024**3):.2f} GiB, "
                f"candidates={stats.candidate_count}, tiles={len(source.tiles)}, "
                f"candidate discovery={stats.candidate_discovery_seconds:.3f}s, "
                f"median-image range scan={stats.range_scan_seconds:.3f}s "
                f"({stats.range_values_inspected:,} values per reduction)."
            )
    available_sources = tuple(source for source in sources if source is not None)
    if not available_sources:
        return []

    range_started = perf_counter()
    velocity_range = _global_velocity_range(available_sources)
    Logger.log(
        "Resolved shared segment AVI velocity range in "
        f"{perf_counter() - range_started:.3f}s: "
        f"[{velocity_range[0]:.3f}, {velocity_range[1]:.3f}] mm/s."
    )
    color_lut = _turbo_lut()
    writer = AviArtifactWriter(output)

    def export_source(source: _MosaicSource) -> str:
        metadata = _video_metadata(source, velocity_range, fps)
        vessel_paths = (
            schema.artery_segments
            if source.vessel == "artery"
            else schema.vein_segments
        )
        metadata["velocity_map_dataset"] = _absolute_h5_path(
            vessel_paths.velocity_map_per_segment
        )
        metadata["segment_dataset"] = _absolute_h5_path(vessel_paths.segments)
        suffix = f"{source.vessel}_segment_velocity_map.avi"
        with writer.open_mjpeg(
            suffix,
            width=source.video_side,
            height=source.video_side,
            fps=fps,
            subfolder=SEGMENT_VELOCITY_MAP_AVI_FOLDER,
            metadata=metadata,
        ) as video:
            render_seconds = 0.0
            for beat_index in source.exported_beat_indexes:
                for time_start in range(
                    0,
                    source.time_count,
                    SEGMENT_VELOCITY_MAP_FRAME_CHUNK_SIZE,
                ):
                    time_stop = min(
                        time_start + SEGMENT_VELOCITY_MAP_FRAME_CHUNK_SIZE,
                        source.time_count,
                    )
                    render_started = perf_counter()
                    frames = _mosaic_frame_chunk(
                        source,
                        np.arange(time_start, time_stop, dtype=np.intp),
                        beat_index,
                        velocity_range,
                        color_lut,
                    )
                    render_seconds += perf_counter() - render_started
                    for frame in frames:
                        video.write_frame(frame)
        performance = video.performance_stats
        output_bytes = video.path.stat().st_size
        Logger.log(
            f"Completed {source.vessel} segment AVI: "
            f"frames={video.frame_count}, size={source.video_side}x{source.video_side}, "
            f"output={output_bytes / (1024**2):.1f} MiB, "
            f"render={render_seconds:.3f}s, "
            f"RGB conversion={performance['rgb_conversion_seconds']:.3f}s, "
            f"JPEG={performance['jpeg_encode_seconds']:.3f}s, "
            f"AVI packing={performance['chunk_pack_seconds']:.3f}s, "
            f"write={performance['file_write_seconds']:.3f}s, "
            f"finalize={performance['finalize_seconds']:.3f}s."
        )
        return str(video.path)

    if len(available_sources) == 1:
        paths = [export_source(available_sources[0])]
    else:
        with ThreadPoolExecutor(
            max_workers=min(2, len(available_sources)),
            thread_name_prefix="segment-avi",
        ) as executor:
            paths = list(executor.map(export_source, available_sources))
    output_bytes = sum(Path(path).stat().st_size for path in paths)
    Logger.log(
        "Segment AVI export totals: "
        f"wall={perf_counter() - export_wall_started:.3f}s, "
        f"process CPU={process_time() - export_cpu_started:.3f}s, "
        f"output={output_bytes / (1024**2):.1f} MiB."
    )
    return paths


def _mosaic_source(vessel: str, segments, packed_value) -> _MosaicSource | None:
    source, _ = _prepare_mosaic_source(vessel, segments, packed_value)
    return source


def _prepare_mosaic_source(
    vessel: str,
    segments,
    packed_value,
) -> tuple[_MosaicSource | None, _MosaicPreparationStats]:
    if packed_value is None:
        return None, _MosaicPreparationStats()
    maps = np.asarray(getattr(packed_value, "data", packed_value), dtype=np.float32)
    if maps.ndim != 6:
        raise ValueError(
            "VelocityMapPerSegment must have shape "
            "(x, y, time, beat, branch, radius)."
        )
    branch_ids = np.asarray(
        getattr(segments, "branch_ids", np.arange(maps.shape[4]) + 1),
        dtype=np.int32,
    ).reshape(-1)
    exported_beat_indexes = tuple(
        range(min(int(maps.shape[3]), SEGMENT_VELOCITY_MAP_EXPORTED_BEATS))
    )

    discovery_started = perf_counter()
    candidate_indexes = _candidate_segment_indexes(segments, maps)
    candidate_discovery_seconds = perf_counter() - discovery_started

    range_started = perf_counter()
    valid_indexes: list[tuple[int, int, int]] = []
    minimum = np.inf
    maximum = -np.inf
    for branch_index, radius_index in candidate_indexes:
        branch_id = (
            int(branch_ids[branch_index])
            if branch_index < branch_ids.size
            else branch_index + 1
        )
        tile_minimum, tile_maximum = _tile_median_image_range(
            maps,
            exported_beat_indexes,
            branch_index,
            radius_index,
        )
        if not np.isfinite(tile_minimum) or not np.isfinite(tile_maximum):
            continue
        valid_indexes.append((branch_index, radius_index, branch_id))
        minimum = min(minimum, tile_minimum)
        maximum = max(maximum, tile_maximum)
    range_scan_seconds = perf_counter() - range_started

    grid_side = max(1, int(math.ceil(math.sqrt(len(valid_indexes)))))
    prefix = "A" if vessel == "artery" else "V"
    tiles = tuple(
        _MosaicTile(
            branch_index=branch_index,
            radius_index=radius_index,
            branch_id=branch_id,
            row=tile_index // grid_side,
            column=tile_index % grid_side,
            label=f"{prefix}{branch_id}/R{radius_index + 1}",
        )
        for tile_index, (branch_index, radius_index, branch_id) in enumerate(
            valid_indexes
        )
    )
    tile_side = min(
        max(int(maps.shape[0]), int(maps.shape[1]), 1),
        max(2, SEGMENT_VELOCITY_MAP_MAX_SIDE // grid_side),
    )
    source = _MosaicSource(
        vessel=vessel,
        maps=maps,
        tiles=tiles,
        grid_side=grid_side,
        tile_side=tile_side,
        exported_beat_indexes=exported_beat_indexes,
        velocity_minimum=float(minimum),
        velocity_maximum=float(maximum),
        render_plan=_render_plan(maps.shape[:2], tiles, tile_side),
    )
    values_per_reduction = (
        int(maps.shape[0])
        * int(maps.shape[1])
        * int(maps.shape[2])
        * len(exported_beat_indexes)
        * len(candidate_indexes)
    )
    return source, _MosaicPreparationStats(
        candidate_count=len(candidate_indexes),
        candidate_discovery_seconds=candidate_discovery_seconds,
        range_scan_seconds=range_scan_seconds,
        range_values_inspected=values_per_reduction,
    )


def _candidate_segment_indexes(segments, maps: np.ndarray) -> list[tuple[int, int]]:
    branch_count = int(maps.shape[4])
    radius_count = int(maps.shape[5])
    masks = getattr(segments, "segment_masks", None)
    if masks is not None:
        mask_values = np.asarray(masks, dtype=bool)
        if (
            mask_values.ndim == 4
            and mask_values.shape[:2] == (radius_count, branch_count)
        ):
            present = np.any(mask_values, axis=(-2, -1))
            return [
                (branch_index, radius_index)
                for branch_index in range(branch_count)
                for radius_index in range(radius_count)
                if present[radius_index, branch_index]
            ]
    return [
        (branch_index, radius_index)
        for branch_index in range(branch_count)
        for radius_index in range(radius_count)
    ]


def _tile_median_image_range(
    maps: np.ndarray,
    beat_indexes: tuple[int, ...],
    branch_index: int,
    radius_index: int,
) -> tuple[float, float]:
    if not beat_indexes:
        return float("inf"), float("-inf")
    frames = np.concatenate(
        [
            maps[:, :, :, beat_index, branch_index, radius_index]
            for beat_index in beat_indexes
        ],
        axis=2,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        median_image = np.nanmedian(frames, axis=2)
    finite = median_image[np.isfinite(median_image)]
    if finite.size == 0:
        return float("inf"), float("-inf")
    return float(np.min(finite)), float(np.max(finite))


def _global_velocity_range(
    sources: tuple[_MosaicSource, ...],
) -> tuple[float, float]:
    minimum = np.inf
    maximum = -np.inf
    for source in sources:
        if np.isfinite(source.velocity_minimum):
            minimum = min(minimum, source.velocity_minimum)
        if np.isfinite(source.velocity_maximum):
            maximum = max(maximum, source.velocity_maximum)
    if not np.isfinite(minimum) or not np.isfinite(maximum):
        return 0.0, 1.0
    if maximum <= minimum:
        padding = max(abs(minimum) * 0.01, 0.5)
        return minimum - padding, maximum + padding
    return minimum, maximum


def _frame_indexes(time_count: int, beat_count: int):
    """Yield time as the fastest-changing movie dimension."""
    for beat_index in range(beat_count):
        for time_index in range(time_count):
            yield time_index, beat_index


def _mosaic_frame(
    source: _MosaicSource,
    time_index: int,
    beat_index: int,
    velocity_range: tuple[float, float],
    color_lut: np.ndarray,
) -> np.ndarray:
    return _mosaic_frame_chunk(
        source,
        np.asarray([time_index], dtype=np.intp),
        beat_index,
        velocity_range,
        color_lut,
    )[0]


def _mosaic_frame_chunk(
    source: _MosaicSource,
    time_indexes: np.ndarray,
    beat_index: int,
    velocity_range: tuple[float, float],
    color_lut: np.ndarray,
) -> np.ndarray:
    """Render a bounded time chunk after gathering all tiles contiguously."""
    times = np.asarray(time_indexes, dtype=np.intp).reshape(-1)
    frames = np.zeros(
        (times.size, source.video_side, source.video_side, 3),
        dtype=np.uint8,
    )
    if times.size == 0 or not source.tiles:
        return frames

    plan = source.render_plan
    display_xy_time_tile = source.maps[
        plan.x_indexes[:, None, None, None],
        plan.y_indexes[None, :, None, None],
        times[None, None, :, None],
        int(beat_index),
        plan.branch_indexes[None, None, None, :],
        plan.radius_indexes[None, None, None, :],
    ]
    display_time_tile_yx = display_xy_time_tile.transpose(2, 3, 1, 0)
    colored_tiles = _colorize_velocity(
        display_time_tile_yx,
        velocity_range,
        color_lut,
    )
    for tile_index, tile in enumerate(source.tiles):
        x_start = tile.column * source.tile_side + plan.x_offset
        y_start = tile.row * source.tile_side + plan.y_offset
        frames[
            :,
            y_start : y_start + plan.display_y_count,
            x_start : x_start + plan.display_x_count,
        ] = colored_tiles[:, tile_index]
    return frames


def _render_plan(
    map_shape: tuple[int, int],
    tiles: tuple[_MosaicTile, ...],
    tile_side: int,
) -> _MosaicRenderPlan:
    x_count, y_count = (int(value) for value in map_shape)
    scale = min(tile_side / x_count, tile_side / y_count, 1.0)
    display_x_count = max(1, int(round(x_count * scale)))
    display_y_count = max(1, int(round(y_count * scale)))
    x_indexes = np.rint(
        np.linspace(0, x_count - 1, display_x_count, dtype=np.float32)
    ).astype(np.intp)
    y_indexes = np.rint(
        np.linspace(0, y_count - 1, display_y_count, dtype=np.float32)
    ).astype(np.intp)
    return _MosaicRenderPlan(
        x_indexes=x_indexes,
        y_indexes=y_indexes,
        branch_indexes=np.asarray(
            [tile.branch_index for tile in tiles],
            dtype=np.intp,
        ),
        radius_indexes=np.asarray(
            [tile.radius_index for tile in tiles],
            dtype=np.intp,
        ),
        x_offset=(tile_side - display_x_count) // 2,
        y_offset=(tile_side - display_y_count) // 2,
        display_x_count=display_x_count,
        display_y_count=display_y_count,
    )


def _colorize_velocity(
    values: np.ndarray,
    velocity_range: tuple[float, float],
    color_lut: np.ndarray,
) -> np.ndarray:
    vmin, vmax = velocity_range
    finite = np.isfinite(values)
    result = np.zeros((*values.shape, 3), dtype=np.uint8)
    if not np.any(finite):
        return result
    scaled = np.clip((values[finite] - vmin) / (vmax - vmin), 0.0, 1.0)
    indexes = np.rint(scaled * (color_lut.shape[0] - 1)).astype(np.uint8)
    result[finite] = color_lut[indexes]
    return result


def _turbo_lut() -> np.ndarray:
    import matplotlib

    matplotlib.use("Agg", force=True)
    colors = matplotlib.colormaps[SEGMENT_VELOCITY_MAP_COLORMAP](
        np.linspace(0.0, 1.0, 256, dtype=np.float32)
    )[:, :3]
    return np.rint(colors * 255.0).astype(np.uint8)


def _video_metadata(
    source: _MosaicSource,
    velocity_range: tuple[float, float],
    fps: float,
) -> dict[str, object]:
    return {
        "title": f"EyeFlow {source.vessel} segment velocity-map mosaic",
        "artifact": "segment_velocity_map",
        "vessel": source.vessel,
        "codec": "MJPEG",
        "fps": float(fps),
        "velocity_unit": "mm/s",
        "velocity_range": [float(value) for value in velocity_range],
        "velocity_range_source": "temporal_median_image",
        "colormap": SEGMENT_VELOCITY_MAP_COLORMAP,
        "source_dimensions": ["x", "y", "time", "beat", "branch", "radius"],
        "frame_order": ["beat", "time"],
        "time_varies_fastest": True,
        "time_samples_per_beat": source.time_count,
        "source_beat_count": source.beat_count,
        "exported_beat_count": len(source.exported_beat_indexes),
        "exported_beat_indexes": list(source.exported_beat_indexes),
        "mosaic_grid": [source.grid_side, source.grid_side],
        "video_size_pixels": [source.video_side, source.video_side],
        "tile_size_pixels": source.tile_side,
        "tile_count": len(source.tiles),
        "tiles": [
            {
                "row": tile.row,
                "column": tile.column,
                "branch_index": tile.branch_index,
                "branch_id": tile.branch_id,
                "radius_index": tile.radius_index,
                "radius_number": tile.radius_index + 1,
                "label": tile.label,
            }
            for tile in source.tiles
        ],
    }


def _absolute_h5_path(path: str | None) -> str | None:
    if path is None:
        return None
    return f"/{path.lstrip('/')}"


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)
