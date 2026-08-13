"""AVI mosaics for rotated per-segment velocity maps."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
import math
import warnings

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from input_output.writers.avi import AviArtifactWriter


SEGMENT_VELOCITY_MAP_AVI_FOLDER = "segment_velocity_map"
SEGMENT_VELOCITY_MAP_FPS = 60.0
SEGMENT_VELOCITY_MAP_COLORMAP = "turbo"
SEGMENT_VELOCITY_MAP_MAX_SIDE = 1024
SEGMENT_VELOCITY_MAP_EXPORTED_BEATS = 1


@dataclass(frozen=True)
class _MosaicTile:
    branch_index: int
    radius_index: int
    branch_id: int
    row: int
    column: int
    label: str


@dataclass(frozen=True)
class _MosaicSource:
    vessel: str
    maps: np.ndarray
    tiles: tuple[_MosaicTile, ...]
    grid_side: int
    tile_side: int

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
    schema = _resolve_output_paths(output_paths)
    sources = (
        _mosaic_source(
            "artery",
            artery_segments,
            packed_segment_maps.get(
                schema.artery_segments.velocity_map_per_segment,
            ),
        ),
        _mosaic_source(
            "vein",
            vein_segments,
            packed_segment_maps.get(
                schema.vein_segments.velocity_map_per_segment,
            ),
        ),
    )
    available_sources = tuple(source for source in sources if source is not None)
    if not available_sources:
        return []

    velocity_range = _global_velocity_range(available_sources)
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
            for time_index, beat_index in _frame_indexes(
                source.time_count,
                min(source.beat_count, SEGMENT_VELOCITY_MAP_EXPORTED_BEATS),
            ):
                video.write_frame(
                    _mosaic_frame(
                        source,
                        time_index,
                        beat_index,
                        velocity_range,
                        color_lut,
                    )
                )
        return str(video.path)

    if len(available_sources) == 1:
        return [export_source(available_sources[0])]
    with ThreadPoolExecutor(
        max_workers=min(2, len(available_sources)),
        thread_name_prefix="segment-avi",
    ) as executor:
        return list(executor.map(export_source, available_sources))


def _mosaic_source(vessel: str, segments, packed_value) -> _MosaicSource | None:
    if packed_value is None:
        return None
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
    valid_indexes: list[tuple[int, int, int]] = []
    for branch_index in range(maps.shape[4]):
        branch_id = (
            int(branch_ids[branch_index])
            if branch_index < branch_ids.size
            else branch_index + 1
        )
        for radius_index in range(maps.shape[5]):
            if _has_finite_value(maps[:, :, :, :, branch_index, radius_index]):
                valid_indexes.append((branch_index, radius_index, branch_id))

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
    return _MosaicSource(
        vessel=vessel,
        maps=maps,
        tiles=tiles,
        grid_side=grid_side,
        tile_side=min(
            max(int(maps.shape[0]), int(maps.shape[1]), 1),
            max(2, SEGMENT_VELOCITY_MAP_MAX_SIDE // grid_side),
        ),
    )


def _global_velocity_range(
    sources: tuple[_MosaicSource, ...],
) -> tuple[float, float]:
    minimum = np.inf
    maximum = -np.inf
    for source in sources:
        for tile in source.tiles:
            values = source.maps[
                :, :, :, :, tile.branch_index, tile.radius_index
            ]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                tile_minimum = float(np.nanmin(values))
                tile_maximum = float(np.nanmax(values))
            if np.isfinite(tile_minimum):
                minimum = min(minimum, tile_minimum)
            if np.isfinite(tile_maximum):
                maximum = max(maximum, tile_maximum)
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
    frame = np.zeros(
        (source.video_side, source.video_side, 3),
        dtype=np.uint8,
    )
    x_count, y_count = source.maps.shape[:2]
    scale = min(source.tile_side / x_count, source.tile_side / y_count, 1.0)
    display_x_count = max(1, int(round(x_count * scale)))
    display_y_count = max(1, int(round(y_count * scale)))
    x_indexes = np.rint(
        np.linspace(0, x_count - 1, display_x_count, dtype=np.float32)
    ).astype(np.int32)
    y_indexes = np.rint(
        np.linspace(0, y_count - 1, display_y_count, dtype=np.float32)
    ).astype(np.int32)
    x_offset = (source.tile_side - display_x_count) // 2
    y_offset = (source.tile_side - display_y_count) // 2
    for tile in source.tiles:
        x_start = tile.column * source.tile_side + x_offset
        y_start = tile.row * source.tile_side + y_offset
        velocity_xy = source.maps[
            :,
            :,
            time_index,
            beat_index,
            tile.branch_index,
            tile.radius_index,
        ]
        display_yx = velocity_xy.T[np.ix_(y_indexes, x_indexes)]
        frame[
            y_start : y_start + display_y_count,
            x_start : x_start + display_x_count,
        ] = _colorize_velocity(display_yx, velocity_range, color_lut)
    return frame


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
        "colormap": SEGMENT_VELOCITY_MAP_COLORMAP,
        "source_dimensions": ["x", "y", "time", "beat", "branch", "radius"],
        "frame_order": ["beat", "time"],
        "time_varies_fastest": True,
        "time_samples_per_beat": source.time_count,
        "source_beat_count": source.beat_count,
        "exported_beat_count": min(
            source.beat_count,
            SEGMENT_VELOCITY_MAP_EXPORTED_BEATS,
        ),
        "exported_beat_indexes": (
            [0] if source.beat_count > 0 else []
        ),
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


def _has_finite_value(values: np.ndarray) -> bool:
    return bool(np.any(np.isfinite(values)))


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
