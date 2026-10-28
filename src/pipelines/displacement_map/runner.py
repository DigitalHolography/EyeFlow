"""Pipeline-context adapter for dense retinal motion estimation."""

from __future__ import annotations

import math
import tempfile
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

import h5py
import numpy as np

from input_output.output_manager import OutputType
from input_output.readers.frame_sequences import FrameSequence
from input_output.readers.h5 import RawH5SourceReader
from input_output.schema import DopplerViewSource
from input_output.schema.doppler_view import (
    DV_ARTERY_MASK_PATH as ARTERY_MASK_PATH,
    DV_LABELED_VESSELS_PATH as LABELED_VESSELS_PATH,
    DV_VEIN_MASK_PATH as VEIN_MASK_PATH,
    DV_VESSEL_MASK_PATH as VESSEL_MASK_PATH,
)

from .calculator import create_retinal_motion_map
from .constants import DEFAULT_REGISTRATION_METHOD, RegistrationMethod
from .parameters import MotionMapConfig

DISPLACEMENT_MAP_STATE = "displacement_map_artifacts"
MAGNITUDE_VIDEO_FILENAME = "displacement_magnitude.avi"
DISPLACEMENT_AVI_SUBFOLDER = "displacement_maps"
DEFAULT_MOMENT_PATH = "moment0"
DEFAULT_FPS = 25.0

MaskMode = Literal["both", "combined", "artery", "vein", "labeled"]

@dataclass(frozen=True, slots=True)
class DisplacementMapPipelineConfig:
    """Input choices owned by the pipeline rather than by the algorithm."""

    moment_path: str = DEFAULT_MOMENT_PATH
    mask_mode: MaskMode = "both"
    fallback_fps: float = DEFAULT_FPS
    registration_method: RegistrationMethod = DEFAULT_REGISTRATION_METHOD


@dataclass(frozen=True, slots=True)
class DisplacementMaskInput:
    name: str
    vessels: tuple[str, ...]
    mask: np.ndarray
    source: str


@dataclass(frozen=True, slots=True)
class DisplacementMapInputs:
    moment: h5py.Dataset
    masks: tuple[DisplacementMaskInput, ...]
    fps: float


@dataclass(slots=True)
class DisplacementMapArtifacts:
    registration_method: str
    field_paths_by_vessel: dict[str, Path]
    temporary_directory: object

    def cleanup(self) -> None:
        cleanup = getattr(self.temporary_directory, "cleanup", None)
        if cleanup is not None:
            cleanup()


def run_displacement_map(
    ctx,
    config: DisplacementMapPipelineConfig | None = None,
) -> None:
    """Estimate vessel-specific fields for downstream use without HDF5 persistence."""

    selected = config or DisplacementMapPipelineConfig()
    inputs = load_displacement_map_inputs(ctx, selected)
    source_filename = inputs.moment.file.filename
    if not source_filename:
        raise ValueError("The HoloDoppler input does not have a filesystem path.")

    temporary_directory = tempfile.TemporaryDirectory(
        prefix="eyeflow-displacement-map-"
    )
    field_paths_by_vessel: dict[str, Path] = {}
    output_videos: list[Path] = []
    try:
        base_config = MotionMapConfig(
            input=Path(source_filename),
            output_dir=Path(temporary_directory.name),
            h5_dataset=inputs.moment.name.lstrip("/"),
            h5_fps=inputs.fps,
            registration_method=selected.registration_method,
            save_field=True,
        )
        sequence = FrameSequence(
            base_config.input, base_config.h5_dataset, base_config.h5_frame_axis,
            base_config.h5_fps, base_config.h5_low_percentile,
            base_config.h5_high_percentile, h5_source=inputs.moment,
            defer_h5_range=True,
        )
        with sequence.cached_for_replay(
            base_config.max_frames, Path(temporary_directory.name)
        ) as replay:
            for mask_input in inputs.masks:
                _run_displacement_mask(
                    ctx, mask_input, inputs, base_config, replay,
                    Path(temporary_directory.name), field_paths_by_vessel, output_videos,
                )
    except Exception:
        temporary_directory.cleanup()
        raise

    previous = ctx.state.get(DISPLACEMENT_MAP_STATE)
    if isinstance(previous, DisplacementMapArtifacts):
        previous.cleanup()
    ctx.state.set(
        DISPLACEMENT_MAP_STATE,
        DisplacementMapArtifacts(
            registration_method=selected.registration_method,
            field_paths_by_vessel=field_paths_by_vessel,
            temporary_directory=temporary_directory,
        ),
    )
    ctx.log("Dense displacement maps prepared for velocity analysis processing.")
    for output_video in output_videos:
        ctx.log(f"Displacement magnitude video written to {output_video}.")


def _run_displacement_mask(
    ctx, mask_input, inputs, base_config, replay, temporary_path,
    field_paths_by_vessel, output_videos,
) -> None:
    output_dir = temporary_path / mask_input.name
    output_dir.mkdir(parents=True, exist_ok=True)
    video_filename = (
        MAGNITUDE_VIDEO_FILENAME
        if len(inputs.masks) == 1
        else f"{mask_input.name}_{MAGNITUDE_VIDEO_FILENAME}"
    )
    output_video = ctx.output.path_for(
        OutputType.AVI, f"{DISPLACEMENT_AVI_SUBFOLDER}/{video_filename}",
    )
    output_video.parent.mkdir(parents=True, exist_ok=True)
    outputs = create_retinal_motion_map(
        replace(base_config, output_dir=output_dir),
        analysis_mask_array=mask_input.mask,
        magnitude_video_path=output_video,
        h5_source=inputs.moment,
        frame_replay=replay,
    )
    field_path = Path(outputs["displacement_field"])
    for vessel in mask_input.vessels:
        field_paths_by_vessel[vessel] = field_path
    output_videos.append(output_video)


def load_displacement_map_inputs(
    ctx,
    config: DisplacementMapPipelineConfig,
) -> DisplacementMapInputs:
    """Resolve the root HD moment and aligned DV vessel mask."""

    ctx.require_inputs("hd", "dv")
    moment = ctx.inputs.hd.as_holodoppler().named_root_moment_dataset(config.moment_path)
    spatial_shape = tuple(int(size) for size in moment.shape[-2:])
    masks = resolve_retina_masks(
        ctx.inputs.dv.h5.h5file,
        spatial_shape,
        config.mask_mode,
        source=ctx.inputs.dv.as_dopplerview(),
    )
    fps = resolve_frame_rate(ctx, config.fallback_fps)
    return DisplacementMapInputs(moment, masks, fps)


def resolve_retina_masks(
    h5file: h5py.File | None,
    spatial_shape: tuple[int, int],
    mode: MaskMode,
    *,
    source: DopplerViewSource | None = None,
) -> tuple[DisplacementMaskInput, ...]:
    """Resolve one shared mask or separate artery and vein masks."""

    if h5file is None:
        raise ValueError("The DopplerView HDF5 input is required.")
    source = source or DopplerViewSource(RawH5SourceReader(h5file=h5file, label="DV"), {})
    if mode == "both":
        masks: list[DisplacementMaskInput] = []
        for vessel in ("artery", "vein"):
            mask, resolved_source = resolve_retina_mask(
                h5file,
                spatial_shape,
                vessel,
                source=source,
            )
            masks.append(
                DisplacementMaskInput(
                    name=vessel,
                    vessels=(vessel,),
                    mask=mask,
                    source=resolved_source,
                )
            )
        return tuple(masks)

    mask, resolved_source = resolve_retina_mask(h5file, spatial_shape, mode, source=source)
    if mode == "artery":
        vessels = ("artery",)
    elif mode == "vein":
        vessels = ("vein",)
    else:
        vessels = ("artery", "vein")
    return (
        DisplacementMaskInput(
            name=mode,
            vessels=vessels,
            mask=mask,
            source=resolved_source,
        ),
    )


def resolve_retina_mask(
    h5file: h5py.File | None,
    spatial_shape: tuple[int, int],
    mode: MaskMode = "combined",
    *,
    source: DopplerViewSource | None = None,
) -> tuple[np.ndarray, str]:
    """Resolve and align one DopplerView retina mask."""

    if h5file is None:
        raise ValueError("The DopplerView HDF5 input is required.")
    source = source or DopplerViewSource(RawH5SourceReader(h5file=h5file, label="DV"), {})
    if mode == "artery":
        mask = _required_mask(source, ARTERY_MASK_PATH)
        source_path = ARTERY_MASK_PATH
    elif mode == "vein":
        mask = _required_mask(source, VEIN_MASK_PATH)
        source_path = VEIN_MASK_PATH
    elif mode == "labeled":
        mask = _required_mask(source, LABELED_VESSELS_PATH)
        source_path = LABELED_VESSELS_PATH
    elif mode == "combined":
        mask, source_path = _combined_mask(source)
    else:
        raise ValueError(f"Unknown displacement-map mask mode: {mode!r}")

    aligned = _align_mask(mask, spatial_shape, source_path)
    if not np.any(aligned):
        raise ValueError(f"DopplerView mask '{source_path}' contains no vessel pixels.")
    return aligned, source_path


def _combined_mask(source: DopplerViewSource) -> tuple[np.ndarray, str]:
    vessel = _optional_mask(source, VESSEL_MASK_PATH)
    if vessel is not None:
        return vessel, VESSEL_MASK_PATH

    labeled = _optional_mask(source, LABELED_VESSELS_PATH)
    if labeled is not None:
        return labeled, LABELED_VESSELS_PATH

    artery = _optional_mask(source, ARTERY_MASK_PATH)
    vein = _optional_mask(source, VEIN_MASK_PATH)
    available = [mask for mask in (artery, vein) if mask is not None]
    if not available:
        raise KeyError(
            "Missing DopplerView retinal vessel mask. Tried: "
            f"{VESSEL_MASK_PATH}, {LABELED_VESSELS_PATH}, "
            f"{ARTERY_MASK_PATH}, {VEIN_MASK_PATH}"
        )
    combined = np.asarray(available[0]) != 0
    for mask in available[1:]:
        if np.asarray(mask).shape != combined.shape:
            raise ValueError("DopplerView artery and vein masks have different shapes.")
        combined |= np.asarray(mask) != 0
    sources = [
        path
        for path, mask in ((ARTERY_MASK_PATH, artery), (VEIN_MASK_PATH, vein))
        if mask is not None
    ]
    return combined, "+".join(sources)


def _required_mask(source: DopplerViewSource, path: str) -> np.ndarray:
    mask = _optional_mask(source, path)
    if mask is None:
        raise KeyError(f"Missing DopplerView mask dataset at '{path}'.")
    return mask


def _optional_mask(source: DopplerViewSource, path: str) -> np.ndarray | None:
    return source.retinal_mask(path)


def _align_mask(
    mask: np.ndarray,
    spatial_shape: tuple[int, int],
    source: str,
) -> np.ndarray:
    value = np.squeeze(np.asarray(mask))
    if value.ndim != 2:
        raise ValueError(
            f"DopplerView mask '{source}' must become 2-D after squeeze, got {value.shape}."
        )
    if value.shape == spatial_shape:
        return np.asarray(value != 0, dtype=bool)
    if value.T.shape == spatial_shape:
        return np.asarray(value.T != 0, dtype=bool)
    raise ValueError(
        f"DopplerView mask '{source}' shape {value.shape} does not match "
        f"HoloDoppler spatial shape {spatial_shape} in either axis order."
    )


def resolve_frame_rate(ctx, fallback: float = DEFAULT_FPS) -> float:
    """Use HoloDoppler timing metadata, falling back for older exports."""

    try:
        dt_seconds = float(ctx.inputs.hd.as_holodoppler().timing().dt_seconds)
        fps = 1.0 / dt_seconds
        if math.isfinite(fps) and fps > 0:
            return fps
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        pass
    fallback = float(fallback)
    if not math.isfinite(fallback) or fallback <= 0:
        raise ValueError("fallback_fps must be finite and greater than zero.")
    ctx.log_warning(f"Could not resolve HoloDoppler timing; using {fallback:g} fps.")
    return fallback


__all__ = [
    "DISPLACEMENT_MAP_STATE",
    "DisplacementMapArtifacts",
    "DisplacementMapPipelineConfig",
    "load_displacement_map_inputs",
    "resolve_retina_mask",
    "resolve_retina_masks",
    "run_displacement_map",
]
