"""Pipeline-context adapter for dense retinal motion estimation."""

from __future__ import annotations

import math
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import h5py
import numpy as np

from input_output.output_manager import OutputType
from input_output.schema import EyeFlowOutputPaths

from .calculator import create_retinal_motion_map
from .constants import DEFAULT_REGISTRATION_METHOD, RegistrationMethod
from .parameters import MotionMapConfig

DISPLACEMENT_MAP_PATH = EyeFlowOutputPaths.active().displacement_map
MAGNITUDE_VIDEO_FILENAME = "displacement_magnitude.mp4"
DEFAULT_MOMENT_PATH = "moment0"
DEFAULT_FPS = 25.0

MaskMode = Literal["combined", "artery", "vein", "labeled"]

ARTERY_MASK_PATH = "segmentation/Retina/artery_mask"
VEIN_MASK_PATH = "segmentation/Retina/vein_mask"
VESSEL_MASK_PATH = "segmentation/Retina/vessel_mask"
LABELED_VESSELS_PATH = "segmentation/Retina/labeled_vessels"


@dataclass(frozen=True, slots=True)
class DisplacementMapPipelineConfig:
    """Input choices owned by the pipeline rather than by the algorithm."""

    moment_path: str = DEFAULT_MOMENT_PATH
    mask_mode: MaskMode = "combined"
    fallback_fps: float = DEFAULT_FPS
    registration_method: RegistrationMethod = DEFAULT_REGISTRATION_METHOD


@dataclass(frozen=True, slots=True)
class DisplacementMapInputs:
    moment: h5py.Dataset
    mask: np.ndarray
    mask_source: str
    fps: float


def run_displacement_map(
    ctx,
    config: DisplacementMapPipelineConfig | None = None,
) -> None:
    """Estimate and persist the dense field plus its magnitude visualization."""

    selected = config or DisplacementMapPipelineConfig()
    inputs = load_displacement_map_inputs(ctx, selected)
    source_filename = inputs.moment.file.filename
    if not source_filename:
        raise ValueError("The HoloDoppler input does not have a filesystem path.")

    output_video = ctx.output.path_for(OutputType.MP4, MAGNITUDE_VIDEO_FILENAME)
    output_video.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="eyeflow-displacement-map-") as tmp_dir:
        algorithm_config = MotionMapConfig(
            input=Path(source_filename),
            output_dir=Path(tmp_dir),
            h5_dataset=inputs.moment.name.lstrip("/"),
            h5_fps=inputs.fps,
            registration_method=selected.registration_method,
            save_field=True,
        )
        outputs = create_retinal_motion_map(
            algorithm_config,
            analysis_mask_array=inputs.mask,
            magnitude_video_path=output_video,
        )
        displacement = np.load(outputs["displacement_field"], mmap_mode="r")
        try:
            ctx.output.h5.write(
                DISPLACEMENT_MAP_PATH,
                displacement,
                axes=["frame", "y", "x", "component"],
                components=["dx", "dy"],
                units="pixels",
                reference="displacement from each frame to the temporal mean",
                source_moment_path=inputs.moment.name,
                source_mask_path=inputs.mask_source,
                registration_method=selected.registration_method,
            )
            ctx.output.h5.flush()
        finally:
            mmap = getattr(displacement, "_mmap", None)
            if mmap is not None:
                mmap.close()

    ctx.log(f"Dense displacement map written to {DISPLACEMENT_MAP_PATH}.")
    ctx.log(f"Displacement magnitude video written to {output_video}.")


def load_displacement_map_inputs(
    ctx,
    config: DisplacementMapPipelineConfig,
) -> DisplacementMapInputs:
    """Resolve the root HD moment and aligned DV vessel mask."""

    ctx.require_inputs("hd", "dv")
    moment = resolve_moment_dataset(ctx.inputs.hd.h5.h5file, config.moment_path)
    spatial_shape = tuple(int(size) for size in moment.shape[-2:])
    mask, mask_source = resolve_retina_mask(
        ctx.inputs.dv.h5.h5file,
        spatial_shape,
        config.mask_mode,
    )
    fps = resolve_frame_rate(ctx, config.fallback_fps)
    return DisplacementMapInputs(moment, mask, mask_source, fps)


def resolve_moment_dataset(
    h5file: h5py.File | None,
    moment_path: str = DEFAULT_MOMENT_PATH,
) -> h5py.Dataset:
    """Return one 3-D moment dataset stored at the HoloDoppler root."""

    if h5file is None:
        raise ValueError("The HoloDoppler HDF5 input is required.")
    normalized = str(moment_path).replace("\\", "/").strip("/")
    if not normalized or "/" in normalized:
        raise ValueError("The HoloDoppler moment must be a root dataset name.")

    candidates = (normalized, "M0") if normalized == "moment0" else (normalized,)
    for candidate in candidates:
        found = h5file.get(candidate)
        if found is None:
            continue
        if not isinstance(found, h5py.Dataset) or found.ndim != 3:
            raise ValueError(
                f"HoloDoppler moment '{candidate}' must be a 3-D dataset, "
                f"got {getattr(found, 'shape', None)}."
            )
        return found
    raise KeyError(
        "Missing HoloDoppler root moment dataset. Tried: "
        + ", ".join(repr(candidate) for candidate in candidates)
    )


def resolve_retina_mask(
    h5file: h5py.File | None,
    spatial_shape: tuple[int, int],
    mode: MaskMode = "combined",
) -> tuple[np.ndarray, str]:
    """Resolve and align one DopplerView retina mask."""

    if h5file is None:
        raise ValueError("The DopplerView HDF5 input is required.")
    if mode == "artery":
        mask = _required_mask(h5file, ARTERY_MASK_PATH)
        source = ARTERY_MASK_PATH
    elif mode == "vein":
        mask = _required_mask(h5file, VEIN_MASK_PATH)
        source = VEIN_MASK_PATH
    elif mode == "labeled":
        mask = _required_mask(h5file, LABELED_VESSELS_PATH)
        source = LABELED_VESSELS_PATH
    elif mode == "combined":
        mask, source = _combined_mask(h5file)
    else:
        raise ValueError(f"Unknown displacement-map mask mode: {mode!r}")

    aligned = _align_mask(mask, spatial_shape, source)
    if not np.any(aligned):
        raise ValueError(f"DopplerView mask '{source}' contains no vessel pixels.")
    return aligned, source


def _combined_mask(h5file: h5py.File) -> tuple[np.ndarray, str]:
    vessel = _optional_mask(h5file, VESSEL_MASK_PATH)
    if vessel is not None:
        return vessel, VESSEL_MASK_PATH

    labeled = _optional_mask(h5file, LABELED_VESSELS_PATH)
    if labeled is not None:
        return labeled, LABELED_VESSELS_PATH

    artery = _optional_mask(h5file, ARTERY_MASK_PATH)
    vein = _optional_mask(h5file, VEIN_MASK_PATH)
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


def _required_mask(h5file: h5py.File, path: str) -> np.ndarray:
    mask = _optional_mask(h5file, path)
    if mask is None:
        raise KeyError(f"Missing DopplerView mask dataset at '{path}'.")
    return mask


def _optional_mask(h5file: h5py.File, path: str) -> np.ndarray | None:
    found = h5file.get(path)
    if found is None:
        return None
    if not isinstance(found, h5py.Dataset):
        raise ValueError(f"DopplerView mask path '{path}' is not a dataset.")
    return np.asarray(found[()])


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
    "DISPLACEMENT_MAP_PATH",
    "DisplacementMapPipelineConfig",
    "load_displacement_map_inputs",
    "resolve_moment_dataset",
    "resolve_retina_mask",
    "run_displacement_map",
]
