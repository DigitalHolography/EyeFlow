"""Create ImageJ-style Sobel gradient artifacts from flat-field moment0."""

from __future__ import annotations

import math
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage

from input_output.output_manager import OutputType
from input_output.writers.avi import MjpegAviWriter


AVI_FILENAME = "spatial_gradient_moment0.avi"
PNG_FILENAME = "spatial_gradient_moment0.png"
DEFAULT_FPS = 25.0
CONTRAST_HIGH_PERCENTILE = 99.9
CONTRAST_GAMMA = 1.0
STATE_KEY = "spatial_gradient_moment0_artifacts"


@dataclass(frozen=True, slots=True)
class SpatialGradientMoment0Artifacts:
    avi_path: Path
    mean_png_path: Path
    gradient_path: Path
    frame_count: int
    display_maximum: float
    temporary_directory: object

    def cleanup(self) -> None:
        cleanup = getattr(self.temporary_directory, "cleanup", None)
        if cleanup is not None:
            cleanup()


def spatial_gradient(frame) -> np.ndarray:
    """Return the magnitude of the horizontal and vertical 3x3 Sobel filters."""

    image = np.asarray(frame, dtype=np.float32)
    if image.ndim != 2:
        raise ValueError(f"A moment0 frame must be 2-D, got shape {image.shape}.")
    # ImageJ's Find Edges command combines the two Sobel responses as
    # sqrt(Gx**2 + Gy**2). Nearest-edge extension matches its border handling.
    finite_image = np.nan_to_num(image, nan=0.0, posinf=0.0, neginf=0.0)
    horizontal = ndimage.sobel(finite_image, axis=1, mode="nearest")
    vertical = ndimage.sobel(finite_image, axis=0, mode="nearest")
    return np.hypot(horizontal, vertical, dtype=np.float32)


def run_spatial_gradient_moment0(ctx) -> SpatialGradientMoment0Artifacts:
    """Export a globally scaled gradient AVI and its temporal-mean PNG."""

    ctx.require_inputs("hd")
    if not ctx.output.available:
        raise ValueError("An output manager is required to export gradient artifacts.")

    moment0ff = ctx.inputs.hd.as_holodoppler().moment0_flat_field_dataset()
    if moment0ff is None:
        raise KeyError(
            "Missing flat-field HoloDoppler moment0 dataset. Tried: "
            "'moment0ff', 'M0FF'."
        )
    frame_count, height, width = (int(size) for size in moment0ff.shape)
    if frame_count <= 0 or height <= 0 or width <= 0:
        raise ValueError(f"moment0ff must have non-empty dimensions, got {moment0ff.shape}.")

    # Preserve the quantitative float32 gradient cube for waveform_velocity.
    # Display-only contrast is applied later and never reaches profile data.
    temporary_directory = tempfile.TemporaryDirectory(
        prefix=".eyeflow-spatial-gradient-"
    )
    gradient_path = Path(temporary_directory.name) / "spatial_gradient_moment0.npy"
    gradient_video = np.lib.format.open_memmap(
        gradient_path,
        mode="w+",
        dtype=np.float32,
        shape=(frame_count, height, width),
    )
    gradient_sum = np.zeros((height, width), dtype=np.float64)
    observed_maximum = 0.0
    for frame_index in range(frame_count):
        gradient = spatial_gradient(moment0ff[frame_index])
        gradient_video[frame_index] = gradient
        gradient_sum += gradient
        frame_maximum = float(np.max(gradient, initial=0.0))
        observed_maximum = max(observed_maximum, frame_maximum)

    mean_gradient = (gradient_sum / float(frame_count)).astype(np.float32)
    gradient_video.flush()
    display_maximum = _contrast_maximum(mean_gradient, observed_maximum)

    fps = resolve_frame_rate(ctx)
    avi_path = ctx.output.path_for(OutputType.AVI, AVI_FILENAME)
    metadata = {
        "title": "EyeFlow spatial gradient of moment0ff",
        "artifact": "spatial_gradient_moment0",
        "algorithm": "Sobel gradient magnitude (ImageJ/Fiji Find Edges)",
        "source_dataset": "/moment0ff",
        "display_range": [0.0, display_maximum],
        "contrast_high_percentile": CONTRAST_HIGH_PERCENTILE,
        "contrast_gamma": CONTRAST_GAMMA,
        "fps": fps,
        "frame_count": frame_count,
    }
    with MjpegAviWriter(
        avi_path,
        width=width,
        height=height,
        fps=fps,
        metadata=metadata,
    ) as video:
        display_sum = np.zeros((height, width), dtype=np.float64)
        for frame_index in range(frame_count):
            display_frame = _display_frame(
                gradient_video[frame_index], display_maximum
            )
            display_sum += display_frame
            video.write_frame(display_frame)

    mean_display_frame = np.rint(display_sum / float(frame_count)).astype(np.uint8)
    mean_png_path = ctx.output.write_png(
        mean_display_frame,
        PNG_FILENAME,
    )
    artifacts = SpatialGradientMoment0Artifacts(
        avi_path=avi_path,
        mean_png_path=mean_png_path,
        gradient_path=gradient_path,
        frame_count=frame_count,
        display_maximum=display_maximum,
        temporary_directory=temporary_directory,
    )
    del gradient_video
    ctx.state.set(STATE_KEY, artifacts)
    ctx.log(
        "Exported moment0 spatial gradients: "
        f"{avi_path} ({frame_count} frames) and {mean_png_path}."
    )
    return artifacts


def _display_frame(gradient: np.ndarray, maximum: float) -> np.ndarray:
    if not math.isfinite(maximum) or maximum <= 0.0:
        return np.zeros(gradient.shape, dtype=np.uint8)
    scaled = np.clip(np.asarray(gradient, dtype=np.float32) / maximum, 0.0, 1.0)
    scaled = np.power(scaled, CONTRAST_GAMMA, dtype=np.float32)
    return np.rint(scaled * 255.0).astype(np.uint8)


def _contrast_maximum(mean_gradient: np.ndarray, fallback: float) -> float:
    values = np.asarray(mean_gradient)
    positive_finite = values[np.isfinite(values) & (values > 0.0)]
    if positive_finite.size:
        maximum = float(
            np.percentile(positive_finite, CONTRAST_HIGH_PERCENTILE)
        )
        if math.isfinite(maximum) and maximum > 0.0:
            return maximum
    return float(fallback)


def resolve_frame_rate(ctx, fallback: float = DEFAULT_FPS) -> float:
    """Resolve the HoloDoppler acquisition rate, with an older-file fallback."""

    try:
        dt_seconds = float(ctx.inputs.hd.as_holodoppler().timing().dt_seconds)
        fps = 1.0 / dt_seconds
        if math.isfinite(fps) and fps > 0.0:
            return fps
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        pass
    fallback = float(fallback)
    if not math.isfinite(fallback) or fallback <= 0.0:
        raise ValueError("fallback fps must be finite and greater than zero.")
    ctx.log_warning(f"Could not resolve HoloDoppler timing; using {fallback:g} fps.")
    return fallback


__all__ = [
    "AVI_FILENAME",
    "CONTRAST_GAMMA",
    "CONTRAST_HIGH_PERCENTILE",
    "PNG_FILENAME",
    "STATE_KEY",
    "SpatialGradientMoment0Artifacts",
    "resolve_frame_rate",
    "run_spatial_gradient_moment0",
    "spatial_gradient",
]
