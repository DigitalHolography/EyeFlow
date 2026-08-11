"""Validated parameter models for displacement-map estimation."""

from dataclasses import dataclass
from pathlib import Path

from .constants import PhotometricMode, RegistrationMethod


@dataclass(frozen=True, slots=True)
class MotionMapConfig:
    """Configuration for one displacement-map calculation."""

    input: Path
    output_dir: Path
    analysis_mask: Path | None = None
    h5_dataset: str = "moment0"
    h5_frame_axis: int = 0
    registration_method: RegistrationMethod = "symmetric_forces_demons"
    iterations: int = 20
    scale: float = 0.5
    registration_initialization: str = "zero"
    registration_metric_radius: int = 4
    registration_learning_rate: float = 1.0
    bspline_mesh_size: int = 8
    temporal_median_window: int = 3
    temporal_alpha: float = 1.0
    photometric_normalization: PhotometricMode = "hybrid"
    illumination_sigma: float = 8.0
    local_contrast_sigma: float = 4.0
    photometric_confidence_floor: float = 0.30
    save_field: bool = True
    max_frames: int | None = None
    invert_analysis_mask: bool = False
    mask_feather_sigma: float = 3.0
    h5_fps: float = 25.0
    h5_low_percentile: float = 1.0
    h5_high_percentile: float = 99.5
    field_sigma: float = 1.0
    update_sigma: float = 0.0
    photometric_clip: float = 3.0
    photometric_dark_percentile: float = 2.0
    structure_edge_weight: float = 0.25
    normalization: str = "global-minmax"
    low_percentile: float = 1.0
    high_percentile: float = 99.5
    fixed_max_px: float = 3.0
    gamma: float = 1.0
    visualization_sigma: float = 0.0
    codec: str = "mp4v"

    def __post_init__(self) -> None:
        if self.iterations <= 0:
            raise ValueError("iterations must be greater than zero")
        if not 0.0 < self.scale <= 1.0:
            raise ValueError("scale must be in (0, 1]")
        if self.registration_initialization not in {"zero", "previous"}:
            raise ValueError("registration_initialization must be 'zero' or 'previous'")
        if self.temporal_median_window not in {1, 3, 5, 7}:
            raise ValueError("temporal_median_window must be one of 1, 3, 5, or 7")
        if not 0.0 < self.temporal_alpha <= 1.0:
            raise ValueError("temporal_alpha must be in (0, 1]")
        if self.illumination_sigma < 0 or self.local_contrast_sigma < 0:
            raise ValueError("photometric sigmas must be non-negative")
        if not 0.0 <= self.photometric_confidence_floor <= 1.0:
            raise ValueError("photometric_confidence_floor must be in [0, 1]")
        if self.registration_metric_radius <= 0:
            raise ValueError("registration_metric_radius must be greater than zero")
        if self.registration_learning_rate <= 0:
            raise ValueError("registration_learning_rate must be greater than zero")
        if self.bspline_mesh_size <= 0:
            raise ValueError("bspline_mesh_size must be greater than zero")
        if self.max_frames is not None and self.max_frames <= 0:
            raise ValueError("max_frames must be greater than zero")


@dataclass(frozen=True, slots=True)
class PhotometricConfig:
    mode: PhotometricMode
    illumination_sigma: float
    local_contrast_sigma: float
    clip: float
    confidence_floor: float
    dark_percentile: float
