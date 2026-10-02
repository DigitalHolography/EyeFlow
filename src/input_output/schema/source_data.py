"""Central source-data models shared across EyeFlow pipelines."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from calculations.topology import OpticDisc


@dataclass(frozen=True)
class HolodopplerTiming:
    """Temporal sampling metadata exported by Holodoppler."""

    sampling_freq: float
    batch_stride: float

    @property
    def dt_seconds(self) -> float:
        return self.batch_stride / self.sampling_freq


@dataclass(frozen=True)
class PixelPitch:
    """Native spatial sampling exported by Holodoppler, in ``(x, y)`` order."""

    x_m: float
    y_m: float

    def __post_init__(self) -> None:
        values = np.asarray((self.x_m, self.y_m), dtype=np.float64)
        if not np.all(np.isfinite(values)) or np.any(values <= 0.0):
            raise ValueError("Holodoppler pixel pitch must contain two positive values.")
        object.__setattr__(self, "x_m", float(values[0]))
        object.__setattr__(self, "y_m", float(values[1]))

    @property
    def xy_m(self) -> tuple[float, float]:
        return self.x_m, self.y_m

    @property
    def isotropic_m(self) -> float:
        """Return the common pitch accepted by today's scalar consumers."""

        if not np.isclose(self.x_m, self.y_m, rtol=1e-6, atol=0.0):
            raise ValueError(
                "EyeFlow currently requires approximately equal Holodoppler "
                f"x/y pixel pitches, got {(self.x_m, self.y_m)} m."
            )
        return float(np.mean((self.x_m, self.y_m)))

    @property
    def isotropic_mm(self) -> float:
        """Return the common pitch in millimetres for profile calculations."""

        return self.isotropic_m * 1e3


@dataclass(frozen=True)
class ImageMaps:
    """Lazy Holodoppler image-map datasets in the aligned analysis frame."""

    moment0: object
    moment2: object


@dataclass(frozen=True)
class VesselMasks:
    """Aligned DopplerView vessel masks used by retinal analysis.

    ``velocity_background`` may retain the original combined vessel support
    when a vessel is intentionally disabled for measurement.
    """

    artery: np.ndarray
    vein: np.ndarray
    labeled: np.ndarray | None = None
    velocity_background: np.ndarray | None = None


@dataclass(frozen=True)
class RetinalSegmentation:
    """Aligned vessel and optic-disc segmentation."""

    vessels: VesselMasks
    optic_disc: OpticDisc


@dataclass(frozen=True)
class HolodopplerMetadata:
    """Physical and temporal Holodoppler acquisition metadata."""

    timing: HolodopplerTiming
    pixel_pitch: PixelPitch


@dataclass(frozen=True)
class DopplerViewMetadata:
    """DopplerView settings and alignment metadata used by EyeFlow."""

    local_background_dist: int
    spatial_axes_swapped_to_match_hd: bool


@dataclass(frozen=True)
class RetinalSourceData:
    """Canonical aligned source information shared by retinal pipelines."""

    image_maps: ImageMaps
    segmentation: RetinalSegmentation
    holodoppler: HolodopplerMetadata
    doppler_view: DopplerViewMetadata

