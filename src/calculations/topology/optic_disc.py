"""Authoritative optic-disc geometry from DopplerView."""

from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real

import numpy as np

from .geometry import (
    AnnulusGeometry,
    annulus_count,
    _validated_center_xy,
    _validated_image_shape,
    image_half_diagonal,
)

DEFAULT_OPTIC_DISC_DIAMETER_MM = 1.91


@dataclass(frozen=True, eq=False)
class OpticDisc:
    """DopplerView optic-disc data in one spatial coordinate frame."""

    mask: np.ndarray | None
    center: tuple[float, float]
    width: float | None
    height: float | None
    is_fallback: bool = False

    def __post_init__(self) -> None:
        center = _validated_center_xy(self.center)
        width = _optional_positive_scalar(self.width, "width")
        height = _optional_positive_scalar(self.height, "height")
        if (width is None) != (height is None):
            raise ValueError("optic-disc width and height must both be present or absent.")

        mask = None
        if self.mask is not None:
            mask = np.asarray(self.mask, dtype=bool)
            if mask.ndim != 2:
                raise ValueError(f"optic-disc mask must be 2-D, got shape {mask.shape}.")
            mask = mask.copy()
            mask.setflags(write=False)

        object.__setattr__(self, "mask", mask)
        object.__setattr__(self, "center", center)
        object.__setattr__(self, "width", width)
        object.__setattr__(self, "height", height)

    @classmethod
    def from_measurements(
        cls,
        mask,
        center,
        width,
        height,
        image_shape: tuple[int, int],
        *,
        pixel_pitch_mm: float | None = None,
    ) -> OpticDisc:
        """Resolve incomplete DopplerView measurements for one image frame.

        A finite, nonempty mask supplies any missing center or dimensions.  If
        the mask is invalid, or neither a usable mask nor complete geometry is
        available, a centered circular disc of diameter 1.91 mm is used.
        The fallback requires the native isotropic pixel pitch in mm/pixel.
        """

        ny, nx = _validated_image_shape(image_shape)
        resolved_mask, invalid_mask = _resolved_source_mask(mask, (ny, nx))
        resolved_center = _valid_center_or_none(center)
        resolved_width = _valid_positive_scalar_or_none(width)
        resolved_height = _valid_positive_scalar_or_none(height)

        if resolved_mask is not None:
            y, x = np.nonzero(resolved_mask)
            if resolved_center is None:
                resolved_center = (float(np.mean(x)), float(np.mean(y)))
            if resolved_width is None:
                resolved_width = float(x.max() - x.min() + 1)
            if resolved_height is None:
                resolved_height = float(y.max() - y.min() + 1)
            return cls(
                resolved_mask,
                resolved_center,
                resolved_width,
                resolved_height,
            )

        if (
            not invalid_mask
            and resolved_center is not None
            and resolved_width is not None
            and resolved_height is not None
        ):
            return cls(
                None,
                resolved_center,
                resolved_width,
                resolved_height,
            )

        return cls.centered_fallback((ny, nx), pixel_pitch_mm=pixel_pitch_mm)

    @classmethod
    def centered_fallback(
        cls,
        image_shape: tuple[int, int],
        *,
        pixel_pitch_mm: float | None,
    ) -> OpticDisc:
        """Return a clipped, image-centered circle with a 1.91 mm diameter."""
        ny, nx = _validated_image_shape(image_shape)
        pitch = _optional_positive_scalar(pixel_pitch_mm, "pixel pitch in mm")
        if pitch is None:
            raise ValueError("Optic-disc fallback requires pixel_pitch_mm in mm/pixel.")
        diameter = _optional_positive_scalar(
            DEFAULT_OPTIC_DISC_DIAMETER_MM / pitch, "fallback diameter in pixels"
        )
        radius = diameter / 2.0
        center_x, center_y = nx / 2.0, ny / 2.0
        y, x = np.ogrid[:ny, :nx]
        fallback_mask = (y - center_y) ** 2 + (x - center_x) ** 2 <= radius**2
        return cls(
            fallback_mask,
            (center_x, center_y),
            diameter,
            diameter,
            is_fallback=True,
        )

    def transposed(self) -> OpticDisc:
        """Return this geometry after swapping its spatial axes."""

        center_x, center_y = self.center
        return OpticDisc(
            mask=None if self.mask is None else self.mask.T,
            center=(center_y, center_x),
            width=self.height,
            height=self.width,
            is_fallback=self.is_fallback,
        )

    def mask_for(self, image_shape: tuple[int, int]) -> np.ndarray:
        """Return the supplied mask or reconstruct its ellipse for ``image_shape``."""

        ny, nx = _validated_image_shape(image_shape)
        if self.mask is not None:
            if self.mask.shape != (ny, nx):
                raise ValueError(
                    "optic-disc mask must have the same orientation and shape as "
                    f"the image: expected {(ny, nx)}, got {self.mask.shape}."
                )
            return self.mask.copy()
        width, height = self._required_dimensions("reconstruct the optic-disc mask")
        center_x, center_y = self.center
        y, x = np.indices((ny, nx), dtype=np.float32)
        x_radius = np.float32(width / 2.0)
        y_radius = np.float32(height / 2.0)
        return ((x - np.float32(center_x)) / x_radius) ** 2 + (
            (y - np.float32(center_y)) / y_radius
        ) ** 2 <= 1.0

    def subtract_from(self, vessel_mask) -> np.ndarray:
        """Return a boolean vessel mask with optic-disc pixels removed."""

        vessel = np.asarray(vessel_mask)
        if vessel.ndim != 2 or vessel.dtype != np.bool_:
            raise TypeError("vessel_mask must be a 2-D boolean array.")
        return vessel & ~self.mask_for(vessel.shape)

    def centered_circle_radius_pixels(
        self,
        *,
        fallback_radius_pixels: float | None = None,
    ) -> int:
        """Return the integer radius used by topology around the optic disc.

        DopplerView dimensions are authoritative when present.  The fallback
        keeps mask-only callers usable when they provide an explicit annulus
        geometry.
        """

        if self.width is not None and self.height is not None:
            radius = min(self.width, self.height) / 2.0
        else:
            radius = _optional_nonnegative_scalar(
                fallback_radius_pixels,
                "fallback radius",
            )
            if radius is None:
                raise ValueError(
                    "optic-disc width and height or a fallback radius are "
                    "required to derive the centered circular topology mask."
                )
        return int(np.ceil(radius))

    def centered_circle_mask_for(
        self,
        image_shape: tuple[int, int],
        *,
        fallback_radius_pixels: float | None = None,
    ) -> np.ndarray:
        """Return the clipped centered circular mask used by topology."""

        ny, nx = _validated_image_shape(image_shape)
        radius = self.centered_circle_radius_pixels(
            fallback_radius_pixels=fallback_radius_pixels,
        )
        center_x, center_y = self.center
        y, x = np.ogrid[:ny, :nx]
        return (y - float(center_y)) ** 2 + (x - float(center_x)) ** 2 <= float(radius**2)

    def subtract_centered_circle_from(
        self,
        vessel_mask,
        *,
        fallback_radius_pixels: float | None = None,
    ) -> np.ndarray:
        """Return a vessel mask with the topology's circular disc removed."""

        vessel = np.asarray(vessel_mask)
        if vessel.ndim != 2 or vessel.dtype != np.bool_:
            raise TypeError("vessel_mask must be a 2-D boolean array.")
        return vessel & ~self.centered_circle_mask_for(
            vessel.shape,
            fallback_radius_pixels=fallback_radius_pixels,
        )

    def annulus_geometry(
        self,
        image_shape: tuple[int, int],
        number_of_radii_in_fov: int | None = None,
        *,
        annulus_size: float | None = None,
        annulus_stride: float | None = None,
        number_of_radii_in_half_diagonal: int = 32,
    ) -> AnnulusGeometry:
        """Derive full-width annuli with half-diagonal sampling density.

        The count includes intermediate annuli: stride = half-diagonal / count
        and size = 2 * stride. Starts extend to the farthest image corner from
        the optic disc, without clipping outer boundaries or dropping tails.
        Explicit native-pixel size/stride and legacy FOV-count calculation
        paths remain available with their existing terminal-boundary behavior.
        """

        for name, value in (("annulus_size", annulus_size), ("annulus_stride", annulus_stride)):
            if value is not None and (
                isinstance(value, bool) or not isinstance(value, Real)
                or not np.isfinite(value) or value <= 0
            ):
                raise ValueError(f"{name} must be finite and strictly positive (native pixels).")
        legacy = annulus_size is None and annulus_stride is None and number_of_radii_in_fov is not None
        count_based = annulus_size is None and annulus_stride is None and not legacy
        if legacy and (isinstance(number_of_radii_in_fov, bool)
                       or not isinstance(number_of_radii_in_fov, Integral) or number_of_radii_in_fov < 1):
            raise ValueError("number_of_radii_in_FOV must be positive and an integer.")
        if count_based and (
            isinstance(number_of_radii_in_half_diagonal, bool)
            or not isinstance(number_of_radii_in_half_diagonal, Integral)
            or number_of_radii_in_half_diagonal < 1
        ):
            raise ValueError("number_of_radii_in_half_diagonal must be positive and an integer.")
        self._required_dimensions("derive annulus geometry")
        ny, nx = _validated_image_shape(image_shape)
        radius_scale = max(image_half_diagonal(ny, nx), 1.0)
        if count_based:
            stride = radius_scale / number_of_radii_in_half_diagonal
            size = 2.0 * stride
            inner = self.centered_circle_radius_pixels() / radius_scale
            cx, cy = self.center
            farthest_radius = float(np.hypot(
                max(abs(cx), abs(nx - 1 - cx)),
                max(abs(cy), abs(ny - 1 - cy)),
            ))
            # The branch-support mask evaluates distances in float32. Round its
            # limit outward so a corner cannot disappear through roundoff.
            outer = float(np.nextafter(
                np.float32(farthest_radius / radius_scale), np.float32(np.inf),
            ))
            count = annulus_count(
                farthest_radius - inner * radius_scale, size, stride, retain_all=True,
            )
            return AnnulusGeometry(
                inner_radius_frac=inner,
                outer_radius_frac=outer,
                ring_width_frac=size / radius_scale,
                ring_count=count,
                segment_length_frac=size / radius_scale,
                ring_stride_frac=stride / radius_scale,
                clip_outer_radius=False,
            )
        size = float(annulus_size) if annulus_size is not None else max(nx, ny) / float(
            number_of_radii_in_fov if legacy else 25
        )
        stride = float(annulus_stride) if annulus_stride is not None else size / (1.0 if legacy else 2.0)
        radial_width = size / radius_scale
        inner = min(self.centered_circle_radius_pixels() / radius_scale, 1.0)
        outer = 1.0
        count = annulus_count((outer - inner) * radius_scale, size, stride)
        return AnnulusGeometry(
            inner_radius_frac=inner,
            outer_radius_frac=outer,
            ring_width_frac=radial_width,
            ring_count=count,
            segment_length_frac=radial_width,
            ring_stride_frac=stride / radius_scale,
        )

    def _required_dimensions(self, purpose: str) -> tuple[float, float]:
        if self.width is None or self.height is None:
            raise ValueError(f"optic-disc width and height are required to {purpose}.")
        return self.width, self.height


def _optional_positive_scalar(value, name: str) -> float | None:
    if value is None:
        return None
    values = np.asarray(value, dtype=np.float64).reshape(-1)
    if values.size != 1 or not np.isfinite(values[0]) or values[0] <= 0:
        raise ValueError(f"optic-disc {name} must be a finite positive scalar.")
    return float(values[0])


def _valid_positive_scalar_or_none(value) -> float | None:
    if value is None:
        return None
    values = np.asarray(value, dtype=np.float64).reshape(-1)
    if values.size != 1 or not np.isfinite(values[0]) or values[0] <= 0:
        return None
    return float(values[0])


def _valid_center_or_none(value) -> tuple[float, float] | None:
    if value is None:
        return None
    values = np.asarray(value, dtype=np.float64).reshape(-1)
    if values.size != 2 or not np.all(np.isfinite(values)):
        return None
    return float(values[0]), float(values[1])


def _resolved_source_mask(
    value,
    image_shape: tuple[int, int],
) -> tuple[np.ndarray | None, bool]:
    if value is None:
        return None, False
    mask = np.asarray(value)
    if mask.ndim != 2:
        return None, True
    if mask.shape != image_shape:
        return None, True
    try:
        if not np.all(np.isfinite(mask)):
            return None, True
    except TypeError:
        return None, True
    mask = np.asarray(mask, dtype=bool)
    if not np.any(mask):
        return None, False
    return mask, False


def _optional_nonnegative_scalar(value, name: str) -> float | None:
    if value is None:
        return None
    values = np.asarray(value, dtype=np.float64).reshape(-1)
    if values.size != 1 or not np.isfinite(values[0]) or values[0] < 0:
        raise ValueError(f"optic-disc {name} must be a finite nonnegative scalar.")
    return float(values[0])
