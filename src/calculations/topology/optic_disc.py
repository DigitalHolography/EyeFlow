"""Authoritative optic-disc geometry from DopplerView."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .geometry import (
    AnnulusGeometry,
    _validated_center_xy,
    _validated_image_shape,
    image_half_diagonal,
)


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
    ) -> OpticDisc:
        """Resolve incomplete DopplerView measurements for one image frame.

        A finite, nonempty mask supplies any missing center or dimensions.  If
        the mask is invalid, or neither a usable mask nor complete geometry is
        available, a centered circular R0 with radius ten percent of the image
        half-diagonal is used.
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

        radius = 0.10 * image_half_diagonal(ny, nx)
        diameter = max(2.0 * radius, float(np.finfo(np.float32).eps))
        integer_radius = int(np.ceil(radius))
        center_x, center_y = nx / 2.0, ny / 2.0
        y, x = np.ogrid[:ny, :nx]
        fallback_mask = (y - center_y) ** 2 + (x - center_x) ** 2 <= float(integer_radius**2)
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
        number_of_radii_in_fov: int = 25,
    ) -> AnnulusGeometry:
        """Derive uniformly spaced annuli from this optic-disc geometry."""

        if number_of_radii_in_fov < 1:
            raise ValueError("number_of_radii_in_fov must be positive.")
        self._required_dimensions("derive annulus geometry")
        ny, nx = _validated_image_shape(image_shape)
        radius_scale = max(image_half_diagonal(ny, nx), 1.0)
        radial_step = max(nx, ny) / float(number_of_radii_in_fov) / radius_scale
        inner = min(self.centered_circle_radius_pixels() / radius_scale, 1.0)
        outer = 1.0
        count = max(1, int(np.ceil((outer - inner) / radial_step)))
        return AnnulusGeometry(
            inner_radius_frac=inner,
            outer_radius_frac=outer,
            ring_width_frac=radial_step,
            ring_count=count,
            segment_length_frac=radial_step,
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
