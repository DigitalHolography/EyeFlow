"""Device-resident resize, masking, rotation and reduction for one segment.

Only small images used for orientation/fitting cross to the CPU mid-pipeline.
The masked movie is reduced on-device and is never downloaded.
"""

from __future__ import annotations


def measure_cross_section_gpu(
    backend,
    sub_stack,
    sub_mask,
    loc_xy,
    optic_disc_center,
    tilt_angle_mask,
    settings,
    side_pixels,
    angle_override=None,
    limits_override=None,
):
    # Local import avoids a cycle with the public dispatcher.
    from .generate_cross_section_signals import (
        _INTERPOLATED_SUBSTACK_SIDE,
        _ROTATED_SUBSTACK_SIDE,
        _cross_section_limits,
        _CrossSectionMeasurement,
        _CrossSectionVelocityMeasurement,
        _interpolated_pixel_size_mm,
        _mean_image_rotation_angle,
        _rotated_profile_sample_count,
        _sample_nanstd_axis0,
    )

    xp, ndi = backend.cupy, backend.ndimage

    def mean(values, axis):
        finite = xp.isfinite(values)
        total = xp.sum(xp.where(finite, values, xp.float32(0)), axis=axis, dtype=xp.float32)
        count = xp.sum(finite, axis=axis)
        result = xp.full(total.shape, xp.nan, dtype=xp.float32)
        xp.divide(total, xp.where(count > 0, count, 1), out=result)
        result[count == 0] = xp.nan
        return result

    def normalized_transform(values, transform, threshold):
        finite = xp.isfinite(values)
        numerator = transform(xp.where(finite, values, xp.float32(0)))
        weights = transform(finite.astype(xp.float32))
        result = xp.full(numerator.shape, xp.nan, dtype=xp.float32)
        keep = weights > threshold if threshold < 0.5 else weights >= threshold
        xp.divide(numerator, xp.where(keep, weights, xp.float32(1)), out=result)
        result[~keep] = xp.nan
        return result

    side = _INTERPOLATED_SUBSTACK_SIDE
    stack = xp.asarray(sub_stack, dtype=xp.float32)
    zoom = (1.0, side / stack.shape[-2], side / stack.shape[-1])
    resized = normalized_transform(
        stack,
        lambda a: ndi.zoom(
            a,
            zoom,
            order=1,
            mode="grid-constant",
            cval=0.0,
            prefilter=False,
            grid_mode=True,
        ),
        1e-6,
    )
    del stack
    mask = ndi.zoom(
        xp.asarray(sub_mask, dtype=xp.float32),
        zoom[-2:],
        order=0,
        mode="grid-constant",
        cval=0.0,
        prefilter=False,
        grid_mode=True,
    ) >= xp.float32(0.5)
    mean_image = mean(resized, axis=0)
    mean_masked = xp.where(mask, mean_image, xp.nan)
    if angle_override is None:
        angle = _mean_image_rotation_angle(
            xp.asnumpy(mean_masked),
            loc_xy,
            optic_disc_center,
            tilt_angle_mask,
            settings,
        )
    else:
        angle = float(angle_override)

    padding = _ROTATED_SUBSTACK_SIDE - side
    before, after = padding // 2, padding - padding // 2

    def pad(values, fill):
        return xp.pad(
            values,
            [(0, 0)] * (values.ndim - 2) + [(before, after)] * 2,
            mode="constant",
            constant_values=fill,
        )

    def rotate(values):
        return ndi.rotate(
            values,
            angle,
            axes=(-2, -1),
            reshape=False,
            order=1,
            mode="constant",
            cval=0.0,
            prefilter=False,
        )

    rotated_mask = rotate(pad(mask.astype(xp.float32), 0)) >= xp.float32(0.5)
    rotated_mean = normalized_transform(pad(mean_image, xp.nan), rotate, 0.5)
    rotated_mean_masked = normalized_transform(pad(mean_masked, xp.nan), rotate, 0.5)
    rotated_mean_masked[~rotated_mask] = xp.nan
    host_mean = xp.asnumpy(rotated_mean)
    host_mean_masked = xp.asnumpy(rotated_mean_masked)
    limits = limits_override
    if limits is None:
        limits = _cross_section_limits(
            host_mean_masked,
            settings,
            pixel_size_mm=_interpolated_pixel_size_mm(settings.pixel_size_mm, side_pixels),
        )
    c1, c2 = limits

    def profiles(values, host_image, retain_movie):
        rotated = normalized_transform(pad(values, xp.nan), rotate, 0.5)
        transverse = mean(rotated, axis=1)
        longitudinal = mean(rotated, axis=2)
        return _CrossSectionVelocityMeasurement(
            raw=xp.asnumpy(mean(transverse[:, c1 : c2 + 1], axis=1)),
            safe_velocity=xp.asnumpy(mean(transverse, axis=1)),
            transverse_profiles=xp.asnumpy(transverse),
            longitudinal_profiles=xp.asnumpy(longitudinal),
            rotated_stack=xp.asnumpy(rotated) if retain_movie else None,
            angle=angle,
            spatial_std=_sample_nanstd_axis0(host_image),
        )

    unmasked = profiles(resized, host_mean, True)
    resized[:, ~mask] = xp.nan
    masked = profiles(resized, host_mean_masked, False)
    return _CrossSectionMeasurement(
        unmasked=unmasked,
        masked=masked,
        rotated_mean=host_mean,
        rotated_mean_masked=host_mean_masked,
        rotated_mask=xp.asnumpy(rotated_mask),
        limits=(c1, c2),
        sample_count=_rotated_profile_sample_count(angle),
    )
