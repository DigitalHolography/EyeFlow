"""Compatibility and immutable storage for run-local continuous gradients."""

from __future__ import annotations

from hashlib import blake2b

import h5py
from input_output.readers.nodes import is_dataset
import numpy as np

from calculations.vessel_segments.measurement import SegmentMeasurementSettings
from calculations.vessel_segments.sampling.models import SegmentSamplingPlan


def continuous_gradient_key(moment0ff, masks, topologies, settings, filter_settings):
    """Identify read-only HD data and all spatial inputs, excluding beat timing.

    Dataset identity is scoped to the open input file. Fingerprint spatial arrays
    so equivalent plans can share results and changed masks/rotations cannot.
    Unsupported direct-call inputs bypass sharing rather than guess identity.
    """
    if not is_dataset(moment0ff) or not all(
        isinstance(plan, SegmentSamplingPlan) for plan in topologies.values()
    ):
        return None
    plans = []
    for name in ("artery", "vein"):
        plan = topologies[name]
        native = plan.native
        plans.append((
            name,
            native.ring_settings,
            native.optic_disc_center_xy,
            _fingerprint(
                masks[name], native.labels, native.branch_ids, native.centerline,
                native.annulus_masks, native.segment_masks,
                native.segment_centers_xy, native.window_bounds_xyxy,
                native.optic_disc_mask, native.delta_radius,
                plan.rotation_degrees, plan.interpolated_masks, plan.rotated_masks,
            ),
        ))
    return (
        "spatial_gradient_moment0.continuous",
        moment0ff.id, tuple(moment0ff.shape), moment0ff.dtype.str,
        tuple(plans), SegmentMeasurementSettings.from_value(settings), filter_settings,
    )


def freeze_continuous_gradients(segments):
    """Protect cached measurement arrays from downstream in-place changes."""
    for result in segments:
        result.segment_signal.setflags(write=False)
        for pair in (result.transverse, result.longitudinal, result.mean_images):
            if pair is not None:
                pair.unmasked.setflags(write=False)
                pair.masked.setflags(write=False)
    return segments


def _fingerprint(*arrays):
    digest = blake2b(digest_size=32)
    for array in arrays:
        if array is None:
            digest.update(b"None")
            continue
        values = np.ascontiguousarray(array)
        digest.update(repr((values.shape, values.dtype.str)).encode("ascii"))
        digest.update(values.view(np.uint8))
    return digest.digest()
