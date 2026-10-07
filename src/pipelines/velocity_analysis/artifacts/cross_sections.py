"""Export per-segment cross-section images."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..models import VelocitySegmentResult


def export_rotated_mean_pngs(
    output,
    result: VelocitySegmentResult,
    vessel_folder: str,
) -> list[Path]:
    """Export every valid unmasked and masked rotated time-mean image."""
    profile = result.profile
    valid_segments = np.asarray(profile.topology.valid_segments, dtype=bool)
    branch_ids = np.asarray(profile.topology.native.branch_ids, dtype=np.int32)
    paths: list[Path] = []
    variants = (
        ("rotated_mean", profile.mean_images.unmasked),
        ("rotated_mean_masked", profile.mean_images.masked),
    )
    for root_folder, variant_images in variants:
        images = np.asarray(variant_images, dtype=np.float32)
        for ring_index in range(images.shape[0]):
            for branch_index, branch_id in enumerate(branch_ids):
                if not valid_segments[ring_index, branch_index]:
                    continue
                path = output.write_png(
                    images[ring_index, branch_index],
                    (
                        f"{root_folder}/{vessel_folder}/"
                        f"ring_{ring_index + 1:03d}_"
                        f"branch_{int(branch_id):03d}.png"
                    ),
                )
                paths.append(Path(path))
    return paths
