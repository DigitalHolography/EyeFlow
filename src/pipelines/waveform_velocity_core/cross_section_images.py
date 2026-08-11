"""Export per-segment cross-section images."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from calculations.blood_flow_velocity import CrossSectionSignalResult


def export_rotated_mean_pngs(
    output,
    result: CrossSectionSignalResult,
    vessel_folder: str,
) -> list[Path]:
    """Export every valid unmasked and masked rotated time-mean image."""
    sample_counts = np.asarray(result.profile_sample_count, dtype=np.int32)
    branch_ids = np.asarray(result.branch_ids, dtype=np.int32)
    paths: list[Path] = []
    variants = (
        ("rotated_mean", result.rotated_mean_images),
        ("rotated_mean_masked", result.rotated_mean_images_masked),
    )
    for root_folder, variant_images in variants:
        images = np.asarray(variant_images, dtype=np.float32)
        for ring_index in range(images.shape[0]):
            for branch_index, branch_id in enumerate(branch_ids):
                if sample_counts[ring_index, branch_index] <= 0:
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
