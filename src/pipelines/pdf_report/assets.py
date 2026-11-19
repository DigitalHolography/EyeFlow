"""Resolve existing raster artifacts for report presentations."""

from __future__ import annotations

from pathlib import Path

import numpy as np

from input_output.reports.pdf_report import load_first_image
from input_output.writers.artifact_names import prefixed_artifact_path
from input_output.writers.png import ReportImagePaths

ANALYSIS_ARTIFACTS = {
    "histogram": lambda vessel: Path(f"histogramVelocity{vessel}.png"),
    "branch_identity": lambda vessel: (
        Path("branch_identity")
        / f"{vessel}_13_substack_boxes_on_labels_with_rings.png"
    ),
    "lumen_diameter": lambda vessel: (
        Path("lumen_diameter") / f"{vessel}_lumen_diameter_distribution.png"
    ),
    "blood_volume_rate": lambda vessel: (
        Path("blood_volume_rate") / f"{vessel}_blood_volume_rate.png"
    ),
}
ANALYSIS_GRAPH_TITLES = {
    "histogram": "Velocity Distribution over Time",
    "branch_identity": "Branch Identity and Sampling Regions",
    "lumen_diameter": "Lumen-Diameter Distribution",
    "blood_volume_rate": "Mask-Derived Blood-Volume Rate",
}


def load_velocity_artifact(
    png_dir: Path | None,
    folder_name: str,
    kind: str,
    vessel: str,
    report_images: ReportImagePaths,
) -> np.ndarray | None:
    """Prefer a registered artifact, with the established filename fallback."""
    candidates = []
    registered = report_images.get((kind, vessel))
    if registered is not None:
        candidates.append(Path(registered))
    if png_dir is not None:
        suffixes = {"vessel_map": f"{vessel}_seg_map_bkg.png"}
        candidates.append(Path(png_dir) / f"{folder_name}_{suffixes[kind]}")
    return load_first_image(candidates)


def load_analysis_artifact(
    png_dir: Path | None,
    folder_name: str,
    kind: str,
    vessel: str,
) -> np.ndarray | None:
    """Load a page-two artifact from its canonical PNG subfolder."""
    if png_dir is None:
        return None
    suffix = ANALYSIS_ARTIFACTS[kind](vessel)
    path = Path(png_dir) / prefixed_artifact_path(suffix, folder_name)
    return load_first_image((path,))


__all__ = [
    "ANALYSIS_GRAPH_TITLES",
    "load_analysis_artifact",
    "load_velocity_artifact",
]
