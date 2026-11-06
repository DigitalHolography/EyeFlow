"""Select compact displacement and spatial-gradient report measurements."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from input_output.readers.nodes import is_dataset, is_group
from input_output.schema import EyeFlowOutputPaths
from input_output.schema.eyeflow_output import processing_path
from input_output.schema.provenance import read_dataset_attrs


@dataclass(frozen=True)
class PhaseSummary:
    label: str
    curve: np.ndarray
    median: float | None
    valid_fraction: float
    branch_count: int
    beat_count: int | None
    unit: str


@dataclass(frozen=True)
class WidthReportData:
    displacement: tuple[PhaseSummary, ...]
    gradient: tuple[PhaseSummary, ...]


def extract_width_report_data(h5file) -> WidthReportData | None:
    """Read only persisted phase summaries; leave absent analyses absent."""
    paths = EyeFlowOutputPaths.active()
    method_group = h5file.get(paths.displacement_profiles_root)
    methods = sorted(method_group.keys()) if is_group(method_group) else []
    displacement = []
    for candidate in methods:
        for vessel in ("Artery", "Vein"):
            dataset = h5file.get(
                f"{paths.displacement_profiles_root}/{candidate}/{vessel}/"
                "Cross_sectional_radial_movement_amplitude_profile"
            )
            series = _phase_summary(dataset, f"{candidate} {vessel}", rank=4)
            if series is not None:
                displacement.append(series)

    gradient = []
    for root, label in (("ProcessingMoments", "Moments"),
                        ("ProcessingBandRatio", "Band ratio")):
        metrics_root = processing_path(paths.spatial_gradient_metrics_root, root)
        for vessel in ("Artery", "Vein"):
            dataset = h5file.get(
                f"{metrics_root}/{vessel}/Transverse/Masked/tk/lumen_size"
            )
            series = _phase_summary(dataset, f"{label} {vessel}", rank=2)
            if series is not None:
                gradient.append(series)

    if not displacement and not gradient:
        return None
    return WidthReportData(tuple(displacement), tuple(gradient))


def _phase_summary(dataset, label: str, *, rank: int) -> PhaseSummary | None:
    """Reduce one bounded phase slab at a time to its finite-sample median."""
    if not is_dataset(dataset) or dataset.ndim != rank or not all(dataset.shape):
        return None
    curve = np.full(dataset.shape[0], np.nan, dtype=float)
    finite_count = 0
    for phase in range(dataset.shape[0]):
        slab = np.asarray(dataset[phase], dtype=float)
        finite = slab[np.isfinite(slab)]
        if finite.size:
            curve[phase] = float(np.median(finite))
            finite_count += finite.size
    valid_phases = curve[np.isfinite(curve)]
    unit = read_dataset_attrs(dataset).get("unit", "pixels")
    if isinstance(unit, bytes):
        unit = unit.decode("utf-8")
    return PhaseSummary(
        label=label,
        curve=curve,
        median=float(np.median(valid_phases)) if valid_phases.size else None,
        valid_fraction=finite_count / dataset.size,
        branch_count=dataset.shape[2] if rank == 4 else dataset.shape[1],
        beat_count=dataset.shape[1] if rank == 4 else None,
        unit=str(unit),
    )
