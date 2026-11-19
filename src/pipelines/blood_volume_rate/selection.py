"""Effective blood-volume-rate products required by a run."""

from __future__ import annotations


def selected_blood_volume_rate_methods(
    selection=None, *, options=(),
) -> frozenset[str]:
    """Resolve explicit and report-required blood-volume-rate products."""

    selected = set(options)
    if selection is not None:
        selected.update(selection.options_for("blood_volume_rate"))
    if selection is not None and selection.pipeline_scheduled("pdf_report"):
        selected.add("masked_edges")
    return frozenset(selected)


__all__ = ["selected_blood_volume_rate_methods"]
