"""Effective blood-volume-rate products required by a run."""

from __future__ import annotations


def selected_blood_volume_rate_methods(ctx) -> frozenset[str]:
    """Include mask geometry whenever the mandatory report page is scheduled."""

    selected = set(ctx.options_for("blood_volume_rate"))
    if ctx.pipeline_scheduled("pdf_report"):
        selected.add("masked_edges")
    return frozenset(selected)


__all__ = ["selected_blood_volume_rate_methods"]
