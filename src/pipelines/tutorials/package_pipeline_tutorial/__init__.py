"""Minimal tutorial for a package-based pipeline.

Keeping the registration here makes the runnable entrypoint easy to find, while
the actual implementation can be split across helper modules next to this file.
"""

from pipeline_engine.imports import pipeline

from .runner import run_package_pipeline_tutorial


@pipeline(
    name="package_pipeline_tutorial",
    description="Tutorial: minimal package-based pipeline skeleton.",
    input_slot="both",
    visibility="hidden",
)
def run(ctx) -> None:
    """Entrypoint called by the pipeline runtime."""

    run_package_pipeline_tutorial(ctx)


__all__ = ["run"]
