"""Minimal tutorial for a package-based pipeline.

The runtime imports this package from the static pipeline catalog. Keeping the
registration here makes the runnable entrypoint easy to find, while the actual
implementation can be split across helper modules next to this file.
"""

from pipelines.imports import pipeline

from .runner import run_package_pipeline_tutorial


@pipeline(
    name="package_pipeline_tutorial",
    description="Tutorial: minimal package-based pipeline skeleton.",
    input_slot="both",
)
def run(ctx) -> None:
    """Entrypoint called by the pipeline runtime."""

    run_package_pipeline_tutorial(ctx)


__all__ = ["run"]
