"""Write Encapsulated PostScript figure artifacts for EyeFlow runs."""

from __future__ import annotations

from pathlib import Path

from .artifact_names import acquisition_stem, prefixed_artifact_path, writer_artifact_path
from .deferred import save_figure_file
from .png import _close_figure


class FigureArtifactWriter:
    """Write stem-prefixed EPS figures for one output namespace."""

    def __init__(self, output, stem: str | None = None) -> None:
        self.output = output
        self.stem = str(stem) if stem else acquisition_stem(output)

    def path(self, suffix: str) -> Path:
        return writer_artifact_path(self.output, _eps_output_type(), suffix, label=self.stem)

    def save_figure(self, fig, suffix: str, *, dpi: int = 150) -> Path:
        return write_eps_file(self.path(suffix), fig, dpi=dpi, close=True)


def write_eps_file(
    path: str | Path, fig, *, dpi: int = 150, stem: str | None = None,
    close: bool = False,
) -> Path:
    """Save a Matplotlib figure as EPS, creating parent directories."""

    target = prefixed_artifact_path(path, stem) if stem is not None else Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        save_figure_file(target, fig, format="eps", dpi=dpi)
    finally:
        if close:
            _close_figure(fig)
    return target


# Format-specific name for callers that do not use the shared writer name.
EpsArtifactWriter = FigureArtifactWriter


def _eps_output_type():
    from input_output.output_manager import OutputType

    return OutputType.EPS


__all__ = ["EpsArtifactWriter", "write_eps_file"]
