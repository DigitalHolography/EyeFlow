"""Write Encapsulated PostScript figure artifacts for EyeFlow runs."""

from __future__ import annotations

from pathlib import Path

from input_output.output_manager import OutputType, artifact_path, output_stem

from .png import _close_figure


class FigureArtifactWriter:
    """Write stem-prefixed EPS figures for one output namespace."""

    def __init__(self, output, stem: str | None = None) -> None:
        self.output = output
        self.stem = str(stem) if stem else output_stem(output)

    def path(self, suffix: str) -> Path:
        return artifact_path(self.output, OutputType.EPS, suffix, stem=self.stem)

    def save_figure(self, fig, suffix: str, *, dpi: int = 150) -> Path:
        return write_eps_file(self.path(suffix), fig, dpi=dpi, close=True)


def write_eps_file(path: str | Path, fig, *, dpi: int = 150, close: bool = False) -> Path:
    """Save a Matplotlib figure as EPS, creating parent directories."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    try:
        fig.savefig(target, format="eps", dpi=dpi)
    finally:
        if close:
            _close_figure(fig)
    return target


# Format-specific name for callers that do not use the shared writer name.
EpsArtifactWriter = FigureArtifactWriter


__all__ = ["EpsArtifactWriter", "write_eps_file"]

