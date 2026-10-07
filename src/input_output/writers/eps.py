"""Write Encapsulated PostScript figure artifacts for EyeFlow runs."""

from __future__ import annotations

from pathlib import Path

from .artifact_names import acquisition_stem, labeled_artifact_path, prefixed_artifact_path


class FigureArtifactWriter:
    """Write stem-prefixed EPS figures for one output namespace."""

    def __init__(self, output, stem: str | None = None) -> None:
        self.output = output
        self.stem = str(stem) if stem else _output_stem(output)
        self.acquisition_stem = acquisition_stem(output, Path(self.stem).name)

    def path(self, suffix: str) -> Path:
        filename = str(labeled_artifact_path(suffix, self.acquisition_stem, self.stem))
        path = self.output.path_for(_eps_output_type(), filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def save_figure(self, fig, suffix: str, *, dpi: int = 150) -> Path:
        path = write_eps_file(self.path(suffix), fig, dpi=dpi)
        _close_figure(fig)
        return path


def write_eps_file(path: str | Path, fig, *, dpi: int = 150, stem: str | None = None) -> Path:
    """Save a Matplotlib figure as EPS, creating parent directories."""

    target = prefixed_artifact_path(path, stem) if stem is not None else Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, format="eps", dpi=dpi)
    return target


# Format-specific name for callers that do not use the shared writer name.
EpsArtifactWriter = FigureArtifactWriter


def _output_stem(output) -> str:
    return acquisition_stem(output)


def _close_figure(fig) -> None:
    import matplotlib.pyplot as plt

    plt.close(fig)


def _eps_output_type():
    from input_output.output_manager import OutputType

    return OutputType.EPS


__all__ = ["EpsArtifactWriter", "write_eps_file"]
