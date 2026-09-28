"""Write Encapsulated PostScript figure artifacts for EyeFlow runs."""

from __future__ import annotations

from pathlib import Path


class FigureArtifactWriter:
    """Write stem-prefixed EPS figures for one output namespace."""

    def __init__(self, output, stem: str | None = None) -> None:
        self.output = output
        self.stem = str(stem) if stem else _output_stem(output)

    def path(self, suffix: str) -> Path:
        filename = f"{self.stem}_{suffix}"
        path = self.output.path_for(_eps_output_type(), filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def save_figure(self, fig, suffix: str, *, dpi: int = 150) -> Path:
        path = write_eps_file(self.path(suffix), fig, dpi=dpi)
        _close_figure(fig)
        return path

    def savefig(self, fig, suffix: str, *, dpi: int = 150) -> Path:
        return self.save_figure(fig, suffix, dpi=dpi)


def write_eps_file(path: str | Path, fig, *, dpi: int = 150) -> Path:
    """Save a Matplotlib figure as EPS, creating parent directories."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(target, format="eps", dpi=dpi)
    return target


# Format-specific name for callers that do not use the shared writer name.
EpsArtifactWriter = FigureArtifactWriter


def _output_stem(output) -> str:
    manager = getattr(output, "manager", None)
    layout = getattr(manager, "layout", None)
    stem = getattr(layout, "stem", None)
    return str(stem or "eyeflow")


def _close_figure(fig) -> None:
    import matplotlib.pyplot as plt

    plt.close(fig)


def _eps_output_type():
    from input_output.output_manager import OutputType

    return OutputType.EPS


__all__ = ["EpsArtifactWriter", "FigureArtifactWriter", "write_eps_file"]

