"""Tests for Encapsulated PostScript figure output."""

from __future__ import annotations

import tempfile
from pathlib import Path

from matplotlib.figure import Figure

from input_output.output_manager import OutputManager
from input_output.writers.eps import EpsArtifactWriter


def test_eps_artifact_writer_creates_stem_prefixed_figure() -> None:
    with tempfile.TemporaryDirectory() as temp_dir:
        output = OutputManager.from_holo(
            Path(temp_dir) / "sample.holo",
            output_root=Path(temp_dir),
        )
        fig = Figure()
        fig.subplots().plot([0.0, 1.0], [1.0, 0.0])

        path = EpsArtifactWriter(output, "sample").save_figure(
            fig,
            "diagnostic.eps",
        )

        assert path == output.layout.ef_dir / "eps" / "sample_diagnostic.eps"
        assert path.read_bytes().startswith(b"%!PS-Adobe")
