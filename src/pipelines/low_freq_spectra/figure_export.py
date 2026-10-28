"""Snapshot reduced spectral products for an isolated background renderer."""

from pathlib import Path

import h5py

from input_output.writers.h5_output import PipelineH5Output
from input_output.schema import EyeFlowOutputPaths
from pipeline_engine.context import PipelineOutput

from .artifacts import export_spectrogram_figures


def snapshot_spectral_figures(output, path: Path) -> None:
    """Copy only reduced spectra and the published mask, on the HDF5 owner."""
    schema = EyeFlowOutputPaths.active()
    spectral_root = schema.spectrograms_frequency_bin_table.split("/", 1)[0]
    with h5py.File(path, "w") as snapshot:
        for name in (spectral_root, schema.segmentation.gaps):
            source = output.h5.get(name)
            if source is None:
                continue
            parent, _, leaf = name.rpartition("/")
            target = snapshot.require_group(parent) if parent else snapshot
            source.file.copy(source, target, name=leaf)


def export_spectral_snapshot(path, manager, *, dt_seconds=None):
    """Child process entry point; no source or result HDF5 handles are shared."""
    with h5py.File(path, "r") as snapshot:
        output = PipelineOutput(manager, PipelineH5Output(snapshot))
        return export_spectrogram_figures(output, dt_seconds=dt_seconds)
