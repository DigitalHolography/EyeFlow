"""HDF5 record encoding and atomic publication of completed local files."""

import os
import shutil
from pathlib import Path
from uuid import uuid4
from time import perf_counter

import h5py

from .h5_storage import _create_dataset
from utils.logger import Logger


def _write_attributes(node, record):
    """Overlay record attributes without discarding native HDF5 payload attrs."""
    for key, value in record.attrs.items():
        if key in node.attrs:
            del node.attrs[key]
        dtype = record.attr_dtypes.get(key)
        if dtype is None:
            node.attrs[key] = value
        else:
            node.attrs.create(key, value, dtype=dtype)


def reconcile_attributes(node, record):
    """Replace authoritative mutable group/root metadata, including removals."""
    for key in tuple(node.attrs):
        if key not in record.attrs:
            del node.attrs[key]
    _write_attributes(node, record)


def write_record(output, path, record, *, metrics=None):
    """Use the same storage policy for progressive output and export snapshots."""
    from pipeline_engine.data.storage import ArrayReference, ValueReference, H5ArrayReference
    import numpy as np

    if record.kind == "group":
        node = output.require_group(path) if path else output
    else:
        parent, _, name = path.rpartition("/")
        group = output.require_group(parent) if parent else output
        value = record.value
        started = perf_counter()
        payload = value.open() if isinstance(value, (ArrayReference, ValueReference, H5ArrayReference)) else value
        prepared = perf_counter()
        try:
            node = _create_dataset(group, name, payload, h5_options=record.options)
        finally:
            if isinstance(payload, np.memmap):
                payload._mmap.close()
        if metrics is not None:
            metrics['output_reference_prepare_seconds'] = metrics.get('output_reference_prepare_seconds',0)+prepared-started
            metrics['output_hdf5_write_seconds'] = metrics.get('output_hdf5_write_seconds',0)+perf_counter()-prepared
    started = perf_counter()
    _write_attributes(node, record)
    if metrics is not None:
        metrics['output_dataset_attribute_seconds'] = metrics.get('output_dataset_attribute_seconds',0)+perf_counter()-started
    return node


def copy_output_node(model, source, target, name=None):
    """Copy a managed subtree into an isolated HDF5 export snapshot."""
    name = name or source.name.rsplit("/", 1)[-1]
    record = model.records[source.path]
    if record.kind == "group":
        node = target.create_group(name)
        for child_name, child in source.items():
            copy_output_node(model, child, node, child_name)
    else:
        node = _create_dataset(
            target, name, model.storage.read(record.value), h5_options=record.options
        )
    _write_attributes(node, record)


class OutputPublisher:
    """Publish a closed complete file without exposing an incomplete filename."""

    def __init__(self, destination, metrics=None, *, cancelled=None):
        self.destination, self.metrics = Path(destination), metrics
        self.cancelled = cancelled

    def publish(self, local):
        started = perf_counter()
        final = self.destination
        incoming = final.with_name(f".ef-{uuid4().hex[:16]}.h5")
        try:
            if not h5py.is_hdf5(local):
                raise OSError("Completed local output is not an HDF5 file")
            final.parent.mkdir(parents=True, exist_ok=True)
            Logger.log_perf(f"START pipeline=run variant=shared stage=hdf5_publication at={started:.6f}")
            shutil.copyfile(local, incoming)
            copied = perf_counter()
            if self.metrics:
                self.metrics.add("h5_publication_copy_seconds", copied-started)
            if self.cancelled is not None and self.cancelled():
                raise OSError("Output publication was cancelled")
            os.replace(incoming, final)
            if self.metrics:
                self.metrics.add("h5_publication_replace_seconds", perf_counter()-copied)
            return final
        finally:
            incoming.unlink(missing_ok=True)
            if self.metrics:
                self.metrics.add("final_publication_seconds", perf_counter()-started)
            Logger.log_perf(f"END pipeline=run variant=shared stage=hdf5_publication at={perf_counter():.6f} wall={perf_counter()-started:.3f}s")
