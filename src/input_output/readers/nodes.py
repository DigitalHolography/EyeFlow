"""Shared dataset/group predicates for HDF5 and managed runtime adapters."""
import h5py


def is_dataset(node):
    return isinstance(node, h5py.Dataset) or getattr(node, "_eyeflow_dataset", False) is True


def is_group(node):
    return isinstance(node, h5py.Group) or getattr(node, "_eyeflow_group", False) is True
