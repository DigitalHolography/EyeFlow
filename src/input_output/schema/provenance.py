"""Resolve workflow provenance without duplicating it on every dataset."""

from __future__ import annotations

import h5py


INHERITED_PROVENANCE_KEYS_ATTR = "inherited_provenance_keys"


def inherited_provenance(h5obj: h5py.Group | h5py.Dataset) -> dict[str, object]:
    """Read only explicitly declared provenance, with nearest groups winning.

    Other group attributes (units, dimensions, metric definitions, etc.) do not
    implicitly apply to children. File-root attributes are not inherited across
    the separate workflow namespaces.
    """
    groups = []
    group = h5obj if isinstance(h5obj, h5py.Group) else h5obj.parent
    while group.name != "/":
        groups.append(group)
        group = group.parent
    result = {}
    for group in reversed(groups):
        for key in group.attrs.get(INHERITED_PROVENANCE_KEYS_ATTR, ()):
            result[key] = group.attrs[key]
    return result


def read_dataset_attrs(dataset: h5py.Dataset) -> dict[str, object]:
    """Return effective attributes for new and legacy EyeFlow datasets.

    Dataset-local values override inherited provenance. Legacy datasets without
    an inheritance declaration retain their original, local-only behavior.
    """
    return {**inherited_provenance(dataset), **dict(dataset.attrs)}
