"""Read the retained processing summary and dataset attributes."""

from __future__ import annotations

import h5py
from input_output.readers.nodes import is_group


INHERITED_PROVENANCE_KEYS_ATTR = "inherited_provenance_keys"
PROCESSING_SUMMARY_KEYS = frozenset({
    "velocity_estimation_method", "velocity_inpainting_radius", "velocity_unit",
})
_WORKFLOW_SOURCE_KEYS = {
    "ProcessingMoments": frozenset({"moment0_source_path", "moment2_source_path"}),
    "ProcessingBandRatio": frozenset({"band_lf_source_path", "band_hf_source_path"}),
}


def processing_summary_keys(root: str) -> frozenset[str]:
    """Keep source paths specific to the workflow that consumes them."""
    return PROCESSING_SUMMARY_KEYS | _WORKFLOW_SOURCE_KEYS.get(root.strip("/"), frozenset())


def inherited_provenance(h5obj: h5py.Group | h5py.Dataset) -> dict[str, object]:
    """Read retained processing summaries, with nearest groups winning.

    Other group attributes (units, dimensions, metric definitions, etc.) do not
    implicitly apply to children. File-root attributes are not inherited across
    the separate workflow namespaces. Processing summaries inherit only the
    fixed summary keys; older explicit declarations remain readable.
    """
    groups = []
    group = h5obj if is_group(h5obj) else h5obj.parent
    while group.name != "/":
        groups.append(group)
        group = group.parent
    result = {}
    for group in reversed(groups):
        if group.name.strip("/") in _WORKFLOW_SOURCE_KEYS:
            result.update({
                key: group.attrs[key]
                for key in processing_summary_keys(group.name)
                if key in group.attrs
            })
        for key in group.attrs.get(INHERITED_PROVENANCE_KEYS_ATTR, ()):
            result[key] = group.attrs[key]
    return result


def read_dataset_attrs(dataset: h5py.Dataset) -> dict[str, object]:
    """Return effective attributes for new and legacy EyeFlow datasets.

    Dataset-local values override the retained processing summary.
    """
    return {**inherited_provenance(dataset), **dict(dataset.attrs)}
