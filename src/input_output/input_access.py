"""Fixed-source input readers shared by EyeFlow pipelines."""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import TYPE_CHECKING

import h5py
import numpy as np

from input_output.schema.base import HolodopplerTiming

if TYPE_CHECKING:
    from collections.abc import Mapping

    from pipeline_engine import PipelineContext


@dataclass(frozen=True)
class ResolvedArray:
    path: str
    value: np.ndarray


def resolve_required_source_array(
    source: h5py.File | None,
    *,
    source_name: str,
    logical_name: str,
    path: str,
) -> ResolvedArray:
    if source is None:
        raise KeyError(f"{source_name} HDF5 input is required for '{logical_name}'.")
    found = source.get(path)
    if not isinstance(found, h5py.Dataset):
        raise KeyError(
            f"Missing {source_name} dataset for '{logical_name}'. Expected: {path}"
        )
    return ResolvedArray(path=path, value=np.asarray(found[()]))


def resolve_holodoppler_timing(
    pipeline_input: PipelineContext,
) -> HolodopplerTiming:
    from input_output.schema import HolodopplerSource

    return HolodopplerSource.from_context(pipeline_input).timing()


def resolve_dt_seconds(pipeline_input: PipelineContext) -> float:
    return resolve_holodoppler_timing(pipeline_input).dt_seconds


def read_first_attr(pipeline_input: PipelineContext, *keys: str):
    for key in keys:
        value = pipeline_input.attrs.get(key, None)
        scalar = _scalar_from_value(value)
        if scalar is not None:
            return scalar
    return None


def read_int_setting(
    pipeline_input: PipelineContext,
    *,
    default: int,
    keys: tuple[str, ...],
) -> int:
    value = read_first_attr(pipeline_input, *keys)
    if value is None:
        return int(default)
    return int(value)


def read_nested_int_setting(
    config: Mapping[str, object],
    section: str,
    key: str,
    *,
    default: int,
) -> int:
    section_value = config.get(section, {})
    if not isinstance(section_value, dict):
        return int(default)
    value = _scalar_from_value(section_value.get(key))
    return int(default) if value is None else int(value)


def _scalar_from_value(value):
    if value is None:
        return None
    array = np.asarray(value).reshape(-1)
    if array.size == 0:
        return None
    scalar = array[0]
    if isinstance(scalar, bytes):
        return scalar.decode("utf-8")
    return scalar.item() if hasattr(scalar, "item") else scalar
