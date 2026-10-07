"""Fixed-source input readers shared by EyeFlow pipelines."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from input_output.schema.source_data import HolodopplerTiming

if TYPE_CHECKING:
    from pipeline_engine import PipelineContext


def resolve_holodoppler_timing(
    pipeline_input: PipelineContext,
) -> HolodopplerTiming:
    from input_output.schema import HolodopplerSource

    return HolodopplerSource.from_context(pipeline_input).timing()


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
