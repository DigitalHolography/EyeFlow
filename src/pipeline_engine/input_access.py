"""Context-aware input conveniences for pipeline runners."""

from __future__ import annotations

from typing import TYPE_CHECKING

from input_output.schema.base import scalar_from_value
from input_output.schema.source_data import HolodopplerTiming

if TYPE_CHECKING:
    from .context import PipelineContext


def resolve_holodoppler_timing(ctx: PipelineContext) -> HolodopplerTiming:
    return ctx.inputs.hd.as_holodoppler().timing()


def read_first_attr(ctx: PipelineContext, *keys: str):
    for key in keys:
        scalar = scalar_from_value(ctx.attrs.get(key, None))
        if scalar is not None:
            return scalar
    return None


def read_int_setting(
    ctx: PipelineContext,
    *,
    default: int,
    keys: tuple[str, ...],
) -> int:
    value = read_first_attr(ctx, *keys)
    return int(default) if value is None else int(value)
