from __future__ import annotations

import h5py
import numpy as np

from input_output import (
    DOPPLER_VIEW_SCHEMA,
    HOLODOPPLER_SCHEMA,
    H5DatasetSpec,
    H5SourceSchema,
    JsonConfigValueSpec,
)
from input_output.input_access import (
    read_nested_int_setting,
    resolve_holodoppler_timing,
    resolve_required_source_array,
)

from pipeline_engine import pipeline, with_attrs


HD_EXAMPLE_DATASET_KEY = "moment0"
DV_EXAMPLE_DATASET_KEY = "retinal_artery_mask"
DV_EXAMPLE_CONFIG_KEY = "local_background_dist"


def _read_schema_dataset(
    source: h5py.File,
    *,
    schema: H5SourceSchema,
    dataset_key: str,
) -> tuple[H5DatasetSpec, str, np.ndarray]:
    spec = schema.dataset(dataset_key)
    resolved = resolve_required_source_array(
        source,
        source_name=schema.label,
        logical_name=spec.key,
        path=spec.path,
    )
    return spec, resolved.path, resolved.value


def _numeric_dataset_summary(values: np.ndarray) -> tuple[np.float32, np.int32]:
    if values.size == 0:
        return np.float32(np.nan), np.int32(0)

    try:
        numeric = np.asarray(values, dtype=float)
    except (TypeError, ValueError):
        return np.float32(np.nan), np.int32(values.size)

    return np.float32(np.nanmean(numeric)), np.int32(numeric.size)


def _read_dv_int_config(ctx, *, key: str) -> tuple[JsonConfigValueSpec, np.int32]:
    spec = DOPPLER_VIEW_SCHEMA.config_value(key)
    if spec.section is None:
        raw_value = spec.read_json_config(ctx.dv_config)
        value = spec.default if raw_value is None else raw_value
        return spec, np.int32(value)

    value = read_nested_int_setting(
        ctx.dv_config,
        spec.section,
        spec.json_key,
        default=int(spec.default or 0),
    )
    return spec, np.int32(value)


@pipeline(
    name="dual_input_tutorial",
    description="Tutorial: read HD and DV inputs simultaneously in one pipeline.",
    dag_produces=["dual_input_summary"],
    input_slot="both",
)
def run(ctx) -> None:
    """
    Tutorial pipeline showing how to consume HD and DV inputs at the same time.

    The incoming `ctx` exposes schema-aware readers:
    - `ctx.hd`: Holodoppler input reader
    - `ctx.dv`: DopplerVision input reader
    - `ctx.ef`: current root-level EyeFlow outputs
    - `ctx.write(path, value, **attrs)`: write one output
    - `ctx.write_many(metrics)`: write a metrics dictionary
    """

    if ctx.hd.h5file is None or ctx.dv.h5file is None:
        raise ValueError("dual_input_tutorial requires both HD and DV inputs.")

    timing = resolve_holodoppler_timing(ctx)
    hd_spec, hd_path, hd_values = _read_schema_dataset(
        ctx.hd,
        schema=HOLODOPPLER_SCHEMA,
        dataset_key=HD_EXAMPLE_DATASET_KEY,
    )
    dv_spec, dv_path, dv_values = _read_schema_dataset(
        ctx.dv,
        schema=DOPPLER_VIEW_SCHEMA,
        dataset_key=DV_EXAMPLE_DATASET_KEY,
    )
    dv_config_spec, local_background_dist = _read_dv_int_config(
        ctx,
        key=DV_EXAMPLE_CONFIG_KEY,
    )

    hd_mean, hd_size = _numeric_dataset_summary(hd_values)
    dv_mean, dv_size = _numeric_dataset_summary(dv_values)
    mean_delta = (
        np.float32(hd_mean - dv_mean)
        if np.isfinite(hd_mean) and np.isfinite(dv_mean)
        else np.float32(np.nan)
    )

    hd_root_keys = sorted(str(key) for key in ctx.hd.keys())
    dv_root_keys = sorted(str(key) for key in ctx.dv.keys())
    shared_root_keys = sorted(set(hd_root_keys) & set(dv_root_keys))

    ctx.write_many(
        {
            "summary/hd_root_group_count": with_attrs(
                np.int32(len(hd_root_keys)),
                {"unit": ["count"]},
            ),
            "summary/dv_root_group_count": with_attrs(
                np.int32(len(dv_root_keys)),
                {"unit": ["count"]},
            ),
            "summary/shared_root_group_count": with_attrs(
                np.int32(len(shared_root_keys)),
                {"unit": ["count"]},
            ),
            "summary/hd_example_size": with_attrs(hd_size, {"unit": ["samples"]}),
            "summary/dv_example_size": with_attrs(dv_size, {"unit": ["samples"]}),
            "summary/hd_example_mean": with_attrs(hd_mean, {"unit": ["a.u."]}),
            "summary/dv_example_mean": with_attrs(dv_mean, {"unit": ["a.u."]}),
            "summary/hd_minus_dv_example_mean": with_attrs(
                mean_delta,
                {"unit": ["a.u."]},
            ),
            "summary/both_inputs_available": with_attrs(
                np.uint8(1),
                {"unit": ["bool"]},
            ),
            "schema/hd_sampling_freq_hz": with_attrs(
                np.float32(timing.sampling_freq),
                {"unit": ["Hz"]},
            ),
            "schema/hd_batch_stride": with_attrs(
                np.float32(timing.batch_stride),
                {"unit": ["frames"]},
            ),
            "schema/hd_dt_seconds": with_attrs(
                np.float32(timing.dt_seconds),
                {"unit": ["s"]},
            ),
            "schema/dv_local_background_dist": with_attrs(
                local_background_dist,
                {"unit": ["pixels"]},
            ),
        }
    )
    ctx.set_attrs(
        {
            "hd_source_file": str(ctx.hd.filename or ""),
            "dv_source_file": str(ctx.dv.filename or ""),
            "hd_schema_label": HOLODOPPLER_SCHEMA.label,
            "dv_schema_label": DOPPLER_VIEW_SCHEMA.label,
            "hd_example_key": hd_spec.key,
            "dv_example_key": dv_spec.key,
            "dv_config_key": dv_config_spec.key,
            "hd_example_path": hd_path,
            "dv_example_path": dv_path,
            "hd_example_dtype": hd_spec.dtype or "",
            "dv_example_dtype": dv_spec.dtype or "",
            "hd_example_dims": list(hd_spec.dims),
            "dv_example_dims": list(dv_spec.dims),
            "shared_root_groups": shared_root_keys,
        }
    )
