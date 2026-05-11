from pipeline_engine.imports import np, pipeline, resolve_holodoppler_timing, with_attrs


HD_EXAMPLE_DATASET_KEY = "moment0"
DV_EXAMPLE_DATASET_KEY = "retinal_artery_mask"
DV_EXAMPLE_CONFIG_KEY = "local_background_dist"


def _numeric_dataset_summary(values) -> tuple[np.float32, np.int32]:
    if values.size == 0:
        return np.float32(np.nan), np.int32(0)

    try:
        numeric = np.asarray(values, dtype=np.float32)
    except (TypeError, ValueError):
        return np.float32(np.nan), np.int32(values.size)

    return np.float32(np.nanmean(numeric)), np.int32(numeric.size)


@pipeline(
    name="dual_input_tutorial",
    description="Tutorial: read HD and DV inputs simultaneously in one pipeline.",
    dag_produces=["dual_input_summary"],
    input_slot="both",
)
def run(ctx) -> None:
    """
    Minimal example for consuming both input files.

    Common operations:
    - `ctx.require_inputs("hd", "dv")` fails early when an input is missing.
    - `ctx.hd.array("moment0")` reads an HD dataset.
    - `ctx.dv.array("retinal_artery_mask")` reads a DV dataset.
    - `ctx.dv.config("local_background_dist")` reads a DV sidecar/H5 setting.
    - `ctx.vars["name"] = value` shares an in-memory value with later pipelines.
    """

    ctx.require_inputs("hd", "dv")

    timing = resolve_holodoppler_timing(ctx)
    hd_values = ctx.hd.array(HD_EXAMPLE_DATASET_KEY)
    dv_values = ctx.dv.array(DV_EXAMPLE_DATASET_KEY)
    local_background_dist = np.int32(ctx.dv.config(DV_EXAMPLE_CONFIG_KEY))

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

    ctx.set_var(
        "dual_input_summary",
        {
            "hd_example_mean": float(hd_mean),
            "dv_example_mean": float(dv_mean),
            "shared_root_groups": shared_root_keys,
        },
    )

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
            "hd_example_key": HD_EXAMPLE_DATASET_KEY,
            "dv_example_key": DV_EXAMPLE_DATASET_KEY,
            "dv_config_key": DV_EXAMPLE_CONFIG_KEY,
            "hd_example_path": ctx.hd.path(HD_EXAMPLE_DATASET_KEY),
            "dv_example_path": ctx.dv.path(DV_EXAMPLE_DATASET_KEY),
            "shared_root_groups": shared_root_keys,
        }
    )
