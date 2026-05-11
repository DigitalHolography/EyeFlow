from pipelines.imports import pipeline


@pipeline(
    name="test2",
    description="Tutorial: consume a ctx.vars value produced by test1.",
    dag_requires=["test1_moment0_summary"],
    input_slot="both",
)
def run(ctx) -> None:
    """
    This is the second half of the tiny DAG tutorial.

    The DAG dependency above means selecting `test2` will also run `test1`.
    The `ctx.vars` dictionary is shared by all pipelines in this one run.
    """

    ctx.require_inputs("hd", "dv")

    moment0_summary = ctx.get_var("test1_moment0_summary")
    if moment0_summary is None:
        raise RuntimeError("test2 expected test1 to store test1_moment0_summary.")

    local_background_dist = ctx.dv.config("local_background_dist", default=2)

    ctx.set_var(
        "test2_context_example",
        {
            "moment0_shape_from_test1": moment0_summary["shape"],
            "dv_local_background_dist": int(local_background_dist),
        },
    )

    # Still no output H5 dataset is created. These variables exist only while
    # this pipeline run is active.
