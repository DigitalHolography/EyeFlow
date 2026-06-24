from pipeline_engine.imports import np, pipeline


@pipeline(
    name="test1",
    description="Tutorial: read HD moment0 and keep a value only in ctx.state.",
    dag_produces=["test1_moment0_summary"],
    input_slot="hd",
)
def run(ctx) -> None:
    """
    This is the first half of a tiny DAG tutorial.

    Select `test2` as a target and the DAG will run this pipeline first because
    `test2` requires the `test1_moment0_summary` key produced here.
    """

    ctx.require_inputs("hd")

    moment0 = ctx.inputs.hd.array("moment0", dtype=np.float32)

    ctx.state.set(
        "test1_moment0_summary",
        {
            "shape": tuple(int(size) for size in moment0.shape),
            "dtype": str(moment0.dtype),
        },
    )

    # Nothing is written to the output H5 because this tutorial only uses
    # ctx.state. Use ctx.output.h5.write(...) for persistent output.
