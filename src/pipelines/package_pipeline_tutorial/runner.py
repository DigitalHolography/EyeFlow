"""Implementation module for the package pipeline tutorial."""


def run_package_pipeline_tutorial(ctx) -> None:
    """No-op pipeline body.

    Add real logic here when the pipeline grows beyond a single file. The
    runtime gives this function the same PipelineContext object as any other
    pipeline, so it can use typed input loaders or ctx.sources for explicit
    H5-path reads and write outputs with ctx.write()/ctx.write_many().

    Use ctx.vars or ctx.set_var(...) for values that later pipelines should
    see but that should not be written to the output H5 file.
    """

    ctx.set_var("package_pipeline_tutorial_ran", True)
    return None
