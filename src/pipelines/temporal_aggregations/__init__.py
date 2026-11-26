"""Export frame-axis summaries of HoloDoppler image videos."""

from calculations.temporal_statistics import temporal_image_statistics
from input_output.schema.holodoppler import HD_BAND_NAME_PATTERN
from pipeline_engine import InputAccess, InputUse, StageGraph, StageSpec, pipeline
from pipeline_engine.context import require_managed_runtime

from .display import png_pixels


def _export_images(ctx):
    ctx.require_inputs("hd")
    videos = ctx.inputs.hd.as_holodoppler().image_videos()
    if not videos or not ctx.output.available:
        ctx.dag.skip("No HoloDoppler image videos or artifact destination.")

    for name, video in videos.items():
        enhance = name in {"moment0", "moment1", "moment2", "color_Doppler_image"} or bool(
            HD_BAND_NAME_PATTERN.fullmatch(name)
        )
        images = temporal_image_statistics(video, display="per_frame" if enhance else "global")
        directory = f"temporal_aggregations/{name}"
        for statistic, image in images.items():
            ctx.output.write_png(png_pixels(image), f"{directory}/temporal_{statistic}.png")


def temporal_aggregations_graph(*, options=()):
    return StageGraph((StageSpec(
        "images", _export_images, artifact=True,
        inputs=(InputUse("hd.image_videos", InputAccess.STREAM),),
    ),))


@pipeline(
    name="temporal_aggregations",
    description="Export temporal mean, standard deviation, and entropy images from HoloDoppler videos.",
    requires=["numpy", "h5py", "skimage"],
    input_slot="hd",
    default_selected=True,
    run_once=True,
    graph_factory=temporal_aggregations_graph,
)
def run(ctx):
    require_managed_runtime(ctx, "temporal_aggregations")
    return temporal_aggregations_graph()
