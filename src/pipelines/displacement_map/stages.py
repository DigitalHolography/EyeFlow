"""Independent vessel products over one shared local frame replay."""

from contextlib import ExitStack
from pathlib import Path

from input_output.readers.frame_sequences import FrameSequence
from pipeline_engine.graphs import PipelineGraph, PipelineNode, PipelineProduct
from pipeline_engine.storage import temporary_directory

from .parameters import MotionMapConfig
from .runner import (
    DISPLACEMENT_MAP_STATE,
    DisplacementMapArtifacts,
    DisplacementMapInputs,
    _run_displacement_mask,
    resolve_frame_rate,
    resolve_retina_masks,
)


def displacement_graph(config):
    scope = ExitStack()
    prepared = {}
    names = ("artery", "vein") if config.mask_mode == "both" else (config.mask_mode,)

    def prepare(ctx):
        ctx.require_inputs("hd", "dv")
        moment = ctx.inputs.hd.as_holodoppler().named_root_moment_dataset(config.moment_path)
        fps = resolve_frame_rate(ctx, config.fallback_fps)
        directory = temporary_directory("eyeflow-displacement-map-")
        base = MotionMapConfig(
            input=Path(moment.file.filename),
            output_dir=Path(directory.name),
            h5_dataset=moment.name.lstrip("/"),
            h5_fps=fps,
            registration_method=config.registration_method,
            save_field=True,
        )
        sequence = FrameSequence(
            base.input,
            base.h5_dataset,
            base.h5_frame_axis,
            base.h5_fps,
            base.h5_low_percentile,
            base.h5_high_percentile,
            h5_source=moment,
            defer_h5_range=True,
        )
        prepared.update(directory=directory, base=base, moment=moment, fps=fps)
        prepared["replay"] = scope.enter_context(
            sequence.cached_for_replay(base.max_frames, Path(directory.name))
        )
        ctx.dag.handoff("prepared", True)

    def vessel(name, ctx):
        ctx.dag.require("prepared")
        moment = prepared["moment"]
        masks = resolve_retina_masks(
            ctx.inputs.dv.h5.h5file, moment.shape[-2:], name, source=ctx.inputs.dv.as_dopplerview()
        )
        inputs = DisplacementMapInputs(moment, masks, prepared["fps"])
        fields, videos = {}, []
        _run_displacement_mask(
            ctx,
            masks[0],
            inputs,
            prepared["base"],
            prepared["replay"],
            Path(prepared["directory"].name),
            fields,
            videos,
            multiple_masks=len(names) > 1,
        )
        previous = ctx.state.get(DISPLACEMENT_MAP_STATE)
        paths = dict(previous.field_paths_by_vessel) if previous else {}
        paths.update(fields)
        ctx.state.set(
            DISPLACEMENT_MAP_STATE,
            DisplacementMapArtifacts(config.registration_method, paths, prepared["directory"]),
        )
        ctx.dag.handoff(f"{name}_fields", {"fields": fields, "videos": videos})
        for video in videos:
            ctx.log(f"Displacement magnitude video written to {video}.")

    def cleanup():
        scope.close()

    from functools import partial

    return PipelineGraph(
        (
            PipelineNode("prepare", prepare, outputs=(PipelineProduct("prepared"),)),
            *(
                PipelineNode(
                    name,
                    partial(vessel, name),
                    requires=("prepared",),
                    outputs=(PipelineProduct(f"{name}_fields"),),
                )
                for name in names
            ),
        ),
        cleanup=cleanup,
    )
