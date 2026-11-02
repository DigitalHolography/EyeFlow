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
        try:
            moment = prepared["moment"]
            masks = resolve_retina_masks(
                ctx.inputs.dv.h5.h5file,
                moment.shape[-2:],
                name,
                source=ctx.inputs.dv.as_dopplerview(),
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
        except Exception as exc:
            ctx.log_warning(f"Skipping displacement {name}: {exc}")
            ctx.dag.skip(str(exc))
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

    def finalize(ctx):
        available = {}
        for name in names:
            product = ctx.dag.get(f"{name}_fields")
            if product is not None:
                available.update(product["fields"])
        if not available:
            raise RuntimeError("No displacement vessel field was produced.")
        ctx.dag.handoff("fields_ready", tuple(sorted(available)))

    def cleanup():
        scope.close()

    from functools import partial

    nodes = [
        PipelineNode("prepare", prepare, outputs=(PipelineProduct("prepared"),)),
        *(
            PipelineNode(
                name,
                partial(vessel, name),
                requires=("prepared",),
                outputs=(PipelineProduct(f"{name}_fields", optional=True),),
            )
            for name in names
        ),
        PipelineNode(
            "finalize",
            finalize,
            optional_requires=tuple(f"{name}_fields" for name in names),
            outputs=(PipelineProduct("fields_ready"),),
        ),
    ]
    return PipelineGraph(tuple(nodes), cleanup=cleanup)


def displacement_profiles_graph(*, retain_segment_maps=False):
    prepared = {}

    def analyze(ctx):
        from pipelines.topology_core.runner import prepared_topologies

        from .runner import DISPLACEMENT_MAP_STATE, DisplacementMapArtifacts
        from .segments import DEFAULT_WORKING_MEMORY_MB, analyze_displacement_vessels
        from .timing import select_displacement_timing

        artifacts = ctx.state.get(DISPLACEMENT_MAP_STATE)
        if not isinstance(artifacts, DisplacementMapArtifacts):
            raise RuntimeError("Displacement fields are unavailable for profile analysis.")
        timing = select_displacement_timing(ctx)
        if timing is None:
            ctx.dag.skip(
                "Detected cardiac timing is unavailable in both velocity workflows."
            )
        if timing.workflow == "frequency_bands":
            ctx.log_warning(
                "Displacement profiles use frequency-band cardiac timing because "
                "valid Doppler-moments timing is unavailable."
            )
        scratch = ctx.runtime.scratch
        working_memory_mb = (
            DEFAULT_WORKING_MEMORY_MB
            if scratch is None
            else min(DEFAULT_WORKING_MEMORY_MB, scratch.memory_budget_bytes / 2**20)
        )
        prepared["analysis"] = analyze_displacement_vessels(
            artifacts.field_paths_by_vessel,
            prepared_topologies(ctx),
            retain_maps=retain_segment_maps,
            working_memory_mb=working_memory_mb,
            array_allocator=(
                scratch.empty_float32
                if scratch is not None
                else None
            ),
        )
        prepared["timing"] = timing
        prepared["registration_method"] = artifacts.registration_method
        ctx.dag.handoff("profiles", True)
        if retain_segment_maps:
            ctx.dag.handoff("segment_maps", True)

    def pack_profiles(_value, _ctx):
        from .profile_outputs import pack_displacement_profile_outputs

        timing = prepared["timing"]
        return pack_displacement_profile_outputs(
            prepared["analysis"],
            prepared["registration_method"],
            timing.boundaries,
            timing_attrs=timing.attrs,
        )

    def pack_segment_maps(_value, _ctx):
        from .profile_outputs import pack_displacement_segment_map_outputs

        timing = prepared["timing"]
        return pack_displacement_segment_map_outputs(
            prepared["analysis"],
            prepared["registration_method"],
            timing.boundaries,
            timing_attrs=timing.attrs,
        )

    products = [PipelineProduct("profiles", pack=pack_profiles)]
    if retain_segment_maps:
        products.append(PipelineProduct("segment_maps", pack=pack_segment_maps))
    return PipelineGraph(
        (
            PipelineNode(
                "analysis",
                analyze,
                outputs=tuple(products),
            ),
        )
    )
