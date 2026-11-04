"""Scientific stages with explicit products and owned local replay resources."""

from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from data_manager import InvocationResource
from input_output.readers.frame_sequences import FrameSequence
from pipeline_engine.storage import ArrayReference

from .parameters import MotionMapConfig
from .runner import DisplacementMapInputs, _run_displacement_mask, resolve_frame_rate, resolve_retina_masks


@dataclass
class RegistrationResources:
    scope: ExitStack
    directory: Path
    base: MotionMapConfig
    moment: object
    fps: float
    replay: object

    def close(self):
        self.scope.close()


def prepare(config, ctx):
    ctx.require_inputs("hd", "dv")
    moment = ctx.inputs.hd.as_holodoppler().named_root_moment_dataset(config.moment_path)
    fps = resolve_frame_rate(ctx, config.fallback_fps)
    # Run-owned backing survives producer termination and replay cleanup.
    directory = ctx.runtime.work_h5.storage.path("")
    directory.mkdir()
    base = MotionMapConfig(input=Path(moment.file.filename), output_dir=directory,
                           h5_dataset=moment.name.lstrip("/"), h5_fps=fps,
                           registration_method=config.registration_method, save_field=True)
    sequence = FrameSequence(base.input, base.h5_dataset, base.h5_frame_axis, base.h5_fps,
                             base.h5_low_percentile, base.h5_high_percentile,
                             h5_source=moment, defer_h5_range=True)
    scope = ExitStack()
    try:
        replay = scope.enter_context(sequence.cached_for_replay(base.max_frames, directory))
        resources = RegistrationResources(scope, directory, base, moment, fps, replay)
        ctx.dag.handoff("prepared", InvocationResource(resources))
    except BaseException:
        scope.close()
        raise


def vessel(name, names, ctx):
    prepared = ctx.dag.require("prepared")
    try:
        masks = resolve_retina_masks(ctx.inputs.dv.h5.h5file, prepared.moment.shape[-2:], name,
                                     source=ctx.inputs.dv.as_dopplerview())
        inputs = DisplacementMapInputs(prepared.moment, masks, prepared.fps)
        fields, videos = {}, []
        _run_displacement_mask(ctx, masks[0], inputs, prepared.base, prepared.replay,
                               prepared.directory, fields, videos, multiple_masks=len(names) > 1)
    except Exception as exc:
        ctx.log_warning(f"Skipping displacement {name}: {exc}")
        ctx.dag.skip(str(exc))
    references = {}
    for vessel_name, path in fields.items():
        array = np.load(path, mmap_mode="r", allow_pickle=False)
        try:
            references[vessel_name] = ArrayReference(str(path), tuple(array.shape), array.dtype)
        finally:
            array._mmap.close()
    ctx.dag.handoff(f"{name}_fields", references)
    for video in videos:
        ctx.log(f"Displacement magnitude video written to {video}.")


def finalize(config, names, ctx):
    available = {}
    for name in names:
        available.update(ctx.dag.get(f"{name}_fields", {}))
    if not available:
        raise RuntimeError("No displacement vessel field was produced.")
    ctx.dag.handoff("registered_fields", {"registration_method": config.registration_method,
                                          "fields": available})


def analyze(retain_segment_maps, ctx):
    from pipelines.topology_core.runner import prepared_topologies
    from .segments import DEFAULT_WORKING_MEMORY_MB, analyze_displacement_vessels
    from .timing import select_displacement_timing

    registration = ctx.dag.require("registration")
    timing = select_displacement_timing(ctx)
    if timing is None:
        ctx.dag.skip("Detected cardiac timing is unavailable in both velocity workflows.")
    if timing.workflow == "frequency_bands":
        ctx.log_warning("Displacement profiles use frequency-band cardiac timing because valid Doppler-moments timing is unavailable.")
    scratch = ctx.runtime.scratch
    memory = DEFAULT_WORKING_MEMORY_MB if scratch is None else min(DEFAULT_WORKING_MEMORY_MB, scratch.memory_budget_bytes / 2**20)
    analysis = analyze_displacement_vessels(registration["fields"], prepared_topologies(ctx),
        retain_maps=retain_segment_maps, working_memory_mb=memory,
        array_allocator=scratch.empty_float32 if scratch is not None else None)
    ctx.dag.handoff("analysis", {"analysis": analysis, "boundaries": timing.boundaries, "timing_attrs": timing.attrs,
                               "registration_method": registration["registration_method"]})


def profiles(ctx):
    from .profile_outputs import pack_displacement_profile_outputs
    value = ctx.dag.require("analysis")
    ctx.dag.handoff("profiles", pack_displacement_profile_outputs(value["analysis"], value["registration_method"],
        value["boundaries"], timing_attrs=value["timing_attrs"]))


def segment_maps(ctx):
    from .profile_outputs import pack_displacement_segment_map_outputs
    value = ctx.dag.require("analysis")
    ctx.dag.handoff("segment_maps", pack_displacement_segment_map_outputs(value["analysis"], value["registration_method"],
        value["boundaries"], timing_attrs=value["timing_attrs"]))


def packed(value, ctx):
    return value
