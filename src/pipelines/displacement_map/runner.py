"""Displacement registration and profile stage graphs with owned replay resources."""

from __future__ import annotations

import math
from contextlib import ExitStack
from dataclasses import dataclass, replace
from functools import partial
from pathlib import Path
from typing import Literal

import h5py
import numpy as np

from data_manager import InvocationResource
from input_output.output_manager import OutputType
from input_output.readers.frame_sequences import FrameSequence
from input_output.readers.h5 import RawH5SourceReader
from input_output.schema import DopplerViewSource
from input_output.schema.doppler_view import (
    DV_ARTERY_MASK_PATH as ARTERY_MASK_PATH,
    DV_LABELED_VESSELS_PATH as LABELED_VESSELS_PATH,
    DV_VEIN_MASK_PATH as VEIN_MASK_PATH,
    DV_VESSEL_MASK_PATH as VESSEL_MASK_PATH,
)
from pipeline_engine.data.storage import ArrayReference
from pipeline_engine.stage_graph import ProductSpec, StageGraph, StageSpec

from .calculator import create_retinal_motion_map
from .constants import DEFAULT_REGISTRATION_METHOD, RegistrationMethod
from .parameters import MotionMapConfig


MAGNITUDE_VIDEO_FILENAME = "displacement_magnitude.avi"


DISPLACEMENT_AVI_SUBFOLDER = "displacement_maps"


DEFAULT_MOMENT_PATH = "moment0"


DEFAULT_FPS = 25.0


MaskMode = Literal["both", "combined", "artery", "vein", "labeled"]


@dataclass(frozen=True, slots=True)
class DisplacementMapPipelineConfig:
    """Input choices owned by the pipeline rather than by the algorithm."""

    moment_path: str = DEFAULT_MOMENT_PATH
    mask_mode: MaskMode = "both"
    fallback_fps: float = DEFAULT_FPS
    registration_method: RegistrationMethod = DEFAULT_REGISTRATION_METHOD


@dataclass(frozen=True, slots=True)
class DisplacementMaskInput:
    name: str
    vessels: tuple[str, ...]
    mask: np.ndarray
    source: str


@dataclass(frozen=True, slots=True)
class DisplacementMapInputs:
    moment: h5py.Dataset
    masks: tuple[DisplacementMaskInput, ...]
    fps: float


def run_displacement_map(ctx, config=None):
    """Return the registration graph; runtime owns its invocation."""
    from pipeline_engine.context import require_managed_runtime
    require_managed_runtime(ctx, "displacement_map")

    return displacement_graph(config or DisplacementMapPipelineConfig())


def run_displacement_profiles(ctx):
    """Return the product-dependent analysis graph."""
    from pipeline_engine.context import require_managed_runtime
    require_managed_runtime(ctx, "displacement_profiles")

    return displacement_profiles_graph(options=ctx.options_for("displacement_profiles"))


def _run_displacement_mask(
    ctx, mask_input, inputs, base_config, replay, temporary_path,
    field_paths_by_vessel, output_videos, *, multiple_masks=False,
) -> None:
    output_dir = temporary_path / mask_input.name
    output_dir.mkdir(parents=True, exist_ok=True)
    video_filename = (
        MAGNITUDE_VIDEO_FILENAME
        if len(inputs.masks) == 1 and not multiple_masks
        else f"{mask_input.name}_{MAGNITUDE_VIDEO_FILENAME}"
    )
    output_video = ctx.output.path_for(
        OutputType.AVI, f"{DISPLACEMENT_AVI_SUBFOLDER}/{video_filename}",
    )
    output_video.parent.mkdir(parents=True, exist_ok=True)
    outputs = create_retinal_motion_map(
        replace(base_config, output_dir=output_dir),
        analysis_mask_array=mask_input.mask,
        magnitude_video_path=output_video,
        h5_source=inputs.moment,
        frame_replay=replay,
        log_label=f"Displacement {mask_input.name}",
    )
    field_path = Path(outputs["displacement_field"])
    for vessel in mask_input.vessels:
        field_paths_by_vessel[vessel] = field_path
    output_videos.append(output_video)


def resolve_retina_masks(
    h5file: h5py.File | None,
    spatial_shape: tuple[int, int],
    mode: MaskMode,
    *,
    source: DopplerViewSource | None = None,
) -> tuple[DisplacementMaskInput, ...]:
    """Resolve one shared mask or separate artery and vein masks."""

    if h5file is None:
        raise ValueError("The DopplerView HDF5 input is required.")
    source = source or DopplerViewSource(RawH5SourceReader(h5file=h5file, label="DV"), {})
    if mode == "both":
        masks: list[DisplacementMaskInput] = []
        for vessel in ("artery", "vein"):
            mask, resolved_source = resolve_retina_mask(
                h5file,
                spatial_shape,
                vessel,
                source=source,
            )
            masks.append(
                DisplacementMaskInput(
                    name=vessel,
                    vessels=(vessel,),
                    mask=mask,
                    source=resolved_source,
                )
            )
        return tuple(masks)

    mask, resolved_source = resolve_retina_mask(h5file, spatial_shape, mode, source=source)
    if mode == "artery":
        vessels = ("artery",)
    elif mode == "vein":
        vessels = ("vein",)
    else:
        vessels = ("artery", "vein")
    return (
        DisplacementMaskInput(
            name=mode,
            vessels=vessels,
            mask=mask,
            source=resolved_source,
        ),
    )


def resolve_retina_mask(
    h5file: h5py.File | None,
    spatial_shape: tuple[int, int],
    mode: MaskMode = "combined",
    *,
    source: DopplerViewSource | None = None,
) -> tuple[np.ndarray, str]:
    """Resolve and align one DopplerView retina mask."""

    if h5file is None:
        raise ValueError("The DopplerView HDF5 input is required.")
    source = source or DopplerViewSource(RawH5SourceReader(h5file=h5file, label="DV"), {})
    if mode == "artery":
        mask = _required_mask(source, ARTERY_MASK_PATH)
        source_path = ARTERY_MASK_PATH
    elif mode == "vein":
        mask = _required_mask(source, VEIN_MASK_PATH)
        source_path = VEIN_MASK_PATH
    elif mode == "labeled":
        mask = _required_mask(source, LABELED_VESSELS_PATH)
        source_path = LABELED_VESSELS_PATH
    elif mode == "combined":
        mask, source_path = _combined_mask(source)
    else:
        raise ValueError(f"Unknown displacement-map mask mode: {mode!r}")

    aligned = _align_mask(mask, spatial_shape, source_path)
    if not np.any(aligned):
        raise ValueError(f"DopplerView mask '{source_path}' contains no vessel pixels.")
    return aligned, source_path


def _combined_mask(source: DopplerViewSource) -> tuple[np.ndarray, str]:
    vessel = _optional_mask(source, VESSEL_MASK_PATH)
    if vessel is not None:
        return vessel, VESSEL_MASK_PATH

    labeled = _optional_mask(source, LABELED_VESSELS_PATH)
    if labeled is not None:
        return labeled, LABELED_VESSELS_PATH

    artery = _optional_mask(source, ARTERY_MASK_PATH)
    vein = _optional_mask(source, VEIN_MASK_PATH)
    available = [mask for mask in (artery, vein) if mask is not None]
    if not available:
        raise KeyError(
            "Missing DopplerView retinal vessel mask. Tried: "
            f"{VESSEL_MASK_PATH}, {LABELED_VESSELS_PATH}, "
            f"{ARTERY_MASK_PATH}, {VEIN_MASK_PATH}"
        )
    combined = np.asarray(available[0]) != 0
    for mask in available[1:]:
        if np.asarray(mask).shape != combined.shape:
            raise ValueError("DopplerView artery and vein masks have different shapes.")
        combined |= np.asarray(mask) != 0
    sources = [
        path
        for path, mask in ((ARTERY_MASK_PATH, artery), (VEIN_MASK_PATH, vein))
        if mask is not None
    ]
    return combined, "+".join(sources)


def _required_mask(source: DopplerViewSource, path: str) -> np.ndarray:
    mask = _optional_mask(source, path)
    if mask is None:
        raise KeyError(f"Missing DopplerView mask dataset at '{path}'.")
    return mask


def _optional_mask(source: DopplerViewSource, path: str) -> np.ndarray | None:
    return source.retinal_mask(path)


def _align_mask(
    mask: np.ndarray,
    spatial_shape: tuple[int, int],
    source: str,
) -> np.ndarray:
    value = np.squeeze(np.asarray(mask))
    if value.ndim != 2:
        raise ValueError(
            f"DopplerView mask '{source}' must become 2-D after squeeze, got {value.shape}."
        )
    if value.shape == spatial_shape:
        return np.asarray(value != 0, dtype=bool)
    if value.T.shape == spatial_shape:
        return np.asarray(value.T != 0, dtype=bool)
    raise ValueError(
        f"DopplerView mask '{source}' shape {value.shape} does not match "
        f"HoloDoppler spatial shape {spatial_shape} in either axis order."
    )


def resolve_frame_rate(ctx, fallback: float = DEFAULT_FPS) -> float:
    """Use HoloDoppler timing metadata, falling back for older exports."""

    try:
        dt_seconds = float(ctx.inputs.hd.as_holodoppler().timing().dt_seconds)
        fps = 1.0 / dt_seconds
        if math.isfinite(fps) and fps > 0:
            return fps
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        pass
    fallback = float(fallback)
    if not math.isfinite(fallback) or fallback <= 0:
        raise ValueError("fallback_fps must be finite and greater than zero.")
    ctx.log_warning(f"Could not resolve HoloDoppler timing; using {fallback:g} fps.")
    return fallback


__all__ = [
    "DisplacementMapPipelineConfig",
    "resolve_retina_mask",
    "resolve_retina_masks",
    "run_displacement_map",
    "run_displacement_profiles",
]


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


def displacement_graph(config=None, *, options=()):
    config = config or DisplacementMapPipelineConfig()
    names = ("artery", "vein") if config.mask_mode == "both" else (config.mask_mode,)
    return StageGraph((
        StageSpec("prepare", partial(prepare, config), products=(ProductSpec("prepared"),)),
        *(StageSpec(name, partial(vessel, name, names), requires=("prepared",),
                    products=(ProductSpec(f"{name}_fields", optional=True),)) for name in names),
        StageSpec("finalize", partial(finalize, config, names),
                  optional_requires=tuple(f"{name}_fields" for name in names),
                  products=(ProductSpec("registered_fields", exposed=True),)),
    ))


def displacement_profiles_graph(*, options=(), retain_segment_maps=None):
    retain_maps = "segment_maps" in options if retain_segment_maps is None else retain_segment_maps
    nodes = [
        StageSpec("analysis", partial(analyze, retain_maps), requires=("registration",),
                  products=(ProductSpec("analysis"),)),
        StageSpec("profiles", profiles, requires=("analysis",),
                  products=(ProductSpec("profiles", pack=packed),)),
    ]
    if retain_maps:
        nodes.append(StageSpec("segment_maps", segment_maps, requires=("analysis",),
                               products=(ProductSpec("segment_maps", pack=packed),)))
    return StageGraph(tuple(nodes))
