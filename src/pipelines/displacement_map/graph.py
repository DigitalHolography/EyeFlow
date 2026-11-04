"""Registration and profile graphs; calculations and resources live in stages."""

from functools import partial

from pipeline_engine.stage_graph import ProductSpec, StageGraph, StageSpec

from .runner import DisplacementMapPipelineConfig
from . import stages


def displacement_graph(config=None, *, options=()):
    config = config or DisplacementMapPipelineConfig()
    names = ("artery", "vein") if config.mask_mode == "both" else (config.mask_mode,)
    return StageGraph((
        StageSpec("prepare", partial(stages.prepare, config), products=(ProductSpec("prepared"),)),
        *(StageSpec(name, partial(stages.vessel, name, names), requires=("prepared",),
                    products=(ProductSpec(f"{name}_fields", optional=True),)) for name in names),
        StageSpec("finalize", partial(stages.finalize, config, names),
                  optional_requires=tuple(f"{name}_fields" for name in names),
                  products=(ProductSpec("registered_fields", exposed=True),)),
    ))


def displacement_profiles_graph(*, options=(), retain_segment_maps=None):
    retain_maps = "segment_maps" in options if retain_segment_maps is None else retain_segment_maps
    nodes = [
        StageSpec("analysis", partial(stages.analyze, retain_maps), requires=("registration",),
                  products=(ProductSpec("analysis"),)),
        StageSpec("profiles", stages.profiles, requires=("analysis",),
                  products=(ProductSpec("profiles", pack=stages.packed),)),
    ]
    if retain_maps:
        nodes.append(StageSpec("segment_maps", stages.segment_maps, requires=("analysis",),
                               products=(ProductSpec("segment_maps", pack=stages.packed),)))
    return StageGraph(tuple(nodes))
