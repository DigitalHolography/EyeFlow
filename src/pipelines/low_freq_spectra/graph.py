"""The seven-stage spectral graph, shared by preflight and execution."""

from pipeline_engine.stage_graph import StageGraph, StageSpec, ProductSpec
from .stages import _packed, _source, _sampling, _calibration, _spatial, _beats, _qc, _figures


def spectral_graph(*, options=()):
    product = lambda name: (ProductSpec(name, pack=_packed),)
    return StageGraph((
        StageSpec("source", _source, products=(ProductSpec("source"),)),
        StageSpec("sampling", _sampling, requires=("source",), products=product("sampling")),
        StageSpec("calibration", _calibration, requires=("source", "sampling"), products=product("calibration")),
        StageSpec("spatial", _spatial, requires=("source",), products=product("spatial")),
        StageSpec("beats", _beats, requires=("spatial",), optional_requires=("calibration",), products=product("beats")),
        StageSpec("qc", _qc, requires=("spatial",), optional_requires=("beats",), products=product("qc")),
        StageSpec("figures", _figures, requires=("spatial",), optional_requires=("calibration", "beats", "qc"), artifact=True),
    ))
