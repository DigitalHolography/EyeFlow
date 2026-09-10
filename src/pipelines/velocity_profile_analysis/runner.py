"""Publish weighted quadratic analysis for artery and vein profiles."""

from time import perf_counter

import h5py

from pipeline_engine.base import DatasetValue
from utils.logger import Logger

from .fitting import analyze_velocity_profiles

SOURCE_PATHS = {
    "Artery": "/Processing/VelocityProfiles/Artery/TransverseVelocityProfileMasked/value",
    "Vein": "/Processing/VelocityProfiles/Vein/TransverseVelocityProfileMasked/value",
}
OUTPUT_ROOT = "/Processing/VelocityProfileAnalysis"


def run_velocity_profile_analysis(ctx) -> dict[str, object]:
    """Fit both vessel classes from bounded slabs of published profiles."""

    datasets = {}
    for vessel, source_path in SOURCE_PATHS.items():
        dataset = ctx.output.h5.get(source_path)
        if not isinstance(dataset, h5py.Dataset):
            raise KeyError(
                f"Required {vessel.lower()} velocity-profile dataset is missing: "
                f"{source_path}. waveform_velocity must publish both vessel profiles "
                "before this analysis."
            )
        datasets[vessel] = dataset

    outputs: dict[str, object] = {}
    for vessel, dataset in datasets.items():
        started = perf_counter()
        results = analyze_velocity_profiles(dataset)
        Logger.log(
            f"Completed {vessel.lower()} weighted velocity-profile analysis in "
            f"{perf_counter() - started:.1f}s."
        )
        attrs = {
            "dimDesc": ["time", "beat", "branch", "radius"],
            "source_path": SOURCE_PATHS[vessel],
            "index_base": 0,
            "model": "a*x^2 + b*x + c",
            "fit_method": "weighted_least_squares",
            "weight_definition": (
                "u=x/(Nx-1); w=0.5 if u<0.25 or u>0.75, else w=1"
            ),
            "integration_method": (
                "unweighted_sum_of_finite_observed_integer_indexes_between_roots"
            ),
            "geometry_policy": "downward_opening_only; fractional_indexes",
            "vessel": vessel.lower(),
        }
        for name, value in results.items():
            outputs[f"{OUTPUT_ROOT}/{vessel}/{name}/value"] = DatasetValue(
                data=value,
                attrs=dict(attrs),
            )
    return outputs


__all__ = ["OUTPUT_ROOT", "SOURCE_PATHS", "run_velocity_profile_analysis"]
