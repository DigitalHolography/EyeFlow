"""Orchestrate velocity-profile analysis products."""

import h5py

from pipeline_engine.base import DatasetValue

from .fitting import analyze_velocity_profiles

SOURCE_PATH = "/Processing/VelocityProfiles/Artery/TransverseVelocityProfileMasked/value"
OUTPUT_ROOT = "/Processing/VelocityProfileAnalysis/Artery"


def run_velocity_profile_analysis(ctx) -> dict[str, object]:
    """Fit individual artery profiles already published in the work HDF5 file."""
    dataset = ctx.output.h5.get(SOURCE_PATH)
    if not isinstance(dataset, h5py.Dataset):
        raise KeyError(
            f"Required velocity-profile dataset is missing: {SOURCE_PATH}. "
            "waveform_velocity must publish velocity_profiles before this analysis."
        )
    results = analyze_velocity_profiles(dataset)
    attrs = {
        "dimDesc": ["time", "beat", "branch", "radius"],
        "source_path": SOURCE_PATH,
        "index_base": 0,
        "model": "a*x^2 + b*x + c",
        "fit_method": "weighted_least_squares",
        "weight_definition": "u=x/(Nx-1); w=0.5 if u<0.25 or u>0.75, else w=1",
        "integration_method": "unweighted_sum_of_finite_observed_integer_indexes_between_roots",
        "geometry_policy": "downward_opening_only; fractional_indexes",
    }
    return {
        f"{OUTPUT_ROOT}/{name}/value": DatasetValue(data=value, attrs=dict(attrs))
        for name, value in results.items()
    }
