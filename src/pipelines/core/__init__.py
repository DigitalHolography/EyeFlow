from .base import MissingPipeline, ProcessPipeline, ProcessResult
from .dag import PipelineDAG, PipelineExecutionPlan


def safe_h5_key(*args, **kwargs):
    from input_output import safe_h5_key as _safe_h5_key

    return _safe_h5_key(*args, **kwargs)


def write_result_h5(*args, **kwargs):
    from input_output import write_result_h5 as _write_result_h5

    return _write_result_h5(*args, **kwargs)


def write_combined_results_h5(*args, **kwargs):
    from input_output import write_combined_results_h5 as _write_combined_results_h5

    return _write_combined_results_h5(*args, **kwargs)


__all__ = [
    "ProcessPipeline",
    "MissingPipeline",
    "ProcessResult",
    "PipelineDAG",
    "PipelineExecutionPlan",
    "safe_h5_key",
    "write_result_h5",
    "write_combined_results_h5",
]
