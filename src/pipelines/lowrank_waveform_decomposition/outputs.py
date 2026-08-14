"""Pack low-rank waveform decomposition endpoints and figure data."""

from __future__ import annotations

import warnings
from collections.abc import Mapping

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import DatasetValue, with_attrs
from pipelines.waveform_velocity_core.regions import (
    QUADRANTS_GROUP_NAME,
    REGION_NAMES,
    normalize_spatial_frame,
    optic_disc_center_xy,
    region_membership,
)

from .calculator import (
    LowRankWaveformDecompositionCalculator,
    ensure_segment_shape,
    mean_subtract,
)

ENDPOINT_METRICS = (
    "A1",
    "A2",
    "R0",
    "R1",
    "R2",
    "rho1",
    "rho2",
    "MPR",
    "Reff",
    "PR",
)

_ENDPOINT_KEY_MAP = {
    "A1": "A1",
    "A2": "A2",
    "R0": "R0",
    "R1": "R1",
    "R2": "R2",
    "rho1": "rho1",
    "rho2": "rho2",
    "MPR": "MPR",
    "Reff": "effective_rank",
    "PR": "participation_ratio",
}

_PER_BEAT_KEY_MAP = {
    "A1": "A1_b_pb",
    "A2": "A2_b_pb",
    "R0": "R0_b_pb",
    "R1": "R1_b_pb",
    "R2": "R2_b_pb",
    "rho1": "rho1_b_pb",
    "rho2": "rho2_b_pb",
    "MPR": "MPR_b_pb",
    "Reff": "effective_rank_b_pb",
    "PR": "participation_ratio_b_pb",
}

_VELOCITY_METRICS = frozenset({"A1", "A2", "R0", "R1", "R2"})

_CROSS_COLUMN_STATS = {
    "mean": "cross_column_mean",
    "lo": "cross_column_lo",
    "hi": "cross_column_hi",
}


def pack_lowrank_waveform_decomposition_outputs(
    velocity_outputs: Mapping[str, object],
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    vein_flag: bool = True,
    include_quadrants: bool = False,
    source_data=None,
    artery_segments=None,
    vein_segments=None,
) -> dict[str, object]:
    """Compute the requested low-rank endpoint and per-beat metric outputs."""
    schema = _resolve_output_paths(output_paths)
    periods = _required_array(velocity_outputs, schema.beat_period_seconds)
    calculator = LowRankWaveformDecompositionCalculator()
    outputs: dict[str, object] = {}

    quadrant_memberships: dict[str, np.ndarray] = {}
    segments_by_vessel = {
        "artery": artery_segments,
        "vein": vein_segments,
    }
    for vessel_name, signal_name, dataset_path in _source_paths(
        schema,
        vein_flag=vein_flag,
    ):
        if dataset_path is None or dataset_path not in velocity_outputs:
            continue

        waveforms = _metric_data(velocity_outputs[dataset_path])
        waveforms = ensure_segment_shape(waveforms, periods)
        prefix = (
            f"{schema.lowrank_waveform_decomposition_root}/"
            f"{vessel_name}/{signal_name}"
        )
        _compute_and_append_outputs(
            outputs,
            prefix,
            calculator,
            waveforms,
            periods,
        )

        if not include_quadrants:
            continue
        membership = quadrant_memberships.get(vessel_name)
        if membership is None:
            membership = _quadrant_membership(
                source_data,
                segments_by_vessel[vessel_name],
                waveforms,
                vessel_name,
            )
            quadrant_memberships[vessel_name] = membership
        for region_index, region_name in enumerate(REGION_NAMES):
            selected = membership[region_index]
            quadrant_waveforms = np.where(
                selected[np.newaxis, np.newaxis, :, :],
                waveforms,
                np.nan,
            )
            quadrant_prefix = (
                f"{schema.lowrank_waveform_decomposition_root}/{vessel_name}/"
                f"{QUADRANTS_GROUP_NAME}/{region_name}/{signal_name}"
            )
            _compute_and_append_outputs(
                outputs,
                quadrant_prefix,
                calculator,
                quadrant_waveforms,
                periods,
            )

    return outputs


def _compute_and_append_outputs(
    outputs: dict[str, object],
    prefix: str,
    calculator: LowRankWaveformDecompositionCalculator,
    waveforms: np.ndarray,
    periods: np.ndarray,
) -> None:
    representation = calculator.compute(waveforms, periods)
    per_beat_panels = calculator.per_beat_svd_panels(waveforms, periods)
    per_beat = calculator._compute_per_beat_endpoints(per_beat_panels)
    _append_requested_outputs(
        outputs,
        prefix,
        representation,
        per_beat,
        per_beat_panels,
        waveforms,
    )


def _quadrant_membership(
    source_data,
    segments,
    waveforms: np.ndarray,
    vessel_name: str,
) -> np.ndarray:
    if source_data is None or segments is None:
        raise RuntimeError(
            f"Quadrant low-rank outputs require {vessel_name} segment geometry."
        )
    branch_ids = np.asarray(segments.branch_ids, dtype=np.int32).reshape(-1)
    labels = np.asarray(segments.labels, dtype=np.int32)
    centers = np.asarray(segments.segment_center_xy, dtype=float)
    center_xy = optic_disc_center_xy(source_data, labels.shape)
    labels, center_xy = normalize_spatial_frame(labels, center_xy)
    membership = region_membership(branch_ids, labels, centers, center_xy)
    if membership.shape[1:] != waveforms.shape[2:]:
        raise ValueError(
            f"{vessel_name.capitalize()} quadrant membership shape "
            f"{membership.shape[1:]} does not match low-rank waveform "
            f"branch/radius shape {waveforms.shape[2:]}."
        )
    return membership


def _source_paths(
    schema: EyeFlowOutputPaths,
    *,
    vein_flag: bool,
) -> tuple[tuple[str, str, str | None], ...]:
    sources: list[tuple[str, str, str | None]] = [
        ("artery", "raw", schema.artery_per_beat.segment_velocity_signal),
        (
            "artery",
            "bandlimited",
            schema.artery_per_beat.segment_velocity_signal_band_limited,
        ),
    ]
    if vein_flag:
        sources.extend(
            (
                ("vein", "raw", schema.vein_per_beat.segment_velocity_signal),
                (
                    "vein",
                    "bandlimited",
                    schema.vein_per_beat.segment_velocity_signal_band_limited,
                ),
            )
        )
    return tuple(sources)


def _append_requested_outputs(
    metrics: dict[str, object],
    prefix: str,
    rep: dict,
    per_beat: Mapping[str, np.ndarray],
    per_beat_panels: Mapping[str, np.ndarray],
    waveforms: np.ndarray,
) -> None:
    acq = rep.get("acq", {})
    for output_name in ENDPOINT_METRICS:
        metrics[f"{prefix}/endpoints/joint/{output_name}"] = _metric(
            acq.get(_ENDPOINT_KEY_MAP[output_name], np.nan),
            unit=_unit_for(output_name),
        )

    for output_name in ENDPOINT_METRICS:
        metrics[f"{prefix}/endpoints/per_beat/{output_name}"] = _metric(
            per_beat.get(_PER_BEAT_KEY_MAP[output_name], np.nan),
            unit=_unit_for(output_name),
            dims=("beat",),
        )

    metrics[f"{prefix}/beat_period/mean"] = _metric(
        acq.get("beat_period_mean", np.nan), unit="s"
    )
    metrics[f"{prefix}/beat_period/std"] = _metric(
        acq.get("beat_period_std", np.nan), unit="s"
    )
    metrics[f"{prefix}/baseline/mu_acq"] = _metric(
        acq.get("mu_acq", np.nan), unit="mm/s"
    )
    metrics[f"{prefix}/baseline/sigma_mu_beat"] = _metric(
        acq.get("sigma_mu_beat", np.nan), unit="mm/s"
    )

    _append_acquisition_level_velocity(metrics, prefix, rep, waveforms)
    _append_waveform_components_joint(metrics, prefix, rep, waveforms)
    _append_waveform_components_per_beat(
        metrics, prefix, per_beat_panels, waveforms
    )
    _append_svd_spectrum(metrics, prefix, rep, per_beat)


def _append_acquisition_level_velocity(
    metrics: dict[str, object],
    prefix: str,
    rep: Mapping[str, object],
    waveforms: np.ndarray,
) -> None:
    group = f"{prefix}/misc/acquisition_level_velocity"
    metrics[f"{group}/cardiac_phase"] = _metric(
        _time_axis(waveforms), dims=("sample",)
    )
    summary = _summarize_time_columns(waveforms, _valid_mask(rep, waveforms))
    metrics[f"{group}/velocity_cross_column_mean"] = _metric(
        summary["mean"], unit="mm/s", dims=("sample",)
    )
    metrics[f"{group}/velocity_cross_column_std"] = _metric(
        summary["std"], unit="mm/s", dims=("sample",)
    )


def _append_waveform_components_joint(
    metrics: dict[str, object],
    prefix: str,
    rep: Mapping[str, object],
    waveforms: np.ndarray,
) -> None:
    group = f"{prefix}/misc/waveform_components_joint"
    valid_mask = _valid_mask(rep, waveforms)
    metrics[f"{group}/cardiac_phase"] = _metric(
        _time_axis(waveforms), dims=("sample",)
    )
    for name, panel in _joint_decomposition_components(rep, waveforms).items():
        _append_component_summary(metrics, f"{group}/{name}", panel, valid_mask)


def _append_waveform_components_per_beat(
    metrics: dict[str, object],
    prefix: str,
    per_beat_panels: Mapping[str, np.ndarray],
    waveforms: np.ndarray,
) -> None:
    group = f"{prefix}/misc/waveform_components_per_beat"
    valid_mask = np.asarray(per_beat_panels.get("valid_mask", []), dtype=bool)
    if valid_mask.shape != waveforms.shape[1:]:
        valid_mask = np.zeros(waveforms.shape[1:], dtype=bool)
    metrics[f"{group}/cardiac_phase"] = _metric(
        _time_axis(waveforms), dims=("sample",)
    )
    for name, panel in _per_beat_decomposition_components(
        per_beat_panels, waveforms
    ).items():
        _append_component_summary(metrics, f"{group}/{name}", panel, valid_mask)


def _append_svd_spectrum(
    metrics: dict[str, object],
    prefix: str,
    rep: Mapping[str, object],
    per_beat: Mapping[str, np.ndarray],
) -> None:
    energy_fraction = np.asarray(rep.get("energy_fraction", []), dtype=float)
    mode = np.arange(1, energy_fraction.size + 1, dtype=int)

    group = f"{prefix}/misc/svd_spectrum_joint"
    metrics[f"{group}/svd_mode_index"] = _metric(mode, dims=("mode",))
    metrics[f"{group}/lambda"] = _metric(energy_fraction, dims=("mode",))

    pb_lambda = _per_beat_energy_fraction(per_beat.get("singular_values_b", []))
    pb_mode = np.arange(1, pb_lambda.shape[1] + 1, dtype=int)
    pb_summary = _summarize_mode_rows(pb_lambda)

    group = f"{prefix}/misc/svd_spectrum_joint_per_beat"
    metrics[f"{group}/svd_mode_index"] = _metric(pb_mode, dims=("mode",))
    metrics[f"{group}/lambda_mean"] = _metric(pb_summary["mean"], dims=("mode",))
    metrics[f"{group}/lambda_lo"] = _metric(pb_summary["lo"], dims=("mode",))
    metrics[f"{group}/lambda_hi"] = _metric(pb_summary["hi"], dims=("mode",))


def _joint_decomposition_components(
    rep: Mapping[str, object],
    waveforms: np.ndarray,
) -> dict[str, np.ndarray]:
    mu = np.asarray(rep.get("mu", np.full(waveforms.shape[1:], np.nan)), dtype=float)
    x = np.asarray(rep.get("x_full", np.full_like(waveforms, np.nan)), dtype=float)
    return {
        "segment_velocity": np.asarray(waveforms, dtype=float),
        "temporal_baseline": np.broadcast_to(mu[None, :, :, :], waveforms.shape),
        "centered_velocity": x,
        "svd_mode1": _joint_mode_component(rep, waveforms, mode=1),
        "svd_mode2": _joint_mode_component(rep, waveforms, mode=2),
    }


def _per_beat_decomposition_components(
    per_beat_panels: Mapping[str, np.ndarray],
    waveforms: np.ndarray,
) -> dict[str, np.ndarray]:
    mu, x = mean_subtract(waveforms)
    return {
        "segment_velocity": np.asarray(waveforms, dtype=float),
        "temporal_baseline": np.broadcast_to(mu[None, :, :, :], waveforms.shape),
        "centered_velocity": x,
        "svd_mode1": _per_beat_mode_component(per_beat_panels, waveforms, mode=1),
        "svd_mode2": _per_beat_mode_component(per_beat_panels, waveforms, mode=2),
    }


def _joint_mode_component(
    rep: Mapping[str, object],
    waveforms: np.ndarray,
    *,
    mode: int,
) -> np.ndarray:
    component = np.full(waveforms.shape, np.nan, dtype=float)
    valid_mask = _valid_mask(rep, waveforms)
    u_panel = np.asarray(rep.get("U_panel", []), dtype=float)
    scores = np.asarray(rep.get("score_panel_bkr", []), dtype=float)
    index = mode - 1
    if (
        u_panel.ndim == 2
        and scores.ndim == 4
        and u_panel.shape[0] == waveforms.shape[0]
        and index < u_panel.shape[1]
        and index < scores.shape[0]
    ):
        panel = u_panel[:, index, None, None, None] * scores[index][None, :, :, :]
        component[:, valid_mask] = panel[:, valid_mask]
    return component


def _per_beat_mode_component(
    per_beat_panels: Mapping[str, np.ndarray],
    waveforms: np.ndarray,
    *,
    mode: int,
) -> np.ndarray:
    component = np.full(waveforms.shape, np.nan, dtype=float)
    valid_mask = np.asarray(per_beat_panels.get("valid_mask", []), dtype=bool)
    u_mode_tb = np.asarray(per_beat_panels.get("u_mode_tb", []), dtype=float)
    scores = np.asarray(per_beat_panels.get("scores_mode_bkr", []), dtype=float)
    index = mode - 1
    if (
        valid_mask.shape == waveforms.shape[1:]
        and u_mode_tb.ndim == 3
        and scores.ndim == 4
        and u_mode_tb.shape[1] == waveforms.shape[0]
        and index < u_mode_tb.shape[0]
        and index < scores.shape[0]
    ):
        panel = (
            u_mode_tb[index, :, :, None, None]
            * scores[index][None, :, :, :]
        )
        component[:, valid_mask] = panel[:, valid_mask]
    return component


def _append_component_summary(
    metrics: dict[str, object],
    group: str,
    panel: np.ndarray,
    valid_mask: np.ndarray,
) -> None:
    summary = _summarize_time_columns(panel, valid_mask)
    for stat, output_name in _CROSS_COLUMN_STATS.items():
        metrics[f"{group}/{output_name}"] = _metric(
            summary[stat],
            unit="mm/s",
            dims=("sample",),
        )


def _summarize_time_columns(
    panel: np.ndarray,
    valid_mask: np.ndarray,
) -> dict[str, np.ndarray]:
    values = np.asarray(panel, dtype=float)
    if (
        values.ndim != 4
        or valid_mask.shape != values.shape[1:]
        or not np.any(valid_mask)
    ):
        n_t = values.shape[0] if values.ndim > 0 else 0
        empty = np.full((n_t,), np.nan, dtype=float)
        return {
            "mean": empty,
            "std": empty.copy(),
            "lo": empty.copy(),
            "hi": empty.copy(),
        }

    columns = values[:, valid_mask]
    mean = _nan_reducer(columns, np.nanmean, axis=1)
    std = _nan_sample_std(columns, axis=1)
    return {
        "mean": mean,
        "std": std,
        "lo": mean - std,
        "hi": mean + std,
    }


def _summarize_mode_rows(values: np.ndarray) -> dict[str, np.ndarray]:
    values = np.asarray(values, dtype=float)
    if values.ndim != 2:
        values = np.empty((0, 0), dtype=float)
    mean = _nan_reducer(values, np.nanmean, axis=0)
    std = _nan_sample_std(values, axis=0)
    return {
        "mean": mean,
        "lo": mean - std,
        "hi": mean + std,
    }


def _per_beat_energy_fraction(singular_values_b) -> np.ndarray:
    s = np.asarray(singular_values_b, dtype=float)
    if s.ndim != 2:
        return np.empty((0, 0), dtype=float)
    energy = s**2
    total = np.nansum(energy, axis=1, keepdims=True)
    return np.divide(
        energy,
        total,
        out=np.full_like(energy, np.nan, dtype=float),
        where=np.isfinite(total) & (total > 0),
    )


def _valid_mask(rep: Mapping[str, object], waveforms: np.ndarray) -> np.ndarray:
    valid_mask = np.asarray(rep.get("valid_column_mask", []), dtype=bool)
    if valid_mask.shape != waveforms.shape[1:]:
        return np.zeros(waveforms.shape[1:], dtype=bool)
    return valid_mask


def _time_axis(waveforms: np.ndarray) -> np.ndarray:
    n_t = int(np.asarray(waveforms).shape[0])
    return np.linspace(0.0, 1.0, n_t, endpoint=False, dtype=float)


def _nan_reducer(values: np.ndarray, reducer, *, axis: int) -> np.ndarray:
    if values.size == 0 or not np.any(np.isfinite(values)):
        shape = tuple(dim for index, dim in enumerate(values.shape) if index != axis)
        return np.full(shape, np.nan, dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.asarray(reducer(values, axis=axis), dtype=float)


def _nan_sample_std(values: np.ndarray, *, axis: int) -> np.ndarray:
    """Sample SD (``ddof=1``); zero when fewer than two finite samples."""
    if values.size == 0 or not np.any(np.isfinite(values)):
        shape = tuple(dim for index, dim in enumerate(values.shape) if index != axis)
        return np.full(shape, np.nan, dtype=float)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        std = np.asarray(np.nanstd(values, axis=axis, ddof=1), dtype=float)
    n = np.sum(np.isfinite(values), axis=axis)
    return np.where(n > 1, std, 0.0)


def _unit_for(metric_name: str) -> str | None:
    return "mm/s" if metric_name in _VELOCITY_METRICS else None


def _metric(
    data,
    *,
    unit: str | None = None,
    dims: tuple[str, ...] | None = None,
):
    attrs: dict[str, object] = {}
    if unit:
        attrs["unit"] = unit
    if dims:
        attrs["dimDesc"] = list(dims)
    return with_attrs(np.asarray(data), attrs) if attrs else np.asarray(data)


def _required_array(metrics: Mapping[str, object], path: str) -> np.ndarray:
    if path not in metrics:
        raise RuntimeError(
            f"Required shared per-beat output '{path}' is unavailable; "
            "check the waveform_velocity DAG dependency."
        )
    return _metric_data(metrics[path])


def _metric_data(value) -> np.ndarray:
    if isinstance(value, DatasetValue):
        value = value.data
    elif isinstance(value, tuple) and len(value) == 2 and isinstance(value[1], dict):
        value = value[0]
    elif hasattr(value, "data") and hasattr(value, "attrs"):
        value = value.data
    return np.asarray(value)


def _resolve_output_paths(
    output_paths: EyeFlowOutputPaths | str | None,
) -> EyeFlowOutputPaths:
    if isinstance(output_paths, EyeFlowOutputPaths):
        return output_paths
    return EyeFlowOutputPaths.active(output_paths)


__all__ = ["pack_lowrank_waveform_decomposition_outputs"]
