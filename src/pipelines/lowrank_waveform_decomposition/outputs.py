"""Pack low-rank waveform decomposition products into EyeFlow paths."""

from __future__ import annotations

from collections.abc import Mapping

import numpy as np

from input_output.schema import EyeFlowOutputPaths
from pipeline_engine import DatasetValue, with_attrs

from .calculator import LowRankWaveformDecompositionCalculator


def pack_lowrank_waveform_decomposition_outputs(
    velocity_outputs: Mapping[str, object],
    output_paths: EyeFlowOutputPaths | str | None = None,
    *,
    include_veins: bool = False,
) -> dict[str, object]:
    """Compute and prefix low-rank products from shared per-beat outputs."""
    schema = _resolve_output_paths(output_paths)
    periods = _required_array(velocity_outputs, schema.beat_period_seconds)
    calculator = LowRankWaveformDecompositionCalculator()
    outputs: dict[str, object] = {}

    sources = [
        (
            "artery/raw",
            schema.artery_per_beat.segment_velocity_signal,
        ),
        (
            "artery/bandlimited",
            schema.artery_per_beat.segment_velocity_signal_band_limited,
        ),
    ]
    if include_veins:
        sources.extend(
            (
                ("vein/raw", schema.vein_per_beat.segment_velocity_signal),
                (
                    "vein/bandlimited",
                    schema.vein_per_beat.segment_velocity_signal_band_limited,
                ),
            )
        )

    for source_name, dataset_path in sources:
        prefix = f"{schema.lowrank_waveform_decomposition_root}/{source_name}"
        if dataset_path is None or dataset_path not in velocity_outputs:
            outputs[f"{prefix}/qc/input_available"] = _metric(
                np.uint8(0), original_class="bool"
            )
            outputs[f"{prefix}/qc/missing_dataset_path"] = str(dataset_path)
            continue

        outputs[f"{prefix}/qc/input_available"] = _metric(
            np.uint8(1), original_class="bool"
        )
        waveforms = _metric_data(velocity_outputs[dataset_path])
        representation = calculator.compute(waveforms, periods)
        per_beat = calculator.compute_per_beat(waveforms, periods)
        _append_representation(
            outputs,
            prefix,
            source_name,
            str(dataset_path),
            representation,
            per_beat,
            calculator,
        )
    return outputs


def _append_representation(
    metrics: dict[str, object],
    prefix: str,
    source_name: str,
    dataset_path: str,
    rep: dict,
    per_beat: dict[str, np.ndarray],
    calculator: LowRankWaveformDecompositionCalculator,
) -> None:
    shape = rep["shape"]
    acq = rep["acq"]
    beatwise = rep["beatwise"]

    metrics[f"{prefix}/config/signal_source"] = source_name
    metrics[f"{prefix}/config/input_dataset_path"] = dataset_path
    metrics[f"{prefix}/config/svd_method"] = (
        "joint (sample, beat*branch*radius) plus per-beat "
        "(sample, branch*radius) SVD"
    )
    metrics[f"{prefix}/config/aggregation"] = (
        "acquisition median over all valid beat-location columns; "
        "beat endpoint median over valid locations"
    )
    metrics[f"{prefix}/config/max_exported_modes"] = calculator.exported_modes
    metrics[f"{prefix}/config/min_valid_samples_fraction"] = (
        calculator.min_valid_samples_fraction
    )
    metrics[f"{prefix}/config/min_valid_columns"] = calculator.min_valid_columns

    for name in (
        "n_t",
        "n_beats",
        "n_branches",
        "n_radii",
        "n_total_columns",
        "n_valid_columns",
    ):
        metrics[f"{prefix}/inputs/{name}"] = int(shape[name])
    metrics[f"{prefix}/inputs/valid_fraction_columns"] = (
        shape["n_valid_columns"] / float(max(1, shape["n_total_columns"]))
    )
    metrics[f"{prefix}/inputs/finite_fraction_per_column_bkr"] = _metric(
        rep["finite_fraction_per_column"],
        dims=("beat", "branch", "radius"),
    )
    metrics[f"{prefix}/inputs/valid_column_mask_bkr"] = _metric(
        rep["valid_column_mask"].astype(np.uint8),
        dims=("beat", "branch", "radius"),
        original_class="bool",
    )
    metrics[f"{prefix}/inputs/valid_columns_per_beat"] = _metric(
        rep["valid_counts_per_beat"], dims=("beat",)
    )
    metrics[f"{prefix}/inputs/valid_fraction_columns_per_beat"] = _metric(
        rep["valid_fraction_per_beat"], dims=("beat",)
    )
    metrics[f"{prefix}/inputs/beat_period_valid_b"] = _metric(
        rep["beat_period_valid"].astype(np.uint8),
        dims=("beat",),
        original_class="bool",
    )

    metrics[f"{prefix}/baseline/mu_bkr"] = _velocity_bkr(rep["mu"])
    metrics[f"{prefix}/baseline/mu_b"] = _velocity_beat(beatwise["mu_b"])
    metrics[f"{prefix}/baseline/abs_mu_b"] = _velocity_beat(
        beatwise["abs_mu_b"]
    )
    metrics[f"{prefix}/baseline/mu_acq"] = _metric(acq["mu_acq"], unit="mm/s")
    metrics[f"{prefix}/baseline/abs_mu_acq"] = _metric(
        acq["abs_mu_acq"], unit="mm/s"
    )
    metrics[f"{prefix}/baseline/sigma_mu_beat"] = _metric(
        acq["sigma_mu_beat"], unit="mm/s"
    )
    metrics[f"{prefix}/baseline/mad_mu_beat"] = _metric(
        acq["mad_mu_beat"], unit="mm/s"
    )

    for name in ("mean", "median", "std"):
        metrics[f"{prefix}/beat_period/{name}"] = _metric(
            acq[f"beat_period_{name}"], unit="s"
        )

    metrics[f"{prefix}/rms/total_pulsatile_rms_bkr"] = _velocity_bkr(
        rep["total_rms_bkr"]
    )
    metrics[f"{prefix}/rms/mean_to_pulsatile_ratio_bkr"] = _metric(
        rep["mean_to_pulsatile_ratio_bkr"],
        dims=("beat", "branch", "radius"),
        formula="abs(mu) / total_pulsatile_rms",
    )
    metrics[f"{prefix}/rms/mean_pulsatile_ratio_bkr"] = _metric(
        rep["mean_pulsatile_ratio_bkr"],
        dims=("beat", "branch", "radius"),
        formula="compatibility alias for mean_to_pulsatile_ratio_bkr",
    )
    for name, unit in (
        ("R0_b", "mm/s"),
        ("TPR_b", "mm/s"),
        ("rho0_b", None),
        ("MPR_b", None),
        ("mpr_b", None),
    ):
        metrics[f"{prefix}/beatwise/{name}"] = _metric(
            beatwise[name], unit=unit, dims=("beat",)
        )
    metrics[f"{prefix}/endpoints/R0"] = _metric(
        acq["R0"],
        unit="mm/s",
        formula="median_bkr(RMS_t(w))",
    )
    metrics[f"{prefix}/endpoints/TPR"] = _metric(
        acq["TPR"],
        unit="mm/s",
        formula="compatibility alias for R0",
    )
    metrics[f"{prefix}/endpoints/rho0"] = _metric(
        acq["rho0"], formula="R0 / R0 = 1"
    )
    metrics[f"{prefix}/endpoints/MPR"] = _metric(
        acq["MPR"], formula="median_bkr(abs(mu) / RMS_t(w))"
    )
    metrics[f"{prefix}/endpoints/mpr"] = _metric(
        acq["mpr"], formula="compatibility alias for MPR"
    )
    for name, unit in (
        ("sigma_R0_beat", "mm/s"),
        ("mad_R0_beat", "mm/s"),
        ("sigma_TPR_beat", "mm/s"),
        ("mad_TPR_beat", "mm/s"),
        ("sigma_rho0_beat", None),
        ("mad_rho0_beat", None),
        ("cv_rho0_beat", None),
        ("sigma_mpr_beat", None),
        ("mad_mpr_beat", None),
        ("cv_mpr_beat", None),
    ):
        metrics[f"{prefix}/variability/{name}"] = _metric(acq[name], unit=unit)

    for key, values in per_beat.items():
        unit = "mm/s" if key.startswith(("A", "R", "TPR")) else None
        metrics[f"{prefix}/per_beat/{key}"] = _metric(
            values, unit=unit, dims=("beat",)
        )

    available = bool(rep.get("svd_available", False))
    metrics[f"{prefix}/qc/svd_available"] = _metric(
        np.uint8(available), original_class="bool"
    )
    metrics[f"{prefix}/qc/svd_reason"] = str(rep.get("svd_reason", "unknown"))
    if not available:
        metrics[f"{prefix}/qc/n_modes"] = 0
        _append_unavailable_mode_outputs(metrics, prefix, shape, calculator)
        return

    metrics[f"{prefix}/qc/n_modes"] = int(rep["n_modes_panel"])
    metrics[f"{prefix}/qc/sign_flips_mode1to2"] = _metric(
        rep["sign_flips"], dims=("mode",)
    )
    for mode in range(1, calculator.exported_modes + 1):
        metrics[f"{prefix}/qc/denominator_floor_rho{mode}"] = _metric(
            np.uint8(not np.isfinite(acq.get(f"rho{mode}", np.nan))),
            original_class="bool",
        )

    metrics[f"{prefix}/decomposition/singular_values"] = _metric(
        rep["s"], unit="mm/s", dims=("mode",)
    )
    metrics[f"{prefix}/decomposition/singular_energy"] = _metric(
        rep["energy"], unit="(mm/s)^2", dims=("mode",)
    )
    metrics[f"{prefix}/decomposition/singular_energy_fraction"] = _metric(
        rep["energy_fraction"], dims=("mode",)
    )
    metrics[f"{prefix}/decomposition/spectrum_mode_count_M"] = np.asarray(
        acq["spectrum_mode_count"], dtype=int
    )
    for name in (
        "effective_rank",
        "participation_ratio",
        "alpha",
        "eta1",
        "eta2",
        "eta12",
    ):
        attrs = (
            {"formula": "sum(energy_2_to_M) / energy_1"}
            if name == "alpha"
            else {}
        )
        metrics[f"{prefix}/decomposition/{name}"] = _metric(
            acq[name], **attrs
        )
    metrics[f"{prefix}/decomposition/G1"] = _metric(
        acq["G1"],
        formula="(singular_value_1 - singular_value_2) / singular_value_1",
    )

    for mode in range(1, calculator.exported_modes + 1):
        _append_mode_outputs(metrics, prefix, mode, rep, acq)


def _append_mode_outputs(
    metrics: dict[str, object],
    prefix: str,
    mode: int,
    rep: dict,
    acq: dict,
) -> None:
    index = mode - 1
    metrics[f"{prefix}/decomposition/u_mode{mode}"] = _metric(
        rep["U_panel"][:, index], dims=("sample",)
    )
    metrics[f"{prefix}/decomposition/scores_mode{mode}_bkr"] = _velocity_bkr(
        rep["score_panel_bkr"][index]
    )
    metrics[f"{prefix}/rms/mode{mode}_amplitude_rms_bkr"] = _velocity_bkr(
        rep["rms_mode_panel"][index]
    )
    metrics[f"{prefix}/residuals/r{mode}_t_bkr"] = _metric(
        rep["residual_t_bkr_panel"][index],
        unit="mm/s",
        dims=("sample", "beat", "branch", "radius"),
    )
    metrics[f"{prefix}/residuals/rms_r{mode}_bkr"] = _velocity_bkr(
        rep["residual_rms_panel"][index]
    )
    for name, unit in (
        (f"A{mode}_b", "mm/s"),
        (f"R{mode}_b", "mm/s"),
        (f"rho{mode}_b", None),
        (f"median_abs_a{mode}_b", "mm/s"),
    ):
        metrics[f"{prefix}/beatwise/{name}"] = _metric(
            rep["beatwise"][name], unit=unit, dims=("beat",)
        )
    for name, unit in (
        (f"A{mode}", "mm/s"),
        (f"R{mode}", "mm/s"),
        (f"rho{mode}", None),
        (f"median_abs_a{mode}", "mm/s"),
    ):
        metrics[f"{prefix}/endpoints/{name}"] = _metric(acq[name], unit=unit)
    for stem, unit in (
        ("sigma_A", "mm/s"),
        ("mad_A", "mm/s"),
        ("cv_A", None),
        ("sigma_R", "mm/s"),
        ("mad_R", "mm/s"),
        ("cv_R", None),
        ("sigma_rho", None),
        ("mad_rho", None),
        ("cv_rho", None),
    ):
        name = f"{stem}{mode}_beat"
        metrics[f"{prefix}/variability/{name}"] = _metric(acq[name], unit=unit)
    for stem, unit in (("A", "mm/s"), ("R", "mm/s")):
        name = f"spatial_mad_{stem}{mode}_median_over_beats"
        metrics[f"{prefix}/variability/{name}"] = _metric(acq[name], unit=unit)


def _append_unavailable_mode_outputs(
    metrics: dict[str, object],
    prefix: str,
    shape: dict,
    calculator: LowRankWaveformDecompositionCalculator,
) -> None:
    n_t = shape["n_t"]
    bkr = (shape["n_beats"], shape["n_branches"], shape["n_radii"])
    for mode in range(1, calculator.exported_modes + 1):
        metrics[f"{prefix}/qc/denominator_floor_rho{mode}"] = _metric(
            np.uint8(1), original_class="bool"
        )
        metrics[f"{prefix}/decomposition/u_mode{mode}"] = _metric(
            np.full((n_t,), np.nan), dims=("sample",)
        )
        metrics[f"{prefix}/decomposition/scores_mode{mode}_bkr"] = _velocity_bkr(
            np.full(bkr, np.nan)
        )
        metrics[f"{prefix}/rms/mode{mode}_amplitude_rms_bkr"] = _velocity_bkr(
            np.full(bkr, np.nan)
        )
        metrics[f"{prefix}/residuals/r{mode}_t_bkr"] = _metric(
            np.full((n_t, *bkr), np.nan),
            unit="mm/s",
            dims=("sample", "beat", "branch", "radius"),
        )
        metrics[f"{prefix}/residuals/rms_r{mode}_bkr"] = _velocity_bkr(
            np.full(bkr, np.nan)
        )
        for name, unit in (
            (f"A{mode}_b", "mm/s"),
            (f"R{mode}_b", "mm/s"),
            (f"rho{mode}_b", None),
            (f"median_abs_a{mode}_b", "mm/s"),
        ):
            metrics[f"{prefix}/beatwise/{name}"] = _metric(
                np.full((bkr[0],), np.nan), unit=unit, dims=("beat",)
            )
        for name, unit in (
            (f"A{mode}", "mm/s"),
            (f"R{mode}", "mm/s"),
            (f"rho{mode}", None),
            (f"median_abs_a{mode}", "mm/s"),
        ):
            metrics[f"{prefix}/endpoints/{name}"] = _metric(np.nan, unit=unit)


def _velocity_bkr(data) -> DatasetValue:
    return _metric(data, unit="mm/s", dims=("beat", "branch", "radius"))


def _velocity_beat(data) -> DatasetValue:
    return _metric(data, unit="mm/s", dims=("beat",))


def _metric(
    data,
    *,
    unit: str | None = None,
    dims: tuple[str, ...] | None = None,
    formula: str | None = None,
    original_class: str | None = None,
):
    attrs: dict[str, object] = {}
    if unit:
        attrs["unit"] = unit
    if dims:
        attrs["dimDesc"] = list(dims)
    if formula:
        attrs["formula"] = formula
    if original_class:
        attrs["original_class"] = original_class
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
