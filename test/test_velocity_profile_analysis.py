"""Weighted artery fits, observed-support sums, and persistent outputs."""

import sys
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pipeline_engine.context import PipelineH5Output
from pipelines.velocity_profile_analysis import fitting
from pipelines.velocity_profile_analysis.runner import (
    OUTPUT_ROOT,
    SOURCE_PATH,
    run_velocity_profile_analysis,
)


def analyze(y, **kwargs):
    return fitting.analyze_velocity_profiles(np.asarray(y)[:, None, None, None, None], **kwargs)


def scalar(result, name):
    return result[name].item()


def test_border_weights_are_quadratic_over_the_complete_domain():
    np.testing.assert_array_equal(
        fitting.border_weights(9),
        [0, 0.4375, 0.75, 0.9375, 1, 0.9375, 0.75, 0.4375, 0],
    )
    np.testing.assert_allclose(fitting.border_weights(6), [0, 0.64, 0.96, 0.96, 0.64, 0])
    np.testing.assert_array_equal(fitting.border_weights(1), [1])
    assert fitting.border_weights(0).size == 0


def test_border_weight_power_is_parameterized():
    np.testing.assert_array_equal(
        fitting.border_weights(5, power=1),
        [0, 0.5, 1, 0.5, 0],
    )
    np.testing.assert_array_equal(
        fitting.border_weights(5, power=4),
        [0, 0.9375, 1, 0.9375, 0],
    )


@pytest.mark.parametrize("power", [0, -1, np.nan, np.inf, -np.inf, True, "2", 1 + 0j])
def test_border_weight_power_must_be_a_finite_positive_real(power):
    with pytest.raises(ValueError, match="finite positive real"):
        fitting.border_weights(9, power=power)
    with pytest.raises(ValueError, match="finite positive real"):
        analyze(np.arange(9.0), weight_power=power)


def test_exact_quadratic_coefficients_vertex_roots_and_unweighted_sums():
    x = np.arange(11.0)
    y = -2 * (x - 1.25) * (x - 8.75)
    r = analyze(y)
    np.testing.assert_allclose([scalar(r, n) for n in ("a", "b", "c")], [-2, 20, -21.875])
    np.testing.assert_allclose(
        [scalar(r, n) for n in ("index_center", "index_left_zero", "index_right_zero")],
        [5, 1.25, 8.75],
    )
    support = (x >= 1.25) & (x <= 8.75)
    assert scalar(r, "Qv") == pytest.approx(y[support].sum())
    assert scalar(r, "Qv_fit") == pytest.approx(y[support].sum())
    assert scalar(r, "Qv_fit") != pytest.approx(
        (y[support] * fitting.border_weights(11)[support]).sum()
    )
    assert scalar(r, "n_fit_samples") == 11
    assert scalar(r, "n_area_samples") == 7
    assert scalar(r, "fit_r_squared") == pytest.approx(1)
    assert scalar(r, "fit_weighted_r_squared") == pytest.approx(1)
    assert scalar(r, "fit_rmse") < 1e-12


@pytest.mark.parametrize("power", [1.0, fitting.DEFAULT_WEIGHT_POWER, 4.0])
def test_weighted_fit_matches_explicit_reference_and_quality_formulas(power):
    x = np.arange(17.0)
    y = 10 - 0.3 * (x - 7) ** 2 + np.sin(x)
    y[[0, 7, 15]] = [np.nan, np.inf, np.nan]
    finite = np.isfinite(y)
    w = fitting.border_weights(len(x), power=power)[finite]
    matrix = np.column_stack((x[finite] ** 2, x[finite], np.ones(finite.sum())))
    expected = np.linalg.lstsq(matrix * np.sqrt(w)[:, None], y[finite] * np.sqrt(w), rcond=None)[0]
    r = analyze(y, weight_power=power)
    np.testing.assert_allclose([scalar(r, n) for n in ("a", "b", "c")], expected, rtol=1e-6)
    residual = y[finite] - matrix @ expected
    rss, wrss = np.sum(residual**2), np.sum(w * residual**2)
    expected_metrics = {
        "fit_rss": rss,
        "fit_rmse": np.sqrt(rss / finite.sum()),
        "fit_weighted_rss": wrss,
        "fit_weighted_rmse": np.sqrt(wrss / w.sum()),
        "fit_r_squared": 1 - rss / np.sum((y[finite] - y[finite].mean()) ** 2),
        "fit_weighted_r_squared": 1
        - wrss / np.sum(w * (y[finite] - np.average(y[finite], weights=w)) ** 2),
    }
    for name, value in expected_metrics.items():
        assert scalar(r, name) == pytest.approx(value, rel=1e-6)
    unweighted = np.linalg.lstsq(matrix, y[finite], rcond=None)[0]
    assert not np.allclose(expected, unweighted)
    left, right = np.sort(np.roots(expected))
    support = finite & (x >= left) & (x <= right)
    assert scalar(r, "Qv") == pytest.approx(np.sum(y[support]), rel=1e-6)
    assert scalar(r, "Qv_fit") == pytest.approx(np.polyval(expected, x[support]).sum(), rel=1e-6)


@pytest.mark.parametrize("kind", ["upward", "flat", "linear", "nearlinear"])
def test_unsuitable_curvature_retains_fit_but_not_geometry(kind):
    x = np.arange(9.0)
    y = {
        "upward": (x - 4) ** 2 - 4,
        "flat": np.full(9, 3.0),
        "linear": x + 1,
        "nearlinear": x + 1 - 1e-16 * x * x,
    }[kind]
    r = analyze(y)
    assert all(np.isfinite(scalar(r, n)) for n in ("a", "b", "c", "fit_rmse"))
    assert all(
        np.isnan(scalar(r, n))
        for n in ("index_center", "index_left_zero", "index_right_zero", "Qv", "Qv_fit")
    )
    assert scalar(r, "n_area_samples") == 0
    if kind == "flat":
        assert np.isnan(scalar(r, "fit_r_squared"))
        assert np.isnan(scalar(r, "fit_weighted_r_squared"))


def test_downward_fit_without_real_roots_retains_vertex():
    r = analyze(-((np.arange(9.0) - 4) ** 2) - 3)
    assert scalar(r, "index_center") == pytest.approx(4)
    for name in ("index_left_zero", "index_right_zero", "Qv", "Qv_fit"):
        assert np.isnan(scalar(r, name))


def test_repeated_root_is_not_mistaken_for_two_distinct_zeros():
    r = analyze(-((np.arange(9.0) - 4) ** 2))
    assert scalar(r, "index_center") == pytest.approx(4)
    assert np.isnan(scalar(r, "index_left_zero"))
    assert np.isnan(scalar(r, "Qv_fit"))


def test_observed_area_preserves_negative_sign():
    y = np.array(
        [
            -2.92171186,
            -1.03488268,
            -2.51977378,
            1.52562512,
            14.71491973,
            -25.66658441,
            -2.36850265,
            1.76512421,
            2.9599399,
        ]
    )
    r = analyze(y)
    x = np.arange(len(y))
    support = (x >= scalar(r, "index_left_zero")) & (x <= scalar(r, "index_right_zero"))
    assert scalar(r, "Qv") < 0
    assert scalar(r, "Qv") == pytest.approx(y[support].sum(), rel=1e-6)
    assert scalar(r, "Qv_fit") > 0


@pytest.mark.parametrize(
    "y,count", [([np.nan] * 5, 0), ([np.nan, 1, np.inf, 2, np.nan], 2), ([1], 1), ([], 0)]
)
def test_insufficient_samples(y, count):
    r = analyze(y)
    assert scalar(r, "n_fit_samples") == count
    assert all(np.isnan(scalar(r, name)) for name in fitting.FLOAT_OUTPUTS)


def test_roots_outside_domain_and_missing_samples_use_only_observed_support():
    x = np.arange(9.0)
    y = -(x + 2.25) * (x - 12.25)
    y[[2, 6]] = np.nan
    r = analyze(y)
    np.testing.assert_allclose(
        [scalar(r, n) for n in ("index_left_zero", "index_right_zero")], [-2.25, 12.25]
    )
    assert scalar(r, "n_area_samples") == 7
    assert scalar(r, "Qv") == pytest.approx(np.nansum(y))
    assert scalar(r, "Qv_fit") == pytest.approx(np.nansum(y))


def test_no_observed_integer_between_roots_returns_nan_areas():
    r = analyze(-(np.arange(9.0) - 0.2) * (np.arange(9.0) - 0.8))
    assert scalar(r, "n_area_samples") == 0
    assert np.isfinite(scalar(r, "index_left_zero"))
    assert np.isnan(scalar(r, "Qv"))
    assert np.isnan(scalar(r, "Qv_fit"))


def test_axis_order_and_time_blocks():
    shape = (2, 2, 3, 2)
    factors = np.arange(1, np.prod(shape) + 1).reshape(shape)
    x = np.arange(9.0)
    v = (-(x - 1.25) * (x - 6.75))[:, None, None, None, None] * factors[None]
    actual = fitting.analyze_velocity_profiles(v, time_block_size=1)
    np.testing.assert_allclose(actual["a"], -factors)
    np.testing.assert_allclose(actual["b"], 8 * factors)
    np.testing.assert_allclose(actual["c"], -8.4375 * factors)
    for value in actual.values():
        assert value.shape == shape


@pytest.mark.parametrize("axis", [1, 2, 3, 4])
def test_empty_axes_preserved(axis):
    shape = [9, 2, 2, 2, 2]
    shape[axis] = 0
    r = fitting.analyze_velocity_profiles(np.empty(shape))
    assert all(value.shape == tuple(shape[1:]) and value.size == 0 for value in r.values())


def test_shared_validity_uses_one_solve_for_multiple_time_profiles():
    x = np.arange(9.0)
    v = np.repeat((-(x - 1) * (x - 7))[:, None, None, None, None], 4, axis=1)
    v[2, :2] = np.nan
    with patch.object(np.linalg, "lstsq", wraps=np.linalg.lstsq) as solve:
        fitting.analyze_velocity_profiles(v)
    assert solve.call_count == 2


def test_rank_deficient_solver_result_returns_nan_with_sample_count():
    with patch.object(np.linalg, "lstsq", return_value=(np.ones((3, 1)), [], 2, [])):
        r = analyze(np.arange(9.0))
    assert scalar(r, "n_fit_samples") == 9
    assert np.isnan(scalar(r, "a"))


def test_hdf5_round_trip_and_bounded_reads(tmp_path):
    v = np.broadcast_to(
        (10 - (np.arange(9.0) - 4) ** 2)[:, None, None, None, None], (9, 3, 1, 2, 1)
    )
    path = tmp_path / "analysis.h5"
    with h5py.File(path, "w") as h5:
        h5.create_dataset(SOURCE_PATH, data=v)
        ctx = SimpleNamespace(output=SimpleNamespace(h5=PipelineH5Output(h5)))
        outputs = run_velocity_profile_analysis(ctx)
        ctx.output.h5.write_many(outputs)
    with h5py.File(path, "r") as h5:
        group = h5[OUTPUT_ROOT]
        assert set(group) == set(fitting.FLOAT_OUTPUTS + fitting.COUNT_OUTPUTS)
        for name in group:
            ds = group[name]["value"]
            assert ds.shape == (3, 1, 2, 1)
            assert ds.dtype == (np.int32 if name in fitting.COUNT_OUTPUTS else np.float32)
            assert list(ds.attrs["dimDesc"]) == ["time", "beat", "branch", "radius"]
            assert ds.attrs["index_base"] == 0
            assert ds.attrs["source_path"] == SOURCE_PATH
            assert ds.attrs["weight_power"] == fitting.DEFAULT_WEIGHT_POWER
            assert ds.attrs["weight_definition"] == (
                "u=x/(Nx-1); d=abs(2*u-1); w=1-d^p"
            )
        assert "/Processing/VelocityProfileAnalysis/Vein" not in h5
        np.testing.assert_array_equal(h5[SOURCE_PATH], v)


def test_no_full_dataset_materialization():
    class SlabOnly:
        shape = (9, 7, 1, 1, 1)
        dtype = np.dtype("float32")

        def __getitem__(self, key):
            assert key[1].stop - key[1].start <= 2
            return np.ones((9, key[1].stop - key[1].start))

        def __array__(self, *args, **kwargs):
            raise AssertionError("full dataset read")

    assert fitting.analyze_velocity_profiles(SlabOnly(), time_block_size=2)["a"].shape == (
        7,
        1,
        1,
        1,
    )


def test_missing_dataset_and_invalid_shape(tmp_path):
    with h5py.File(tmp_path / "missing.h5", "w") as h5:
        ctx = SimpleNamespace(output=SimpleNamespace(h5=PipelineH5Output(h5)))
        with pytest.raises(KeyError, match="Required velocity-profile dataset"):
            run_velocity_profile_analysis(ctx)
        h5.create_dataset(SOURCE_PATH, data=np.ones((4, 3)))
        with pytest.raises(ValueError, match="shape"):
            run_velocity_profile_analysis(ctx)
