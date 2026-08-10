"""Tests for the native low-rank waveform decomposition pipeline."""

from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import pipelines  # noqa: E402
from input_output.schema import EyeFlowOutputPaths  # noqa: E402
from input_output.writers.h5 import write_value_dataset  # noqa: E402
from pipeline_engine import PIPELINE_REGISTRY, DatasetValue, PipelineDAG  # noqa: E402
from pipelines.lowrank_waveform_decomposition.calculator import (  # noqa: E402
    LowRankWaveformDecompositionCalculator,
)
from pipelines.lowrank_waveform_decomposition.outputs import (  # noqa: E402
    pack_lowrank_waveform_decomposition_outputs,
)
from pipelines.lowrank_waveform_decomposition.runner import (  # noqa: E402
    LOWRANK_WAVEFORM_OUTPUTS_STATE,
    run_lowrank_waveform_decomposition,
)
from pipelines.waveform_velocity_core import runner as core_runner  # noqa: E402


class _State:
    def __init__(self, values=None):
        self.values = dict(values or {})

    def get(self, key, default=None):
        return self.values.get(key, default)

    def set(self, key, value):
        self.values[key] = value


def _data(value):
    return value.data if isinstance(value, DatasetValue) else np.asarray(value)


def _synthetic_rank_two_block() -> tuple[np.ndarray, np.ndarray]:
    sample = np.linspace(0.0, 2.0 * np.pi, 48, endpoint=False)
    mode1 = np.sin(sample)
    mode2 = np.cos(2.0 * sample)
    coefficients1 = np.asarray(
        [
            [[1.0, 1.2], [0.8, 1.1]],
            [[1.1, 1.3], [0.9, 1.2]],
            [[0.9, 1.1], [0.7, 1.0]],
            [[1.2, 1.4], [1.0, 1.3]],
        ]
    )
    coefficients2 = np.asarray(
        [
            [[0.30, -0.20], [0.10, -0.25]],
            [[0.25, -0.15], [0.15, -0.20]],
            [[0.35, -0.25], [0.05, -0.30]],
            [[0.20, -0.10], [0.20, -0.15]],
        ]
    )
    baseline = 3.0 + 0.1 * np.arange(4)[:, None, None]
    block = (
        baseline[None]
        + mode1[:, None, None, None] * coefficients1[None]
        + mode2[:, None, None, None] * coefficients2[None]
    )
    return block.astype(np.float32), np.asarray([[0.8, 0.82, 0.79, 0.81]])


class LowRankWaveformDecompositionTests(unittest.TestCase):
    def test_pipeline_is_discoverable_and_depends_on_waveform_velocity(self) -> None:
        pipelines.load_pipeline_catalog()
        descriptor = PIPELINE_REGISTRY["lowrank_waveform_decomposition"]
        self.assertEqual(("waveform_velocity",), descriptor.dag_requires)
        self.assertEqual(("veins",), tuple(option.name for option in descriptor.options))
        self.assertFalse(descriptor.options[0].default_enabled)
        self.assertEqual(
            "Processing/Metrics/lowrank_waveform_decomposition",
            EyeFlowOutputPaths.active().lowrank_waveform_decomposition_root,
        )

        plan = PipelineDAG(PIPELINE_REGISTRY.values()).resolve_targets(
            ["lowrank_waveform_decomposition"]
        )
        self.assertEqual(
            (
                "waveform_velocity_core",
                "waveform_velocity",
                "lowrank_waveform_decomposition",
            ),
            plan.names,
        )

    def test_joint_svd_recovers_rank_two_panel_and_endpoints(self) -> None:
        block, periods = _synthetic_rank_two_block()
        result = LowRankWaveformDecompositionCalculator().compute(block, periods)

        self.assertTrue(result["svd_available"])
        self.assertEqual((48, 2), result["U_panel"].shape)
        self.assertEqual((2, 4, 2, 2), result["score_panel_bkr"].shape)
        self.assertEqual((2, 48, 4, 2, 2), result["residual_t_bkr_panel"].shape)
        self.assertAlmostEqual(float(np.sum(result["energy_fraction"][:2])), 1.0)
        self.assertGreater(result["acq"]["A1"], result["acq"]["A2"])
        self.assertLess(result["acq"]["rho2"], 1e-6)
        self.assertGreater(result["acq"]["rho1"], result["acq"]["rho2"])
        np.testing.assert_allclose(result["beatwise"]["mu_b"], [3.0, 3.1, 3.2, 3.3])

        valid = result["valid_column_mask"]
        expected_r0 = float(np.median(result["total_rms_bkr"][valid]))
        expected_mpr = float(
            np.median(result["mean_to_pulsatile_ratio_bkr"][valid])
        )
        expected_a1 = float(np.median(result["rms_mode_panel"][0][valid]))
        singular_values = result["s"]
        expected_alpha = float(
            np.sum(singular_values[1:] ** 2) / singular_values[0] ** 2
        )
        expected_g1 = float(
            (singular_values[0] - singular_values[1]) / singular_values[0]
        )
        self.assertAlmostEqual(result["acq"]["R0"], expected_r0)
        self.assertAlmostEqual(result["acq"]["TPR"], expected_r0)
        self.assertAlmostEqual(result["acq"]["rho0"], 1.0)
        self.assertAlmostEqual(result["acq"]["MPR"], expected_mpr)
        self.assertAlmostEqual(result["acq"]["A1"], expected_a1)
        self.assertAlmostEqual(result["acq"]["alpha"], expected_alpha)
        self.assertAlmostEqual(result["acq"]["G1"], expected_g1)
        self.assertEqual(result["acq"]["spectrum_mode_count"], result["s"].size)

    def test_per_beat_svd_exports_one_endpoint_per_beat(self) -> None:
        block, periods = _synthetic_rank_two_block()
        result = LowRankWaveformDecompositionCalculator().compute_per_beat(
            block, periods
        )

        self.assertEqual((4,), result["A1_b_pb"].shape)
        self.assertTrue(np.all(np.isfinite(result["A1_b_pb"])))
        self.assertTrue(np.all(result["rho2_b_pb"] < 1e-6))
        self.assertTrue(np.all(result["alpha_b_pb"] > 0))
        self.assertTrue(np.all(result["G1_b_pb"] > 0))
        np.testing.assert_allclose(result["rho0_b_pb"], 1.0)
        np.testing.assert_allclose(result["R0_b_pb"], result["TPR_b_pb"])

    def test_acquisition_uses_one_median_over_all_valid_columns(self) -> None:
        sample = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
        amplitudes = np.asarray([[[1.0, 2.0, 3.0]], [[100.0, 1.0, 1.0]]])
        block = 5.0 + np.sin(sample)[:, None, None, None] * amplitudes[None]
        block[:, 1, 0, 1:] = np.nan

        result = LowRankWaveformDecompositionCalculator().compute(
            block, [0.8, 0.8]
        )

        valid_rms = result["total_rms_bkr"][result["valid_column_mask"]]
        expected = float(np.median(valid_rms))
        two_stage = float(np.median(result["beatwise"]["R0_b"]))
        self.assertAlmostEqual(result["acq"]["R0"], expected)
        self.assertNotAlmostEqual(result["acq"]["R0"], two_stage)

    def test_outputs_use_schema_root_and_veins_are_opt_in(self) -> None:
        schema = EyeFlowOutputPaths.active()
        block, periods = _synthetic_rank_two_block()
        inputs = {
            schema.beat_period_seconds: periods,
            schema.artery_per_beat.segment_velocity_signal: block,
            schema.artery_per_beat.segment_velocity_signal_band_limited: block * 0.9,
            schema.vein_per_beat.segment_velocity_signal: block * 0.5,
            schema.vein_per_beat.segment_velocity_signal_band_limited: block * 0.4,
        }

        artery_only = pack_lowrank_waveform_decomposition_outputs(inputs)
        root = schema.lowrank_waveform_decomposition_root
        endpoint = f"{root}/artery/raw/endpoints/A1"
        self.assertIn(endpoint, artery_only)
        self.assertIn(f"{root}/artery/raw/endpoints/R0", artery_only)
        self.assertIn(f"{root}/artery/raw/endpoints/MPR", artery_only)
        self.assertNotIn(f"{root}/artery/raw/endpoints/mpr_prime", artery_only)
        self.assertGreater(float(_data(artery_only[endpoint])), 0)
        self.assertFalse(any(key.startswith(f"{root}/vein/") for key in artery_only))

        with_veins = pack_lowrank_waveform_decomposition_outputs(
            inputs, include_veins=True
        )
        self.assertIn(f"{root}/vein/raw/endpoints/A1", with_veins)
        self.assertIn(f"{root}/artery/raw/per_beat/A1_b_pb", with_veins)

    def test_packed_outputs_write_to_h5_with_dimensions_and_formulas(self) -> None:
        schema = EyeFlowOutputPaths.active()
        block, periods = _synthetic_rank_two_block()
        outputs = pack_lowrank_waveform_decomposition_outputs(
            {
                schema.beat_period_seconds: periods,
                schema.artery_per_beat.segment_velocity_signal: block,
                schema.artery_per_beat.segment_velocity_signal_band_limited: block,
            }
        )
        root = schema.lowrank_waveform_decomposition_root

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "lowrank.h5"
            with h5py.File(path, "w") as h5file:
                for key, value in outputs.items():
                    write_value_dataset(h5file, key, value)
            with h5py.File(path, "r") as h5file:
                residual = h5file[f"{root}/artery/raw/residuals/r2_t_bkr"]
                alpha = h5file[f"{root}/artery/raw/decomposition/alpha"]
                g1 = h5file[f"{root}/artery/raw/decomposition/G1"]
                self.assertEqual((48, 4, 2, 2), residual.shape)
                self.assertEqual(
                    ["sample", "beat", "branch", "radius"],
                    [str(value) for value in residual.attrs["dimDesc"]],
                )
                self.assertEqual(
                    "sum(energy_2_to_M) / energy_1",
                    alpha.attrs["formula"],
                )
                self.assertEqual(
                    "(singular_value_1 - singular_value_2) / singular_value_1",
                    g1.attrs["formula"],
                )
                self.assertNotIn("unit", g1.attrs)

    def test_too_few_columns_keeps_baseline_and_marks_svd_unavailable(self) -> None:
        sample = np.linspace(0.0, 2.0 * np.pi, 32, endpoint=False)
        block = 2.0 + np.sin(sample)[:, None, None, None] * np.ones((1, 1, 1, 2))
        result = LowRankWaveformDecompositionCalculator().compute(block, [0.8])

        self.assertFalse(result["svd_available"])
        self.assertEqual("too_few_valid_columns", result["svd_reason"])
        self.assertAlmostEqual(result["acq"]["mu_acq"], 2.0)
        self.assertGreater(result["acq"]["R0"], 0)
        self.assertEqual(result["acq"]["R0"], result["acq"]["TPR"])

    def test_runner_consumes_shared_state_and_defaults_to_arteries(self) -> None:
        schema = EyeFlowOutputPaths.active()
        block, periods = _synthetic_rank_two_block()
        state = _State(
            {
                core_runner.VELOCITY_PER_BEAT_OUTPUTS_STATE: {
                    schema.beat_period_seconds: periods,
                    schema.artery_per_beat.segment_velocity_signal: block,
                    schema.artery_per_beat.segment_velocity_signal_band_limited: block,
                }
            }
        )
        ctx = SimpleNamespace(
            state=state,
            options_for=lambda name: frozenset(),
        )

        outputs = run_lowrank_waveform_decomposition(ctx)

        self.assertTrue(outputs)
        self.assertIs(outputs, state.get(LOWRANK_WAVEFORM_OUTPUTS_STATE))
        self.assertFalse(any("/vein/" in key for key in outputs))

    def test_lowrank_schedule_forces_per_beat_segment_foundation(self) -> None:
        ctx = SimpleNamespace(
            pipeline_scheduled=lambda name: name
            == "lowrank_waveform_decomposition",
            options_for=lambda name: frozenset(),
        )

        self.assertTrue(core_runner._per_beat_required(ctx))
        self.assertTrue(core_runner._segments_required(ctx))


if __name__ == "__main__":
    unittest.main()
