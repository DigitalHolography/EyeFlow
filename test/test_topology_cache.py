"""Tests for run-scoped prepared-topology reuse."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.topology import (  # noqa: E402
    PreparedTopology,
    AnnulusGeometry,
    OpticDisc,
    prepare_topologies,
    run_topology_cache,
    topology_source_id,
)


class TopologyCacheTests(unittest.TestCase):
    def test_caches_final_joint_windows_and_logs_hits(self) -> None:
        from dataclasses import replace
        from utils.logger import Logger
        vessel = np.ones((9, 9), bool)
        disc = np.zeros_like(vessel)
        settings = AnnulusGeometry(0., 1., 1., 1)
        initial = [_prepared_topology(vessel.shape, marker=float(i)) for i in (1, 2)]
        for prepared in initial:
            prepared.topology.window_side_pixels = 1
        final = [
            replace(prepared, topology=SimpleNamespace(
                **{**vars(prepared.topology), "window_side_pixels": 3},
            ))
            for prepared in initial
        ]
        cache = {}
        with patch("calculations.topology.workflow.prepare_topology", side_effect=initial) as prepare, patch(
            "calculations.topology.workflow._resize_prepared_topology", side_effect=final,
        ) as resize, patch.object(Logger, "log") as log:
            first = prepare_topologies(
                {"artery": vessel, "vein": vessel}, _disc(disc), settings, source_id="scan", cache=cache,
            )
            second = prepare_topologies(
                {"artery": vessel, "vein": vessel}, _disc(disc), settings, source_id="scan", cache=cache,
            )
        self.assertEqual(2, prepare.call_count)
        self.assertEqual(2, resize.call_count)
        for name in first:
            self.assertIs(first[name], second[name])
        messages = [call.args[0] for call in log.call_args_list]
        self.assertTrue(any("Topology cache hit: artery" in message for message in messages))
        self.assertTrue(any("Topology cache hit: vein" in message for message in messages))

    def test_cache_is_owned_by_one_run_state(self) -> None:
        first_state: dict[str, object] = {}
        second_state: dict[str, object] = {}

        first = run_topology_cache(first_state)

        self.assertIs(first, run_topology_cache(first_state))
        self.assertIsNot(first, run_topology_cache(second_state))

    def test_identical_inputs_reuse_topology_and_changed_mask_invalidates_it(
        self,
    ) -> None:
        vessel_mask = np.zeros((9, 9), dtype=bool)
        vessel_mask[3:6, 2:7] = True
        changed_mask = vessel_mask.copy()
        changed_mask[2, 4] = True
        optic_disc_mask = np.zeros_like(vessel_mask)
        optic_disc_mask[4, 4] = True
        settings = AnnulusGeometry(0.0, 1.0, 1.0, 1)
        first_prepared = _prepared_topology(vessel_mask.shape, marker=1.0)
        changed_prepared = _prepared_topology(vessel_mask.shape, marker=2.0)
        cache = {}

        with patch(
            "calculations.topology.workflow.prepare_topology",
            side_effect=(first_prepared, changed_prepared),
        ) as prepare:
            first = prepare_topologies(
                {"artery": vessel_mask},
                _disc(optic_disc_mask),
                settings,
                source_id="scan-a",
                cache=cache,
            )
            reused = prepare_topologies(
                {"artery": vessel_mask.copy()},
                _disc(optic_disc_mask.copy()),
                settings,
                source_id="scan-a",
                cache=cache,
            )
            changed = prepare_topologies(
                {"artery": changed_mask},
                _disc(optic_disc_mask),
                settings,
                source_id="scan-a",
                cache=cache,
            )

        self.assertIs(first["artery"], reused["artery"])
        self.assertIs(changed_prepared, changed["artery"])
        self.assertEqual(2, prepare.call_count)

    def test_source_identity_includes_both_input_files(self) -> None:
        self.assertEqual(
            "scan_HD.h5|scan_DV.h5",
            topology_source_id("scan_HD.h5", "scan_DV.h5"),
        )

    def test_changed_other_vessel_mask_does_not_invalidate_topology(self) -> None:
        artery = np.zeros((9, 9), dtype=bool)
        artery[3:6, 2:4] = True
        vein = np.zeros_like(artery)
        vein[3:6, 6] = True
        changed_vein = vein.copy()
        changed_vein[2, 6] = True
        disc = np.zeros_like(artery)
        disc[4, 4] = True
        settings = AnnulusGeometry(0.0, 1.0, 1.0, 1)
        prepared = [
            _prepared_topology(artery.shape, marker=float(index))
            for index in range(3)
        ]
        cache = {}

        with patch(
            "calculations.topology.workflow.prepare_topology",
            side_effect=prepared,
        ) as prepare:
            first = prepare_topologies(
                {"artery": artery, "vein": vein},
                _disc(disc),
                settings,
                source_id="scan-a",
                cache=cache,
                window_side_pixels=3,
            )
            changed = prepare_topologies(
                {"artery": artery, "vein": changed_vein},
                _disc(disc),
                settings,
                source_id="scan-a",
                cache=cache,
                window_side_pixels=3,
            )

        self.assertIs(first["artery"], changed["artery"])
        self.assertIsNot(first["vein"], changed["vein"])
        self.assertEqual(3, prepare.call_count)

    def test_changed_explicit_center_invalidates_cached_topology(self) -> None:
        vessel = np.zeros((9, 9), dtype=bool)
        vessel[3:6, 2:7] = True
        disc = np.zeros_like(vessel)
        settings = AnnulusGeometry(0.0, 1.0, 1.0, 1)
        prepared = [
            _prepared_topology(vessel.shape, marker=float(index))
            for index in range(2)
        ]
        cache = {}
        with patch(
            "calculations.topology.workflow.prepare_topology",
            side_effect=prepared,
        ) as prepare:
            first = prepare_topologies(
                {"artery": vessel}, _disc(disc, (3.0, 4.0)), settings,
                source_id="scan", cache=cache,
            )
            second = prepare_topologies(
                {"artery": vessel}, _disc(disc, (4.0, 4.0)), settings,
                source_id="scan", cache=cache,
            )
        self.assertIsNot(first["artery"], second["artery"])
        self.assertEqual(2, prepare.call_count)


def _prepared_topology(
    shape: tuple[int, int],
    *,
    marker: float,
) -> PreparedTopology:
    labels = np.zeros(shape, dtype=np.int32)
    labels[3:6, 3:6] = 1
    topology = SimpleNamespace(
        window_side_pixels=3,
        annulus_masks=np.ones((1, *shape), dtype=bool),
        branch_ids=np.asarray([1], dtype=np.int32),
        labels=labels,
    )
    return PreparedTopology(
        topology=topology,
        rotation_degrees=np.asarray([[marker]], dtype=np.float32),
        interpolated_masks=np.ones((1, 1, 3, 3), dtype=bool),
        rotated_masks=np.ones((1, 1, 5, 5), dtype=bool),
    )


def _disc(mask: np.ndarray, center=(4.0, 4.0)) -> OpticDisc:
    return OpticDisc(mask, center, None, None)


if __name__ == "__main__":
    unittest.main()
