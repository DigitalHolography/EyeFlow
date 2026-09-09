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
    SegmentRingSettings,
    prepare_topologies,
    run_topology_cache,
    topology_source_id,
)


class TopologyCacheTests(unittest.TestCase):
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
        settings = SegmentRingSettings(0.0, 1.0, 1.0, 1)
        first_prepared = _prepared_topology(vessel_mask.shape, marker=1.0)
        changed_prepared = _prepared_topology(vessel_mask.shape, marker=2.0)
        cache = {}

        with patch(
            "calculations.topology.workflow.prepare_topology",
            side_effect=(first_prepared, changed_prepared),
        ) as prepare:
            first = prepare_topologies(
                {"artery": vessel_mask},
                optic_disc_mask,
                settings,
                source_id="scan-a",
                cache=cache,
            )
            reused = prepare_topologies(
                {"artery": vessel_mask.copy()},
                optic_disc_mask.copy(),
                settings,
                source_id="scan-a",
                cache=cache,
            )
            changed = prepare_topologies(
                {"artery": changed_mask},
                optic_disc_mask,
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


if __name__ == "__main__":
    unittest.main()
