"""Tests for spatial topology persisted beside segment waveforms."""

from __future__ import annotations

import sys
import unittest
from types import SimpleNamespace
from unittest.mock import patch
from pathlib import Path

import h5py
import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from calculations.blood_flow_velocity import CrossSectionSignalSettings  # noqa: E402
from calculations.blood_flow_velocity.analysis_preparation.segments.generate_cross_section_signals import (  # noqa: E402
    _fill_cross_section_signals,
)
from input_output.schema import EyeFlowOutputPaths  # noqa: E402
from input_output.writers.h5 import write_value_dataset  # noqa: E402
from pipelines.waveform_shape_metrics.velocity.topology import (  # noqa: E402
    pack_segment_topology_outputs,
)


class SegmentCenterTests(unittest.TestCase):
    def test_cross_section_centroids_are_recorded_in_waveform_axis_order(self) -> None:
        labels = np.zeros((7, 7), dtype=np.int32)
        labels[2, 1] = 1
        labels[5, 4] = 2
        branches = SimpleNamespace(
            labels=labels,
            branch_ids=np.asarray([1, 2], dtype=np.int32),
        )
        sections = np.zeros((2, 7, 7), dtype=bool)
        sections[0, 2, 1] = True
        sections[1, 5, 4] = True
        velocity = np.zeros((3, 7, 7), dtype=np.float32)
        raw = np.full((2, 2, 3), np.nan, dtype=np.float32)
        safe = np.full_like(raw, np.nan)
        centers = np.full((2, 2, 2), np.nan, dtype=np.float32)
        settings = CrossSectionSignalSettings(3.0, False, 0.5, False, 0.01)

        signal = np.arange(3, dtype=np.float32)
        with patch(
            "calculations.blood_flow_velocity.analysis_preparation.segments."
            "generate_cross_section_signals._tilt_angle",
            return_value=0.0,
        ), patch(
            "calculations.blood_flow_velocity.analysis_preparation.segments."
            "generate_cross_section_signals._cross_section_velocity",
            return_value=(signal, signal),
        ):
            _fill_cross_section_signals(
                raw,
                safe,
                centers,
                velocity,
                sections,
                branches,
                np.asarray([3.0, 3.0], dtype=np.float32),
                settings,
            )

        np.testing.assert_array_equal(centers[0, 0], [1.0, 2.0])
        np.testing.assert_array_equal(centers[1, 1], [4.0, 5.0])
        self.assertTrue(np.all(np.isnan(centers[0, 1])))
        self.assertTrue(np.all(np.isnan(centers[1, 0])))
        np.testing.assert_array_equal(raw[0, 0], signal)
        np.testing.assert_array_equal(raw[1, 1], signal)


class TopologyOutputTests(unittest.TestCase):
    def test_topology_values_and_h5_attributes_preserve_index_contracts(self) -> None:
        optic_disc_mask = np.zeros((5, 7), dtype=bool)
        optic_disc_mask[2, 3] = True
        source = _source(optic_disc_mask=optic_disc_mask)
        artery = _segment_result(
            labels=np.asarray(
                [
                    [0, 5, 5, 0, 0, 0, 0],
                    [0, 5, 5, 0, 9, 9, 0],
                    [0, 0, 0, 0, 9, 9, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0],
                ],
                dtype=np.int32,
            ),
            branch_ids=[5, 9],
            centers=[[[1.5, 0.5], [2.0, 1.0]], [[4.5, 1.5], [np.nan, np.nan]]],
            radius_count=2,
        )
        vein = _segment_result(
            labels=np.zeros((5, 7), dtype=np.int32),
            branch_ids=[],
            centers=np.empty((0, 2, 2), dtype=np.float32),
            radius_count=2,
        )

        metrics = pack_segment_topology_outputs(source, artery, vein)
        schema = EyeFlowOutputPaths.active()

        with h5py.File("topology.h5", "w", driver="core", backing_store=False) as h5:
            for path, value in metrics.items():
                write_value_dataset(h5, path, value)

            np.testing.assert_array_equal(
                h5[schema.topology.artery.branch_label_map][()],
                artery.labels,
            )
            np.testing.assert_array_equal(
                h5[schema.topology.artery.branch_ids][()],
                [5, 9],
            )
            centers_ds = h5[schema.topology.artery.segment_center_xy]
            np.testing.assert_allclose(centers_ds[()], artery.segment_center_xy)
            self.assertEqual(
                list(centers_ds.attrs["dimDesc"]),
                ["branch", "radius", "coordinate"],
            )
            self.assertEqual(
                list(centers_ds.attrs["coordinate_order"]),
                ["x", "y"],
            )
            self.assertEqual(2, centers_ds.attrs["waveform_branch_axis"])
            self.assertEqual(3, centers_ds.attrs["waveform_radius_axis"])

            center_ds = h5[schema.topology.optic_disc.center_xy]
            np.testing.assert_array_equal(center_ds[()], [3.0, 2.0])
            self.assertEqual(list(center_ds.attrs["coordinate_order"]), ["x", "y"])

            mask_ds = h5[schema.topology.optic_disc.mask]
            np.testing.assert_array_equal(mask_ds[()].astype(bool), optic_disc_mask)
            self.assertEqual("bool", mask_ds.attrs["original_class"])
            self.assertEqual("dopplerview_segmentation", mask_ds.attrs["source"])

            vein_centers = h5[schema.topology.vein.segment_center_xy]
            self.assertEqual((0, 2, 2), vein_centers.shape)

    def test_optic_disc_mask_is_reconstructed_when_source_mask_is_missing(self) -> None:
        source = _source(optic_disc_mask=None)
        empty = _segment_result(
            labels=np.zeros((5, 7), dtype=np.int32),
            branch_ids=[],
            centers=np.empty((0, 2, 2), dtype=np.float32),
            radius_count=2,
        )

        metrics = pack_segment_topology_outputs(source, empty, empty)
        path = EyeFlowOutputPaths.active().topology.optic_disc.mask
        mask, attrs = metrics[path]

        self.assertTrue(mask[2, 3])
        self.assertEqual(
            "reconstructed_from_dopplerview_center_width_height",
            attrs["source"],
        )

    def test_topology_packing_rejects_branch_waveform_index_drift(self) -> None:
        source = _source(optic_disc_mask=None)
        segments = _segment_result(
            labels=np.asarray(
                [[1, 0, 2], [0, 0, 0]],
                dtype=np.int32,
            ),
            branch_ids=[1, 2],
            centers=[[[0.0, 0.0]], [[2.0, 0.0]]],
            radius_count=1,
        )
        segments.velocity = np.full((1, 1, 3), np.nan, dtype=np.float32)

        with self.assertRaisesRegex(ValueError, "segment waveform has 1 branches"):
            pack_segment_topology_outputs(source, segments, segments)


def _source(*, optic_disc_mask):
    return SimpleNamespace(
        retinal_artery_mask=np.zeros((5, 7), dtype=bool),
        optic_disc_mask=optic_disc_mask,
        optic_disc_center=np.asarray([3.0, 2.0], dtype=np.float32),
        optic_disc_width=np.float32(4.0),
        optic_disc_height=np.float32(2.0),
    )


def _segment_result(*, labels, branch_ids, centers, radius_count):
    ids = np.asarray(branch_ids, dtype=np.int32)
    return SimpleNamespace(
        labels=np.asarray(labels, dtype=np.int32),
        branch_ids=ids,
        segment_center_xy=np.asarray(centers, dtype=np.float32).reshape(
            ids.size,
            radius_count,
            2,
        ),
        velocity=np.full((radius_count, max(ids.size, 1), 3), np.nan, dtype=np.float32),
    )


if __name__ == "__main__":
    unittest.main()
