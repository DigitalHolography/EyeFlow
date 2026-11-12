"""Tests for per-segment velocity-map mosaic AVI export."""

from __future__ import annotations

import json
import struct
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.output_manager import OutputType  # noqa: E402
from input_output.schema import EyeFlowOutputPaths  # noqa: E402
from pipeline_engine.base import DatasetValue  # noqa: E402
from pipelines.velocity_analysis.artifacts.segment_map_video import (  # noqa: E402
    export_segment_velocity_map_avis,
)


class SegmentVelocityMapAviTests(unittest.TestCase):
    def test_two_movies_share_global_range_and_embed_tile_labels(self) -> None:
        schema = EyeFlowOutputPaths.active()
        artery_maps = np.full((2, 2, 2, 3, 1, 1), 2.0, dtype=np.float32)
        vein_maps = np.full((2, 2, 2, 3, 1, 1), 8.0, dtype=np.float32)
        packed = {
            schema.artery_segments.velocity_map_per_segment: DatasetValue(artery_maps),
            schema.vein_segments.velocity_map_per_segment: DatasetValue(vein_maps),
        }
        artery = SimpleNamespace(branch_ids=np.asarray([3], dtype=np.int32))
        vein = SimpleNamespace(branch_ids=np.asarray([7], dtype=np.int32))

        with tempfile.TemporaryDirectory() as temp_dir:
            output = _FakeOutput(Path(temp_dir))
            paths = export_segment_velocity_map_avis(
                output,
                artery,
                vein,
                packed,
            )

            self.assertEqual(2, len(paths))
            metadata = [_comment_metadata(Path(path).read_bytes()) for path in paths]
            self.assertEqual([[2.0, 8.0], [2.0, 8.0]], [m["velocity_range"] for m in metadata])
            self.assertEqual(["A3/R1"], [tile["label"] for tile in metadata[0]["tiles"]])
            self.assertEqual(["V7/R1"], [tile["label"] for tile in metadata[1]["tiles"]])
            self.assertTrue(all(m["fps"] == 60.0 for m in metadata))
            self.assertTrue(
                all(
                    m["velocity_range_source"] == "temporal_median_image"
                    for m in metadata
                )
            )
            self.assertEqual([3, 3], [m["source_beat_count"] for m in metadata])
            self.assertEqual([1, 1], [m["exported_beat_count"] for m in metadata])
            self.assertEqual([[0], [0]], [m["exported_beat_indexes"] for m in metadata])
            self.assertEqual([2, 2], [_avi_frame_count(Path(path)) for path in paths])
            self.assertEqual(
                "/ProcessingMoments/VelocityMapPerSegment/Artery",
                metadata[0]["velocity_map_dataset"],
            )
            self.assertEqual(
                "/Segmentation/Vein/Segments",
                metadata[1]["segment_dataset"],
            )
            self.assertTrue(
                all("segment_velocity_map" in Path(path).parts for path in paths)
            )


def _comment_metadata(contents: bytes) -> dict[str, object]:
    position = contents.index(b"ICMT")
    payload_size = struct.unpack_from("<I", contents, position + 4)[0]
    payload = contents[position + 8 : position + 8 + payload_size].rstrip(b"\x00")
    return json.loads(payload.decode("utf-8"))


def _avi_frame_count(path: Path) -> int:
    contents = path.read_bytes()
    position = contents.index(b"avih")
    fields = struct.unpack_from("<14I", contents, position + 8)
    return int(fields[4])


class _FakeOutput:
    available = True

    def __init__(self, root: Path) -> None:
        self.root = root
        self.manager = SimpleNamespace(layout=SimpleNamespace(stem="scan"))

    def path_for(self, output_type: OutputType, filename: str | None = None) -> Path:
        return self.root / output_type.value / (filename or "scan")


if __name__ == "__main__":
    unittest.main()
