"""Tests for the streaming Motion-JPEG AVI writer."""

from __future__ import annotations

from io import BytesIO
import json
import struct
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from PIL import Image

SRC_DIR = Path(__file__).resolve().parents[1] / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from input_output.output_manager import OutputType  # noqa: E402
from input_output.writers.avi import AviArtifactWriter  # noqa: E402


class AviWriterTests(unittest.TestCase):
    def test_mjpeg_avi_has_index_metadata_and_decodable_frames(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            output = _FakeOutput(Path(temp_dir))
            artifact_writer = AviArtifactWriter(output, "scan")
            metadata = {
                "title": "segment mosaic",
                "tiles": [{"label": "A3/R2"}],
            }
            with artifact_writer.open_mjpeg(
                "artery.avi",
                width=8,
                height=6,
                fps=60,
                subfolder="segment_velocity_map",
                metadata=metadata,
                jpeg_quality=100,
            ) as video:
                video.write_frame(np.zeros((6, 8, 3), dtype=np.uint8))
                video.write_frame(np.full((6, 8, 3), 180, dtype=np.uint8))

            self.assertEqual(
                Path(temp_dir)
                / "avi"
                / "segment_velocity_map"
                / "scan_artery.avi",
                video.path,
            )
            contents = video.path.read_bytes()
            self.assertEqual(b"RIFF", contents[:4])
            self.assertEqual(b"AVI ", contents[8:12])
            self.assertIn(b"MJPG", contents)
            self.assertIn(b"idx1", contents)
            self.assertEqual(2, _main_header_frame_count(contents))
            self.assertEqual(metadata, _comment_metadata(contents))

            frames = _jpeg_frames(contents)
            self.assertEqual(2, len(frames))
            self.assertEqual((6, 8, 3), frames[0].shape)
            self.assertLess(float(np.mean(frames[0])), 2.0)
            self.assertAlmostEqual(180.0, float(np.mean(frames[1])), delta=2.0)


def _main_header_frame_count(contents: bytes) -> int:
    position = contents.index(b"avih")
    payload_size = struct.unpack_from("<I", contents, position + 4)[0]
    fields = struct.unpack_from("<14I", contents, position + 8)
    if payload_size != 56:
        raise AssertionError(payload_size)
    return int(fields[4])


def _comment_metadata(contents: bytes) -> dict[str, object]:
    position = contents.index(b"ICMT")
    payload_size = struct.unpack_from("<I", contents, position + 4)[0]
    payload = contents[position + 8 : position + 8 + payload_size].rstrip(b"\x00")
    return json.loads(payload.decode("utf-8"))


def _jpeg_frames(contents: bytes) -> list[np.ndarray]:
    position = contents.index(b"movi") + 4
    index_position = contents.index(b"idx1", position)
    frames: list[np.ndarray] = []
    while position < index_position:
        chunk_id = contents[position : position + 4]
        chunk_size = struct.unpack_from("<I", contents, position + 4)[0]
        payload_start = position + 8
        payload = contents[payload_start : payload_start + chunk_size]
        if chunk_id == b"00dc":
            with Image.open(BytesIO(payload)) as image:
                frames.append(np.asarray(image.convert("RGB")))
        position = payload_start + chunk_size + chunk_size % 2
    return frames


class _FakeOutput:
    available = True

    def __init__(self, root: Path) -> None:
        self.root = root
        self.manager = SimpleNamespace(layout=SimpleNamespace(stem="scan"))

    def path_for(self, output_type: OutputType, filename: str | None = None) -> Path:
        return self.root / output_type.value / (filename or "scan")


if __name__ == "__main__":
    unittest.main()
