"""One representative AVI output contract."""

from io import BytesIO
import json
import struct

import numpy as np
from PIL import Image

from input_output.output_manager import OutputManager
from input_output.runs.layout import HoloRunLayout
from input_output.writers.avi import AviArtifactWriter


def test_mjpeg_avi_has_metadata_index_and_decodable_frames(tmp_path):
    manager = OutputManager(HoloRunLayout.from_holo(tmp_path / "scan.holo"))
    metadata = {"title": "segment mosaic", "tiles": [{"label": "A3/R2"}]}
    with AviArtifactWriter(manager).open_mjpeg(
        "artery.avi", width=8, height=6, fps=60,
        subfolder="segment_velocity_map", metadata=metadata, jpeg_quality=100,
    ) as video:
        video.write_frame(np.zeros((6, 8, 3), dtype=np.uint8))
        video.write_frame(np.full((6, 8, 3), 180, dtype=np.uint8))

    assert video.path == manager.layout.ef_dir / "avi/segment_velocity_map/scan_artery.avi"
    contents = video.path.read_bytes()
    assert contents[:4] == b"RIFF"
    assert contents[8:12] == b"AVI "
    assert b"MJPG" in contents
    assert b"idx1" in contents
    header = contents.index(b"avih")
    assert struct.unpack_from("<I", contents, header + 8 + 4 * 4)[0] == 2

    comment = contents.index(b"ICMT")
    comment_size = struct.unpack_from("<I", contents, comment + 4)[0]
    payload = contents[comment + 8:comment + 8 + comment_size].rstrip(b"\x00")
    assert json.loads(payload.decode("utf-8")) == metadata

    position = contents.index(b"movi") + 4
    frames = []
    while contents[position:position + 4] == b"00dc":
        size = struct.unpack_from("<I", contents, position + 4)[0]
        with Image.open(BytesIO(contents[position + 8:position + 8 + size])) as image:
            frames.append(np.asarray(image.convert("RGB")))
        position += 8 + size + size % 2
    assert len(frames) == 2
    assert frames[0].shape == (6, 8, 3)
    assert float(np.mean(frames[0])) < 2.0
    assert abs(float(np.mean(frames[1])) - 180.0) < 2.0
