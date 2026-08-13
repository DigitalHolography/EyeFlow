"""Streaming Motion-JPEG AVI artifact writer."""

from __future__ import annotations

from io import BytesIO
import json
from pathlib import Path
import struct
from typing import BinaryIO, Mapping

import numpy as np
from PIL import Image


_AVI_HAS_INDEX = 0x10
_AVI_KEY_FRAME = 0x10
_UINT32_MAX = (1 << 32) - 1


class AviArtifactWriter:
    """Create stem-prefixed AVI artifacts in an output namespace."""

    def __init__(self, output, stem: str | None = None) -> None:
        self.output = output
        self.stem = str(stem) if stem else _output_stem(output)

    def path(self, suffix: str, *, subfolder: str | None = None) -> Path:
        filename = f"{self.stem}_{suffix}"
        if subfolder:
            filename = f"{subfolder}/{filename}"
        path = self.output.path_for(_avi_output_type(), filename)
        path.parent.mkdir(parents=True, exist_ok=True)
        return path

    def open_mjpeg(
        self,
        suffix: str,
        *,
        width: int,
        height: int,
        fps: float,
        subfolder: str | None = None,
        metadata: Mapping[str, object] | None = None,
        jpeg_quality: int = 82,
    ) -> MjpegAviWriter:
        """Open a streaming MJPEG writer for one AVI artifact."""
        return MjpegAviWriter(
            self.path(suffix, subfolder=subfolder),
            width=width,
            height=height,
            fps=fps,
            metadata=metadata,
            jpeg_quality=jpeg_quality,
        )


class MjpegAviWriter:
    """Write RGB frames to an AVI 1.0 file using the MJPEG codec."""

    def __init__(
        self,
        path: str | Path,
        *,
        width: int,
        height: int,
        fps: float,
        metadata: Mapping[str, object] | None = None,
        jpeg_quality: int = 82,
    ) -> None:
        if width <= 0 or height <= 0:
            raise ValueError("AVI width and height must be positive.")
        if width > 32767 or height > 32767:
            raise ValueError("AVI width and height must fit the 16-bit stream rectangle.")
        if fps <= 0:
            raise ValueError("AVI frame rate must be positive.")
        if not 1 <= jpeg_quality <= 100:
            raise ValueError("JPEG quality must be between 1 and 100.")
        self.path = Path(path)
        self.width = int(width)
        self.height = int(height)
        self.fps = float(fps)
        self.metadata = dict(metadata or {})
        self.jpeg_quality = int(jpeg_quality)
        self._file: BinaryIO | None = None
        self._frame_index: list[tuple[int, int]] = []
        self._maximum_frame_size = 0
        self._hdrl_position = 12
        self._hdrl_size = 0
        self._movi_size_position = 0
        self._movi_type_position = 0

    def __enter__(self) -> MjpegAviWriter:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._file = self.path.open("w+b")
        self._write_container_header()
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        if exc_type is not None:
            self._discard_partial_file()
            return False
        try:
            self.close()
        except Exception:
            self._discard_partial_file()
            raise
        return False

    @property
    def frame_count(self) -> int:
        return len(self._frame_index)

    def write_frame(self, frame) -> None:
        """JPEG-encode and append one RGB frame."""
        file = self._required_file()
        image = _rgb_uint8(frame)
        if image.shape[:2] != (self.height, self.width):
            raise ValueError(
                "AVI frame shape does not match the configured canvas: "
                f"expected {(self.height, self.width, 3)}, got {image.shape}."
            )
        encoded = _jpeg_bytes(image, quality=self.jpeg_quality)
        chunk_position = file.tell()
        relative_offset = chunk_position - self._movi_type_position
        if relative_offset > _UINT32_MAX:
            raise ValueError("AVI 1.0 output exceeded its 4 GiB offset limit.")
        file.write(_chunk(b"00dc", encoded))
        self._frame_index.append((relative_offset, len(encoded)))
        self._maximum_frame_size = max(self._maximum_frame_size, len(encoded))

    def close(self) -> Path:
        """Finalize the AVI index and patch its frame-count headers."""
        file = self._required_file()
        movi_end = file.tell()
        index_payload = b"".join(
            struct.pack(
                "<4sIII",
                b"00dc",
                _AVI_KEY_FRAME,
                offset,
                size,
            )
            for offset, size in self._frame_index
        )
        file.write(_chunk(b"idx1", index_payload))
        file_end = file.tell()
        riff_size = file_end - 8
        movi_size = movi_end - self._movi_size_position - 4
        if max(riff_size, movi_size) > _UINT32_MAX:
            raise ValueError("AVI 1.0 output exceeded its 4 GiB size limit.")

        file.seek(4)
        file.write(struct.pack("<I", riff_size))
        file.seek(self._movi_size_position)
        file.write(struct.pack("<I", movi_size))
        final_hdrl = self._header_list(
            total_frames=self.frame_count,
            suggested_buffer_size=self._maximum_frame_size,
        )
        if len(final_hdrl) != self._hdrl_size:
            raise RuntimeError("AVI header size changed while finalizing output.")
        file.seek(self._hdrl_position)
        file.write(final_hdrl)
        file.flush()
        file.close()
        self._file = None
        return self.path

    def _write_container_header(self) -> None:
        file = self._required_file()
        file.write(b"RIFF\x00\x00\x00\x00AVI ")
        header = self._header_list(total_frames=0, suggested_buffer_size=0)
        self._hdrl_size = len(header)
        file.write(header)
        file.write(_info_list(self.metadata))
        file.write(b"LIST")
        self._movi_size_position = file.tell()
        file.write(b"\x00\x00\x00\x00")
        self._movi_type_position = file.tell()
        file.write(b"movi")

    def _header_list(
        self,
        *,
        total_frames: int,
        suggested_buffer_size: int,
    ) -> bytes:
        rate_scale = 1000
        rate = max(1, int(round(self.fps * rate_scale)))
        microseconds_per_frame = max(1, int(round(1_000_000.0 / self.fps)))
        maximum_bytes_per_second = int(round(suggested_buffer_size * self.fps))
        main_header = struct.pack(
            "<14I",
            microseconds_per_frame,
            maximum_bytes_per_second,
            0,
            _AVI_HAS_INDEX,
            total_frames,
            0,
            1,
            suggested_buffer_size,
            self.width,
            self.height,
            0,
            0,
            0,
            0,
        )
        stream_header = struct.pack(
            "<4s4sIHHIIIIIIIIhhhh",
            b"vids",
            b"MJPG",
            0,
            0,
            0,
            0,
            rate_scale,
            rate,
            0,
            total_frames,
            suggested_buffer_size,
            _UINT32_MAX,
            0,
            0,
            0,
            self.width,
            self.height,
        )
        stream_format = struct.pack(
            "<IiiHH4sIiiII",
            40,
            self.width,
            self.height,
            1,
            24,
            b"MJPG",
            self.width * self.height * 3,
            0,
            0,
            0,
            0,
        )
        stream_list = _list_chunk(
            b"strl",
            _chunk(b"strh", stream_header) + _chunk(b"strf", stream_format),
        )
        return _list_chunk(b"hdrl", _chunk(b"avih", main_header) + stream_list)

    def _required_file(self) -> BinaryIO:
        if self._file is None:
            raise RuntimeError("MJPEG AVI writer is not open.")
        return self._file

    def _discard_partial_file(self) -> None:
        if self._file is not None:
            self._file.close()
            self._file = None
        self.path.unlink(missing_ok=True)


def _chunk(fourcc: bytes, payload: bytes) -> bytes:
    padding = b"\x00" if len(payload) % 2 else b""
    return fourcc + struct.pack("<I", len(payload)) + payload + padding


def _list_chunk(list_type: bytes, payload: bytes) -> bytes:
    return b"LIST" + struct.pack("<I", len(payload) + 4) + list_type + payload


def _info_list(metadata: Mapping[str, object]) -> bytes:
    if not metadata:
        return b""
    comment = json.dumps(metadata, ensure_ascii=True, separators=(",", ":"))
    payload = _text_chunk(b"INAM", str(metadata.get("title", "EyeFlow AVI")))
    payload += _text_chunk(b"ISFT", "EyeFlow")
    payload += _text_chunk(b"ICMT", comment)
    return _list_chunk(b"INFO", payload)


def _text_chunk(fourcc: bytes, value: str) -> bytes:
    return _chunk(fourcc, value.encode("utf-8") + b"\x00")


def _jpeg_bytes(image: np.ndarray, *, quality: int) -> bytes:
    output = BytesIO()
    Image.fromarray(image, mode="RGB").save(
        output,
        format="JPEG",
        quality=quality,
        subsampling=2,
    )
    return output.getvalue()


def _rgb_uint8(frame) -> np.ndarray:
    image = np.asarray(frame)
    if image.ndim == 2:
        image = np.repeat(image[:, :, None], 3, axis=2)
    if image.ndim != 3 or image.shape[2] != 3:
        raise ValueError("AVI frames must be two-dimensional grayscale or RGB arrays.")
    if image.dtype != np.uint8:
        image = np.clip(image, 0, 255).astype(np.uint8)
    return np.ascontiguousarray(image)


def _output_stem(output) -> str:
    manager = getattr(output, "manager", None)
    layout = getattr(manager, "layout", None)
    stem = getattr(layout, "stem", None)
    if stem is None:
        stem = getattr(getattr(output, "layout", None), "stem", None)
    return str(stem or "eyeflow")


def _avi_output_type():
    from input_output.output_manager import OutputType

    return OutputType.AVI
