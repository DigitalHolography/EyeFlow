"""Frame-sequence and image-mask readers for input files."""

from __future__ import annotations

import math
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Iterator

import numpy as np
from .writers.h5 import normalize_h5_path

try:
    import cv2
except ImportError:
    cv2 = None

try:
    import h5py
except ImportError:
    h5py = None

class FrameSequence:
    def __init__(
        self,
        path: Path,
        h5_dataset: str,
        h5_frame_axis: int,
        h5_fps: float,
        h5_low_percentile: float,
        h5_high_percentile: float,
        *,
        h5_source: h5py.Dataset | None = None,
    ) -> None:
        self.path = path
        self.h5_source = h5_source
        self.is_h5 = h5_source is not None or path.suffix.lower() in {
            ".h5", ".hdf5", ".hdf"
        }
        self.h5_dataset = normalize_h5_path(h5_dataset)
        if h5_source is not None and normalize_h5_path(h5_source.name) != self.h5_dataset:
            raise ValueError(
                f"Supplied HDF5 dataset {h5_source.name!r} does not match "
                f"the requested path {self.h5_dataset!r}."
            )
        self.h5_frame_axis = int(h5_frame_axis)
        self.h5_low_percentile = float(h5_low_percentile)
        self.h5_high_percentile = float(h5_high_percentile)
        self.h5_value_min = 0.0
        self.h5_value_max = 1.0

        if self.is_h5:
            if h5py is None:
                raise RuntimeError("h5py est requis pour lire un fichier HDF5.")
            with self._open_h5_dataset() as dataset:
                if not isinstance(dataset, h5py.Dataset) or dataset.ndim != 3:
                    raise RuntimeError(
                        f"Le dataset {self.h5_dataset} doit Ãªtre 3-D, forme trouvÃ©e : "
                        f"{getattr(dataset, 'shape', None)}"
                    )
                shape = tuple(int(v) for v in dataset.shape)
            spatial_axes = [axis for axis in range(3) if axis != self.h5_frame_axis]
            self.frame_count = shape[self.h5_frame_axis]
            self.height = shape[spatial_axes[0]]
            self.width = shape[spatial_axes[1]]
            self.fps = float(h5_fps)
            self.h5_value_min, self.h5_value_max = self._estimate_h5_range()
        else:
            capture = cv2.VideoCapture(str(self.path))
            if not capture.isOpened():
                raise RuntimeError(f"Impossible d'ouvrir : {self.path}")
            self.frame_count = int(capture.get(cv2.CAP_PROP_FRAME_COUNT))
            self.width = int(capture.get(cv2.CAP_PROP_FRAME_WIDTH))
            self.height = int(capture.get(cv2.CAP_PROP_FRAME_HEIGHT))
            self.fps = float(capture.get(cv2.CAP_PROP_FPS))
            capture.release()
            if self.width <= 0 or self.height <= 0:
                raise RuntimeError("Dimensions vidÃ©o invalides.")
            if not math.isfinite(self.fps) or self.fps <= 0:
                self.fps = 25.0

    @contextmanager
    def _open_h5_dataset(self) -> Iterator[h5py.Dataset]:
        if self.h5_source is not None:
            yield self.h5_source
            return
        assert h5py is not None
        with h5py.File(self.path, "r") as handle:
            if self.h5_dataset not in handle:
                available = ", ".join(sorted(handle.keys()))
                raise RuntimeError(
                    f"Dataset HDF5 introuvable : {self.h5_dataset}. "
                    f"ClÃ©s disponibles : {available}"
                )
            yield handle[self.h5_dataset]

    def _slice_h5_frame(self, dataset: Any, index: int) -> np.ndarray:
        selector: list[Any] = [slice(None), slice(None), slice(None)]
        selector[self.h5_frame_axis] = int(index)
        frame = np.asarray(dataset[tuple(selector)], dtype=np.float32)
        frame = np.squeeze(frame)
        if frame.ndim != 2:
            raise RuntimeError(f"Une frame HDF5 doit Ãªtre 2-D, forme trouvÃ©e : {frame.shape}")
        return frame

    def _estimate_h5_range(self) -> tuple[float, float]:
        sample_count = min(max(self.frame_count, 1), 32)
        indices = np.linspace(0, max(0, self.frame_count - 1), sample_count).astype(int)
        samples: list[np.ndarray] = []
        with self._open_h5_dataset() as dataset:
            for index in np.unique(indices):
                frame = self._slice_h5_frame(dataset, int(index))
                values = frame[np.isfinite(frame)]
                if not values.size:
                    continue
                step = max(1, int(math.ceil(math.sqrt(values.size / 4096.0))))
                samples.append(values[::step].astype(np.float32, copy=False))
        if not samples:
            raise RuntimeError("Le dataset HDF5 ne contient aucune valeur finie.")
        values = np.concatenate(samples)
        low = float(np.percentile(values, self.h5_low_percentile))
        high = float(np.percentile(values, self.h5_high_percentile))
        if high <= low + 1e-12:
            low = float(values.min())
            high = float(values.max())
        if high <= low + 1e-12:
            high = low + 1.0
        return low, high

    def _h5_to_bgr_u8(self, frame: np.ndarray) -> np.ndarray:
        clean = np.nan_to_num(
            frame,
            nan=self.h5_value_min,
            posinf=self.h5_value_max,
            neginf=self.h5_value_min,
        )
        normalized = np.clip(
            (clean - self.h5_value_min) / (self.h5_value_max - self.h5_value_min),
            0.0,
            1.0,
        )
        gray = np.round(normalized * 255.0).astype(np.uint8)
        return cv2.cvtColor(gray, cv2.COLOR_GRAY2BGR)

    def iter_frames(self, max_frames: int | None = None) -> Iterator[np.ndarray]:
        limit = self.frame_count if max_frames is None else min(self.frame_count, max_frames)
        if self.is_h5:
            with self._open_h5_dataset() as dataset:
                for index in range(limit):
                    yield self._h5_to_bgr_u8(self._slice_h5_frame(dataset, index))
            return

        capture = cv2.VideoCapture(str(self.path))
        if not capture.isOpened():
            raise RuntimeError(f"Impossible d'ouvrir : {self.path}")
        try:
            count = 0
            while count < limit:
                ok, frame = capture.read()
                if not ok or frame is None:
                    break
                yield frame
                count += 1
        finally:
            capture.release()


def load_binary_mask(path: Path, shape: tuple[int, int], invert: bool) -> np.ndarray:
    mask = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
    if mask is None:
        raise RuntimeError(f"Impossible de lire le masque : {path}")
    if mask.shape != shape:
        mask = cv2.resize(mask, (shape[1], shape[0]), interpolation=cv2.INTER_NEAREST)
    mask = (mask > 0).astype(np.uint8) * 255
    return cv2.bitwise_not(mask) if invert else mask
