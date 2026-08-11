"""Temporal filtering for dense displacement fields."""

import numpy as np


def componentwise_median(fields: list[np.ndarray]) -> np.ndarray:
    return np.median(np.stack(fields, axis=0), axis=0).astype(np.float32)


class CenteredMedianBuffer:
    def __init__(self, window: int) -> None:
        self.radius = window // 2
        self.fields: list[np.ndarray] = []
        self.start_index = 0
        self.seen = 0
        self.next_output = 0

    def _field_at(self, index: int) -> np.ndarray:
        return self.fields[index - self.start_index]

    def _emit(self, center: int, total: int) -> np.ndarray:
        indices = [
            min(max(index, 0), total - 1)
            for index in range(center - self.radius, center + self.radius + 1)
        ]
        median = componentwise_median([self._field_at(index) for index in indices])
        self.next_output += 1
        keep_from = max(0, self.next_output - self.radius)
        drop_count = keep_from - self.start_index
        if drop_count > 0:
            del self.fields[:drop_count]
            self.start_index = keep_from
        return median

    def push(self, field: np.ndarray) -> list[np.ndarray]:
        self.fields.append(field)
        self.seen += 1
        outputs: list[np.ndarray] = []
        while self.next_output + self.radius < self.seen:
            outputs.append(self._emit(self.next_output, self.seen))
        return outputs

    def finish(self) -> list[np.ndarray]:
        outputs: list[np.ndarray] = []
        while self.next_output < self.seen:
            outputs.append(self._emit(self.next_output, self.seen))
        return outputs
