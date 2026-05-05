from __future__ import annotations

from dataclasses import dataclass

from domain.blood_flow_velocity import PerBeatAnalysisInput


@dataclass(frozen=True)
class WaveformShapeMetricsContext:
    per_beat_analysis: PerBeatAnalysisInput
    dopplerview_analysis: dict[str, object]
    attrs: dict[str, object]


@dataclass
class DopplerViewStepContext:
    cache: dict[str, object]
    holodoppler_config: dict[str, object]
    dopplerview_config: dict[str, object]

    def require(self, key: str):
        if key not in self.cache:
            raise RuntimeError(f"Missing required context key: '{key}'")
        return self.cache[key]

    def set(self, key: str, value) -> None:
        self.cache[key] = value
