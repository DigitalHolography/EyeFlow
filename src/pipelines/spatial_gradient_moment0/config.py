"""Typed configuration for spatial-gradient preprocessing experiments.

The editable JSON file next to this module is the single source of experimental
parameters.  Dataclass defaults deliberately mirror the pre-experiment EyeFlow
pipeline so an installed package still has a safe fallback if the JSON file is
unavailable.
"""

from __future__ import annotations

import copy
import json
import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping


CONFIG_ENVIRONMENT_VARIABLE = "EYEFLOW_SPATIAL_GRADIENT_CONFIG"
DEFAULT_CONFIG_PATH = Path(__file__).with_name("preprocessing_config.json")


@dataclass
class CropConfig:
    enabled: bool = True
    source: str = "velocity_geometry"


@dataclass
class RotationConfig:
    enabled: bool = True
    interpolation: str = "existing_bilinear"


@dataclass
class Gaussian3DConfig:
    enabled: bool = True
    sigma_x: float = 1.0
    sigma_y: float = 1.0
    sigma_t: float = 2.0
    truncate: float = 4.0


@dataclass
class TemporalFilterConfig:
    enabled: bool = True
    method: str = "median"
    window: int = 17
    sigma: float = 2.0
    truncate: float = 4.0


@dataclass
class SpatialFilterConfig:
    enabled: bool = False
    method: str = "gaussian"
    kernel_x: int = 3
    kernel_y: int = 7
    sigma_x: float = 1.0
    sigma_y: float = 3.0
    truncate: float = 4.0


@dataclass
class GradientConfig:
    enabled: bool = True
    method: str = "sobel"
    direction: str = "magnitude"
    signed: bool = False
    derivative_sigma_t: float = 0.0
    derivative_sigma_x: float = 1.0
    derivative_sigma_y: float = 1.0
    truncate: float = 4.0


@dataclass
class ClaheConfig:
    enabled: bool = False
    clip_limit: float = 0.01
    tile_x: int = 8
    tile_y: int = 8


@dataclass
class DogConfig:
    enabled: bool = False
    sigma_small_x: float = 1.0
    sigma_small_y: float = 1.0
    sigma_small_t: float = 0.0
    sigma_large_x: float = 4.0
    sigma_large_y: float = 4.0
    sigma_large_t: float = 0.0
    truncate: float = 4.0


@dataclass
class NonLocalMeansConfig:
    enabled: bool = False
    patch_size: int = 5
    patch_distance: int = 6
    h: float = 0.8
    fast_mode: bool = True


@dataclass
class BilateralConfig:
    enabled: bool = False
    sigma_color: float = 0.05
    sigma_spatial: float = 3.0
    window_size: int = 7


@dataclass
class LongitudinalFilterConfig:
    enabled: bool = False
    method: str = "mean"
    kernel_y: int = 7
    trim_fraction: float = 0.1


@dataclass
class PostGradientFilterConfig:
    enabled: bool = True
    method: str = "gaussian"
    sigma_x: float = 0.0
    sigma_y: float = 0.0
    sigma_t: float = 2.0
    kernel_x: int = 1
    kernel_y: int = 1
    kernel_t: int = 1
    truncate: float = 4.0


@dataclass
class OutputConfig:
    save_debug_avi: bool = True
    save_intermediates: bool = False
    save_final_float_tiff: bool = False
    save_previews: bool = True
    diagnostic_frames: list[int] = field(default_factory=lambda: [0])
    directory_name: str = "spatial_gradient_experiments"


@dataclass
class PreprocessingConfig:
    """One fully resolved preprocessing experiment."""

    name: str = "legacy_existing"
    pipeline: list[str] = field(
        default_factory=lambda: [
            "crop",
            "rotate",
            "temporal_filter",
            "gradient",
            "post_gradient_filter",
        ]
    )
    crop: CropConfig = field(default_factory=CropConfig)
    rotation: RotationConfig = field(default_factory=RotationConfig)
    gaussian3d: Gaussian3DConfig = field(default_factory=Gaussian3DConfig)
    temporal_filter: TemporalFilterConfig = field(default_factory=TemporalFilterConfig)
    spatial_filter: SpatialFilterConfig = field(default_factory=SpatialFilterConfig)
    gradient: GradientConfig = field(default_factory=GradientConfig)
    clahe: ClaheConfig = field(default_factory=ClaheConfig)
    dog: DogConfig = field(default_factory=DogConfig)
    non_local_means: NonLocalMeansConfig = field(default_factory=NonLocalMeansConfig)
    bilateral: BilateralConfig = field(default_factory=BilateralConfig)
    longitudinal_filter: LongitudinalFilterConfig = field(
        default_factory=LongitudinalFilterConfig
    )
    post_gradient_filter: PostGradientFilterConfig = field(
        default_factory=PostGradientFilterConfig
    )
    output: OutputConfig = field(default_factory=OutputConfig)

    @classmethod
    def from_dict(cls, values: Mapping[str, Any]) -> "PreprocessingConfig":
        _reject_unknown_keys(values, _PREPROCESSING_FIELDS, "preprocessing")
        kwargs: dict[str, Any] = {}
        for key, config_type in _NESTED_CONFIG_TYPES.items():
            section = values.get(key, {})
            if not isinstance(section, Mapping):
                raise ValueError(f"preprocessing.{key} must be an object.")
            kwargs[key] = _dataclass_from_mapping(
                config_type, section, path=f"preprocessing.{key}"
            )
        if "name" in values:
            kwargs["name"] = str(values["name"])
        if "pipeline" in values:
            pipeline = values["pipeline"]
            if not isinstance(pipeline, list) or not all(
                isinstance(stage, str) for stage in pipeline
            ):
                raise ValueError("preprocessing.pipeline must be a list of stage names.")
            kwargs["pipeline"] = list(pipeline)
        return cls(**kwargs)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class SweepConfig:
    enabled: bool = False
    presets: list[str] = field(default_factory=list)
    parameters: dict[str, list[Any]] = field(default_factory=dict)

    @classmethod
    def from_dict(cls, values: Mapping[str, Any]) -> "SweepConfig":
        _reject_unknown_keys(values, {"enabled", "presets", "parameters"}, "sweep")
        presets = values.get("presets", [])
        parameters = values.get("parameters", {})
        if not isinstance(presets, list) or not all(isinstance(item, str) for item in presets):
            raise ValueError("sweep.presets must be a list of preset names.")
        if not isinstance(parameters, Mapping):
            raise ValueError("sweep.parameters must be an object of parameter lists.")
        normalized_parameters: dict[str, list[Any]] = {}
        for path, candidates in parameters.items():
            if not isinstance(candidates, list):
                raise ValueError(f"Sweep parameter {path!r} must contain a list.")
            if not candidates:
                raise ValueError(f"Sweep parameter {path!r} has no candidate values.")
            normalized_parameters[str(path)] = copy.deepcopy(candidates)
        return cls(
            enabled=bool(values.get("enabled", False)),
            presets=list(presets),
            parameters=normalized_parameters,
        )


@dataclass
class ConfigurationBundle:
    path: Path
    active_preset: str
    preprocessing: PreprocessingConfig
    sweep: SweepConfig
    base: dict[str, Any]
    presets: dict[str, dict[str, Any]]
    overrides: dict[str, Any]

    def resolve_preset(self, name: str) -> PreprocessingConfig:
        if name not in self.presets:
            choices = ", ".join(sorted(self.presets))
            raise ValueError(f"Unknown preprocessing preset {name!r}; choose from: {choices}.")
        merged = _deep_merge(self.base, self.presets[name])
        merged = _deep_merge(merged, self.overrides)
        merged["name"] = name
        return PreprocessingConfig.from_dict(merged)


_NESTED_CONFIG_TYPES = {
    "crop": CropConfig,
    "rotation": RotationConfig,
    "gaussian3d": Gaussian3DConfig,
    "temporal_filter": TemporalFilterConfig,
    "spatial_filter": SpatialFilterConfig,
    "gradient": GradientConfig,
    "clahe": ClaheConfig,
    "dog": DogConfig,
    "non_local_means": NonLocalMeansConfig,
    "bilateral": BilateralConfig,
    "longitudinal_filter": LongitudinalFilterConfig,
    "post_gradient_filter": PostGradientFilterConfig,
    "output": OutputConfig,
}
_PREPROCESSING_FIELDS = {"name", "pipeline", *_NESTED_CONFIG_TYPES}


def load_configuration(path: str | Path | None = None) -> ConfigurationBundle:
    """Load the central JSON file and resolve its active named preset."""

    config_path = _resolve_config_path(path)
    if not config_path.is_file():
        if path is not None or os.environ.get(CONFIG_ENVIRONMENT_VARIABLE):
            raise FileNotFoundError(f"Spatial-gradient configuration not found: {config_path}")
        default = PreprocessingConfig()
        return ConfigurationBundle(
            path=config_path,
            active_preset=default.name,
            preprocessing=default,
            sweep=SweepConfig(),
            base=default.to_dict(),
            presets={default.name: {}},
            overrides={},
        )
    try:
        payload = json.loads(config_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid JSON in {config_path}: {exc}") from exc
    if not isinstance(payload, Mapping):
        raise ValueError(f"Spatial-gradient configuration must be a JSON object: {config_path}")
    _reject_unknown_keys(
        payload,
        {"schema_version", "active_preset", "preprocessing", "presets", "overrides", "sweep"},
        "configuration",
    )
    if payload.get("schema_version", 1) != 1:
        raise ValueError("Unsupported spatial-gradient configuration schema_version.")
    base = payload.get("preprocessing", {})
    presets = payload.get("presets", {})
    overrides = payload.get("overrides", {})
    if not all(isinstance(value, Mapping) for value in (base, presets, overrides)):
        raise ValueError("preprocessing, presets, and overrides must be JSON objects.")
    normalized_presets: dict[str, dict[str, Any]] = {}
    for name, preset in presets.items():
        if not isinstance(preset, Mapping):
            raise ValueError(f"Preset {name!r} must be a JSON object.")
        normalized_presets[str(name)] = copy.deepcopy(dict(preset))
    active_preset = str(payload.get("active_preset", "legacy_existing"))
    if active_preset not in normalized_presets:
        choices = ", ".join(sorted(normalized_presets))
        raise ValueError(f"Unknown active_preset {active_preset!r}; choose from: {choices}.")
    merged = _deep_merge(dict(base), normalized_presets[active_preset])
    merged = _deep_merge(merged, dict(overrides))
    merged["name"] = active_preset
    return ConfigurationBundle(
        path=config_path,
        active_preset=active_preset,
        preprocessing=PreprocessingConfig.from_dict(merged),
        sweep=SweepConfig.from_dict(payload.get("sweep", {})),
        base=copy.deepcopy(dict(base)),
        presets=normalized_presets,
        overrides=copy.deepcopy(dict(overrides)),
    )


def config_with_parameter(
    config: PreprocessingConfig,
    path: str,
    value: Any,
) -> PreprocessingConfig:
    """Return a config copy with one dotted path replaced (used by sweeps)."""

    payload = config.to_dict()
    parts = path.split(".")
    if not parts or any(not part for part in parts):
        raise ValueError(f"Invalid sweep parameter path: {path!r}.")
    cursor: Any = payload
    for part in parts[:-1]:
        if not isinstance(cursor, dict) or part not in cursor:
            raise ValueError(f"Unknown sweep parameter path: {path!r}.")
        cursor = cursor[part]
    if not isinstance(cursor, dict) or parts[-1] not in cursor:
        raise ValueError(f"Unknown sweep parameter path: {path!r}.")
    cursor[parts[-1]] = copy.deepcopy(value)
    return PreprocessingConfig.from_dict(payload)


def _resolve_config_path(path: str | Path | None) -> Path:
    selected = path or os.environ.get(CONFIG_ENVIRONMENT_VARIABLE) or DEFAULT_CONFIG_PATH
    return Path(selected).expanduser().resolve()


def _deep_merge(base: Mapping[str, Any], overlay: Mapping[str, Any]) -> dict[str, Any]:
    merged = copy.deepcopy(dict(base))
    for key, value in overlay.items():
        if isinstance(value, Mapping) and isinstance(merged.get(key), Mapping):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def _dataclass_from_mapping(config_type, values: Mapping[str, Any], *, path: str):
    known = set(config_type.__dataclass_fields__)
    _reject_unknown_keys(values, known, path)
    return config_type(**dict(values))


def _reject_unknown_keys(values: Mapping[str, Any], known: set[str], path: str) -> None:
    unknown = sorted(str(key) for key in values if key not in known)
    if unknown:
        raise ValueError(f"Unknown {path} setting(s): {', '.join(unknown)}.")


__all__ = [
    "CONFIG_ENVIRONMENT_VARIABLE",
    "DEFAULT_CONFIG_PATH",
    "ConfigurationBundle",
    "PreprocessingConfig",
    "SweepConfig",
    "config_with_parameter",
    "load_configuration",
]
