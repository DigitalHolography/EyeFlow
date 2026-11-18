"""Application identity, resource discovery, and persisted settings."""

from __future__ import annotations

import json
import os
import re
import sys
from collections.abc import Iterable, Iterator, Mapping, MutableMapping
from copy import deepcopy
from importlib import metadata as importlib_metadata
from numbers import Integral
from pathlib import Path
from typing import TYPE_CHECKING, Any, Generic, Literal, TypeVar, cast

if TYPE_CHECKING:
    from pipeline_engine.planning.base import PipelineDescriptor

from velocity_calibration import (
    DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
    validate_band_ratio_frequency_scale_hz,
)

APP_NAME = "EyeFlow"
SETTINGS_FILENAME = "settings.json"
DEFAULT_SETTINGS_FILENAME = "default_settings.json"
# Application fallback; calculations own their independent API default.
DEFAULT_BAND_LIMITED_HARMONICS_SETTING = 12
DEFAULT_COPY_HD_PASSTHROUGH = False
DEFAULT_NUMBER_OF_RADII_IN_HALF_DIAGONAL = 32
LAST_RUN_LOG_FILENAME = "last_EF_log.txt"
PIPELINES_DIR_ENV = "EYEFLOW_PIPELINES_DIR"
VERSION_PATTERN = re.compile(r'^version\s*=\s*"([^"]+)"\s*$')
INVALID_PATH_CHARS_PATTERN = re.compile(r'[<>:"/\\|?*]+')
VelocityEstimationMethod = Literal["doppler_moments", "frequency_bands"]
DEFAULT_VELOCITY_ESTIMATION_METHOD: VelocityEstimationMethod = "doppler_moments"
VELOCITY_ESTIMATION_METHODS: frozenset[VelocityEstimationMethod] = frozenset(
    {DEFAULT_VELOCITY_ESTIMATION_METHOD, "frequency_bands"}
)
# Public workflow names map to estimator/provenance identifiers. Extend this map
# when adding a workflow; "both" selects every registered workflow.
VELOCITY_WORKFLOW_METHODS = {
    "doppler_moments": "doppler_moments",
    "bandratio": "frequency_bands",
}
DEFAULT_VELOCITY_WORKFLOW_SELECTION = "both"


def validate_velocity_workflow_selection(value: object) -> str:
    # Retain compatibility with configurations using the former estimator name.
    if value == "frequency_bands":
        return "bandratio"
    if isinstance(value, str) and (value == "both" or value in VELOCITY_WORKFLOW_METHODS):
        return value
    allowed = ", ".join((*VELOCITY_WORKFLOW_METHODS, "both"))
    raise ValueError(f"Invalid velocity_estimation_method {value!r}. Expected one of: {allowed}.")


def selected_velocity_methods(value: object = DEFAULT_VELOCITY_WORKFLOW_SELECTION) -> tuple[str, ...]:
    selection = validate_velocity_workflow_selection(value)
    if selection == "both":
        return tuple(VELOCITY_WORKFLOW_METHODS.values())
    return (VELOCITY_WORKFLOW_METHODS[selection],)


PIPELINE_NAME_MIGRATIONS: dict[str, str] = {
    "retinal_velocity": "velocity",
    "waveform_velocity": "velocity_analysis",
}


def validate_band_limited_harmonics_setting(value: object) -> int:
    """Validate a JSON/run setting before output-directory replacement.

    Normalize integer settings to native integers for JSON serialization.
    Calculation argument validation belongs to calculations.math.bandlimiting.
    """
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
        raise ValueError("band_limited_harmonics must be a non-negative integer, excluding DC.")
    return int(value)


def validate_input_backend(value: object) -> str:
    if value not in ("auto", "direct", "replay"):
        raise ValueError("input_backend must be auto, direct or replay")
    return value


def validate_input_unknown_transport(value: object) -> str:
    if value not in ("coordinator", "worker"):
        raise ValueError("input_unknown_transport must be coordinator or worker")
    return value


def validate_copy_hd_passthrough(value: object) -> bool:
    if not isinstance(value, bool):
        raise ValueError("copy_hd_passthrough must be true or false.")
    return value


def validate_velocity_estimation_method(value: object) -> VelocityEstimationMethod:
    """Return a supported velocity estimator name or raise a clear error."""

    if isinstance(value, str) and value in VELOCITY_ESTIMATION_METHODS:
        return cast(VelocityEstimationMethod, value)
    allowed = ", ".join(sorted(VELOCITY_ESTIMATION_METHODS))
    raise ValueError(
        "Invalid velocity_estimation_method "
        f"{value!r}. Expected one of: {allowed}."
    )


def _read_version_from_pyproject(pyproject_path: Path) -> str | None:
    try:
        lines = pyproject_path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return None

    for line in lines:
        match = VERSION_PATTERN.match(line)
        if match:
            return match.group(1)
    return None


def app_version() -> str | None:
    env_version = os.getenv("EYEFLOW_VERSION", "").strip()
    if env_version:
        return env_version

    try:
        return importlib_metadata.version(APP_NAME)
    except importlib_metadata.PackageNotFoundError:
        pass

    for root in _resource_roots():
        version = _read_version_from_pyproject(root / "pyproject.toml")
        if version:
            return version
    return None


def app_display_name() -> str:
    version = app_version()
    return f"{APP_NAME} {version}" if version else APP_NAME


def _settings_subdir_name() -> str:
    version = app_version()
    if not version:
        return APP_NAME

    safe_version = INVALID_PATH_CHARS_PATTERN.sub("-", version).rstrip(" .")
    return safe_version or APP_NAME


def default_settings_path() -> Path:
    appdata = os.getenv("APPDATA") or os.getenv("LOCALAPPDATA")
    if appdata:
        base_dir = Path(appdata)
    else:
        xdg_config = os.getenv("XDG_CONFIG_HOME")
        base_dir = Path(xdg_config) if xdg_config else Path.home() / ".config"
    return base_dir / APP_NAME / _settings_subdir_name() / SETTINGS_FILENAME


def _resource_roots() -> list[Path]:
    roots: list[Path] = []
    frozen_root = getattr(sys, "_MEIPASS", None)
    if frozen_root:
        roots.append(Path(frozen_root))
    if getattr(sys, "frozen", False):
        roots.append(Path(sys.executable).resolve().parent)
    roots.append(Path(__file__).resolve().parents[1])
    roots.append(Path.cwd())
    return roots


def runtime_pipelines_path() -> Path:
    env_path = os.getenv(PIPELINES_DIR_ENV)
    if env_path:
        return Path(env_path).expanduser()

    if getattr(sys, "frozen", False):
        return Path(sys.executable).resolve().parent / "pipelines"

    return Path(__file__).resolve().parent / "pipelines"


def default_settings_template_path() -> Path | None:
    env_path = os.getenv("EYEFLOW_DEFAULT_SETTINGS")
    if env_path:
        candidate = Path(env_path).expanduser()
        if candidate.is_file():
            return candidate

    for root in _resource_roots():
        candidate = root / DEFAULT_SETTINGS_FILENAME
        if candidate.is_file():
            return candidate
    return None


def _load_settings_file(path: Path) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, OSError, json.JSONDecodeError):
        return {}
    return payload if isinstance(payload, dict) else {}


def normalize_named_visibility(
    item_names: Iterable[str], stored_visibility: Mapping[str, bool] | None
) -> tuple[dict[str, bool], bool]:
    ordered_names = list(dict.fromkeys(item_names))
    clean_stored = {
        name: value
        for name, value in (stored_visibility or {}).items()
        if isinstance(name, str) and isinstance(value, bool)
    }

    if not clean_stored:
        return {name: True for name in ordered_names}, bool(ordered_names)

    visibility: dict[str, bool] = {}
    changed = False
    for name in ordered_names:
        if name in clean_stored:
            visibility[name] = clean_stored[name]
        else:
            visibility[name] = False
            changed = True

    if set(clean_stored) != set(visibility):
        changed = True
    return visibility, changed


def normalize_pipeline_visibility(
    pipeline_names: Iterable[str],
    stored_visibility: Mapping[str, bool] | None,
    *,
    missing_defaults: Mapping[str, bool] | None = None,
) -> tuple[dict[str, bool], bool]:
    names = list(dict.fromkeys(pipeline_names))
    migrated_visibility, names_migrated = _migrate_pipeline_names(
        stored_visibility
    )
    visibility, changed = normalize_named_visibility(names, migrated_visibility)
    changed = changed or names_migrated
    stored_names = {
        name
        for name, value in migrated_visibility.items()
        if isinstance(name, str) and isinstance(value, bool)
    }
    if stored_names:
        for name in names:
            if name in stored_names:
                continue
            default = bool((missing_defaults or {}).get(name, False))
            if visibility.get(name) != default:
                visibility[name] = default
                changed = True
    return visibility, changed


def normalize_pipeline_options(
    options_by_pipeline: Mapping[str, Iterable[object]],
    stored_options: Mapping[str, Mapping[str, bool]] | None,
) -> tuple[dict[str, dict[str, bool]], bool]:
    """Normalize persisted option selections against the current catalog."""

    migrated_options, names_migrated = _migrate_pipeline_names(stored_options)
    clean_stored: dict[str, dict[str, bool]] = {}
    for pipeline_name, values in migrated_options.items():
        if not isinstance(pipeline_name, str) or not isinstance(values, Mapping):
            continue
        clean_stored[pipeline_name] = {
            name: enabled
            for name, enabled in values.items()
            if isinstance(name, str) and isinstance(enabled, bool)
        }

    normalized: dict[str, dict[str, bool]] = {}
    for pipeline_name, options in options_by_pipeline.items():
        pipeline_values: dict[str, bool] = {}
        stored_pipeline = clean_stored.get(pipeline_name, {})
        for option in options:
            option_name = str(getattr(option, "name", option)).strip()
            if not option_name:
                continue
            default_enabled = bool(getattr(option, "default_enabled", True))
            pipeline_values[option_name] = stored_pipeline.get(
                option_name,
                default_enabled,
            )
        if pipeline_values:
            normalized[pipeline_name] = pipeline_values

    changed = names_migrated or normalized != clean_stored
    return normalized, changed


def _migrate_pipeline_names(
    values: Mapping[str, Any] | None,
) -> tuple[dict[str, Any], bool]:
    """Move persisted pipeline settings to their current identifiers."""

    migrated = dict(values or {})
    changed = False
    for old_name, new_name in PIPELINE_NAME_MIGRATIONS.items():
        if old_name not in migrated:
            continue
        if new_name not in migrated:
            migrated[new_name] = migrated[old_name]
        del migrated[old_name]
        changed = True
    return migrated, changed


class AppSettingsStore:
    def __init__(
        self,
        path: Path | None = None,
        default_template_path: Path | None = None,
    ) -> None:
        self.path = path or default_settings_path()
        self.default_template_path = (
            default_template_path
            if default_template_path is not None
            else (default_settings_template_path() if path is None else None)
        )

    def load_defaults(self) -> dict[str, Any]:
        if self.default_template_path is None:
            return {}
        return _load_settings_file(self.default_template_path)

    def initialize_from_defaults(self) -> bool:
        if self.path.exists():
            return False
        defaults = self.load_defaults()
        if not defaults:
            return False
        self.save(defaults)
        return True

    def load(self) -> dict[str, Any]:
        if not self.path.exists():
            return self.load_defaults()
        return _load_settings_file(self.path)

    def save(self, settings: Mapping[str, Any]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        tmp_path = self.path.with_suffix(".tmp")
        settings = dict(settings)
        tmp_path.write_text(
            json.dumps(settings, indent=2, sort_keys=True),
            encoding="utf-8",
        )
        tmp_path.replace(self.path)

    def import_file(self, path: Path) -> None:
        try:
            settings = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise ValueError(
                f"Invalid JSON at line {exc.lineno}, column {exc.colno}: {exc.msg}"
            ) from exc
        if not isinstance(settings, dict):
            raise TypeError("The configuration must contain a JSON object.")
        count = settings.get(
            "number_of_radii_in_half_diagonal", DEFAULT_NUMBER_OF_RADII_IN_HALF_DIAGONAL,
        )
        if isinstance(count, bool) or not isinstance(count, Integral) or count < 1:
            raise ValueError("number_of_radii_in_half_diagonal must be a positive integer.")
        validate_band_limited_harmonics_setting(
            settings.get("band_limited_harmonics", DEFAULT_BAND_LIMITED_HARMONICS_SETTING)
        )
        validate_band_ratio_frequency_scale_hz(
            settings.get(
                "band_ratio_frequency_scale_hz",
                DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ,
            )
        )
        validate_copy_hd_passthrough(
            settings.get("copy_hd_passthrough", DEFAULT_COPY_HD_PASSTHROUGH)
        )
        validate_velocity_workflow_selection(
            settings.get("velocity_estimation_method", DEFAULT_VELOCITY_WORKFLOW_SELECTION)
        )
        validate_input_backend(settings.get("input_backend", "auto"))
        validate_input_unknown_transport(settings.get("input_unknown_transport", "coordinator"))
        self.save(settings)



_SettingValue = TypeVar("_SettingValue")


def _settings_snapshot(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {name: _settings_snapshot(item) for name, item in value.items()}
    return deepcopy(value)


class SettingsSection(MutableMapping[str, _SettingValue], Generic[_SettingValue]):
    """A live view into a settings group, retained across imports and reloads."""

    def __init__(self, settings: UserSettings, *path: str) -> None:
        self._settings = settings
        self._path = path

    def _mapping(self, *, create: bool = False) -> dict[str, Any]:
        values = self._settings._data()
        for name in self._path:
            child = values.get(name)
            if not isinstance(child, dict):
                if not create:
                    return {}
                child = values[name] = {}
            values = child
        return values

    def __getitem__(self, name: str) -> _SettingValue:
        return self._mapping()[name]

    def __setitem__(self, name: str, value: _SettingValue) -> None:
        self._mapping(create=True)[name] = self._prepare(name, value)

    def _prepare(self, name: str, value: _SettingValue) -> Any:
        return deepcopy(value)

    def __delitem__(self, name: str) -> None:
        del self._mapping()[name]

    def __iter__(self) -> Iterator[str]:
        return iter(self._mapping())

    def __len__(self) -> int:
        return len(self._mapping())

    def replace(self, values: Mapping[str, _SettingValue]) -> None:
        # Snapshot first: values may itself be this live view.
        replacement = {name: self._prepare(name, value) for name, value in values.items()}
        target = self._mapping(create=True)
        target.clear()
        target.update(replacement)


class BooleanSettings(SettingsSection[bool]):
    """Named Boolean selections for one settings group."""

    def _prepare(self, name: str, value: bool) -> bool:
        if not isinstance(name, str) or not isinstance(value, bool):
            raise ValueError("Setting selections require string names and Boolean values.")
        return value

    def selected_names(self) -> tuple[str, ...]:
        return tuple(name for name, enabled in self.items() if enabled)


class PipelineVisibilitySettings(BooleanSettings):
    """User-selected pipeline targets; execution dependencies remain run-local."""

    def normalize(self, pipelines: Iterable[PipelineDescriptor]) -> bool:
        pipelines = tuple(pipelines)
        values, changed = normalize_pipeline_visibility(
            (item.name for item in pipelines), self,
            missing_defaults={item.name: item.default_selected for item in pipelines},
        )
        self.replace(values)
        return changed

    def selected_targets(self, pipelines: Iterable[PipelineDescriptor]) -> tuple[str, ...]:
        return tuple(item.name for item in pipelines if self.get(item.name, False))


class PipelineOptionsSettings(SettingsSection[Mapping[str, bool]]):
    """Pipeline option groups, each exposing a live Boolean settings view."""

    def __getitem__(self, name: str) -> BooleanSettings:
        self._mapping()[name]  # Keep normal mapping semantics for missing pipelines.
        return BooleanSettings(self._settings, *self._path, name)

    def _prepare(self, name: str, value: Mapping[str, bool]) -> dict[str, bool]:
        if not isinstance(name, str) or not isinstance(value, Mapping):
            raise ValueError("Pipeline options require a pipeline name and a Boolean mapping.")
        values = dict(value)
        if any(not isinstance(key, str) or not isinstance(enabled, bool)
               for key, enabled in values.items()):
            raise ValueError("Pipeline options require string names and Boolean values.")
        return values

    def setdefault(self, name: str, default: Mapping[str, bool] | None = None) -> BooleanSettings:
        if name not in self:
            self[name] = {} if default is None else default
        return self[name]

    def normalize(self, pipelines: Iterable[PipelineDescriptor]) -> bool:
        values, changed = normalize_pipeline_options(
            {item.name: item.options for item in pipelines if item.options}, self,
        )
        self.replace(values)
        return changed

    def selected_options(self, pipelines: Iterable[PipelineDescriptor]) -> dict[str, tuple[str, ...]]:
        pipelines = tuple(pipelines)
        values, _ = normalize_pipeline_options(
            {item.name: item.options for item in pipelines if item.options}, self,
        )
        return {
            item.name: tuple(option.name for option in item.options
                             if values[item.name].get(option.name, False))
            for item in pipelines if item.options
        }


class UserSettings(AppSettingsStore, MutableMapping[str, Any]):
    """In-memory user settings with explicit persistence and structured groups.

    Application code imports ``user_settings`` below. Explicit instances with
    a custom path support isolated tests/tools. Reading is lazy; importing this
    module never loads or writes the user's configuration. Edits are shared
    immediately in this process; ``save()`` persists them and ``reload()`` reads
    external changes. Consumers read properties and persist edits with save().
    """

    def __init__(self, path: Path | None = None,
                 default_template_path: Path | None = None, *,
                 values: Mapping[str, Any] | None = None) -> None:
        super().__init__(path, default_template_path)
        self._values: dict[str, Any] | None = None
        self._pipeline_visibility = PipelineVisibilitySettings(self, "pipeline_visibility")
        self._pipeline_options = PipelineOptionsSettings(self, "pipeline_options")
        if values is not None:
            self._replace_values(values)

    def _data(self) -> dict[str, Any]:
        if self._values is None:
            self._replace_values(super().load())
        return cast(dict[str, Any], self._values)

    def _replace_values(self, values: Mapping[str, Any]) -> None:
        self._values = deepcopy(dict(values))
        raw_options = self._values.get("pipeline_options", {})
        options = {
            name: {key: value for key, value in values.items()
                   if isinstance(key, str) and isinstance(value, bool)}
            for name, values in raw_options.items()
            if isinstance(name, str) and isinstance(values, Mapping)
        } if isinstance(raw_options, Mapping) else {}
        raw_visibility = self._values.get("pipeline_visibility", {})
        visibility = {
            name: value for name, value in raw_visibility.items()
            if isinstance(name, str) and isinstance(value, bool)
        } if isinstance(raw_visibility, Mapping) else {}
        legacy = options.get("displacement_map", {})
        if legacy:
            options.setdefault("displacement_profiles", {"segment_maps": legacy.get("segment_maps", False)})
            if legacy.get("profiles") is True:
                visibility.setdefault("displacement_profiles", True)
        if "pipeline_options" in self._values:
            self._values["pipeline_options"] = options
        if visibility or "pipeline_visibility" in self._values:
            self._values["pipeline_visibility"] = visibility

    def load(self) -> dict[str, Any]:
        """Return an independent snapshot of the shared settings."""
        return deepcopy(self._data())

    def save(self, settings: Mapping[str, Any] | None = None) -> None:
        values = self.load() if settings is None else _settings_snapshot(settings)
        super().save(values)
        self._replace_values(values)

    def initialize_from_defaults(self) -> bool:
        if self.path.exists() or not self.load_defaults():
            return False
        self.save()
        return True

    def reload(self) -> None:
        self._replace_values(super().load())

    def for_run(self) -> UserSettings:
        """Return a validated in-memory copy for one run and its worker processes."""
        settings = UserSettings(path=self.path, values=self.load())
        for name in (
            "band_ratio_frequency_scale_hz", "band_limited_harmonics", "copy_hd_passthrough",
            "velocity_estimation_method", "input_backend", "input_unknown_transport",
            "number_of_radii_in_half_diagonal",
        ):
            settings[name] = getattr(settings, name)
        return settings

    def __getitem__(self, name: str) -> Any:
        value = self._data()[name]
        if name == "pipeline_visibility":
            return self.pipeline_visibility
        if name == "pipeline_options":
            return self.pipeline_options
        return value

    def __setitem__(self, name: str, value: Any) -> None:
        if name == "pipeline_visibility":
            self.pipeline_visibility.replace(value)
        elif name == "pipeline_options":
            self.pipeline_options.replace(value)
        else:
            self._data()[name] = deepcopy(value)

    def __delitem__(self, name: str) -> None:
        del self._data()[name]

    def __iter__(self) -> Iterator[str]:
        return iter(self._data())

    def __len__(self) -> int:
        return len(self._data())

    @property
    def pipeline_visibility(self) -> PipelineVisibilitySettings:
        self._data()
        return self._pipeline_visibility

    @property
    def pipeline_options(self) -> PipelineOptionsSettings:
        self._data()
        return self._pipeline_options

    @property
    def ui_mode(self) -> str:
        mode = self._data().get("ui_mode", "minimal")
        return mode if mode in {"minimal", "advanced"} else "minimal"

    @ui_mode.setter
    def ui_mode(self, value: str) -> None:
        self._data()["ui_mode"] = "advanced" if value == "advanced" else "minimal"

    @property
    def velocity_estimation_method(self) -> str:
        return validate_velocity_workflow_selection(
            self._data().get("velocity_estimation_method", DEFAULT_VELOCITY_WORKFLOW_SELECTION)
        )

    @velocity_estimation_method.setter
    def velocity_estimation_method(self, value: str) -> None:
        self._data()["velocity_estimation_method"] = validate_velocity_workflow_selection(value)

    @property
    def input_backend(self) -> str:
        return validate_input_backend(self._data().get("input_backend", "auto"))

    @input_backend.setter
    def input_backend(self, value: str) -> None:
        self._data()["input_backend"] = validate_input_backend(value)

    @property
    def input_unknown_transport(self) -> str:
        return validate_input_unknown_transport(
            self._data().get("input_unknown_transport", "coordinator")
        )

    @input_unknown_transport.setter
    def input_unknown_transport(self, value: str) -> None:
        self._data()["input_unknown_transport"] = validate_input_unknown_transport(value)

    @property
    def band_limited_harmonics(self) -> int:
        return validate_band_limited_harmonics_setting(
            self._data().get("band_limited_harmonics", DEFAULT_BAND_LIMITED_HARMONICS_SETTING)
        )

    @band_limited_harmonics.setter
    def band_limited_harmonics(self, value: int) -> None:
        self._data()["band_limited_harmonics"] = validate_band_limited_harmonics_setting(value)

    @property
    def band_ratio_frequency_scale_hz(self) -> float:
        return validate_band_ratio_frequency_scale_hz(
            self._data().get("band_ratio_frequency_scale_hz", DEFAULT_BAND_RATIO_FREQUENCY_SCALE_HZ)
        )

    @band_ratio_frequency_scale_hz.setter
    def band_ratio_frequency_scale_hz(self, value: float) -> None:
        self._data()["band_ratio_frequency_scale_hz"] = validate_band_ratio_frequency_scale_hz(value)

    @property
    def copy_hd_passthrough(self) -> bool:
        return validate_copy_hd_passthrough(
            self._data().get("copy_hd_passthrough", DEFAULT_COPY_HD_PASSTHROUGH)
        )

    @copy_hd_passthrough.setter
    def copy_hd_passthrough(self, value: bool) -> None:
        self._data()["copy_hd_passthrough"] = validate_copy_hd_passthrough(value)

    @property
    def last_input_directory(self) -> Path | None:
        value = self._data().get("last_input_directory")
        if not isinstance(value, str) or not value.strip():
            return None
        directory = Path(value)
        try:
            return directory if directory.is_absolute() and directory.is_dir() else None
        except (OSError, ValueError):
            return None

    @last_input_directory.setter
    def last_input_directory(self, value: Path) -> None:
        self._data()["last_input_directory"] = str(value.expanduser().resolve())

    @property
    def number_of_radii_in_half_diagonal(self) -> int:
        count = self._data().get(
            "number_of_radii_in_half_diagonal", DEFAULT_NUMBER_OF_RADII_IN_HALF_DIAGONAL,
        )
        if isinstance(count, bool) or not isinstance(count, Integral) or count < 1:
            raise ValueError("number_of_radii_in_half_diagonal must be a positive integer.")
        return int(count)

    @number_of_radii_in_half_diagonal.setter
    def number_of_radii_in_half_diagonal(self, value: int) -> None:
        self._data()["number_of_radii_in_half_diagonal"] = value


# One process-wide dataset, shared by GUI, CLI, and any other consumer.
user_settings = UserSettings()
