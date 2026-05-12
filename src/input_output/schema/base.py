"""Pydantic models for declared HDF5 source schemas."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field, field_validator


class H5DatasetSpec(BaseModel):
    """One named HDF5 dataset in an upstream source file."""

    model_config = ConfigDict(frozen=True)

    key: str
    path: str
    dtype: str | None = None
    dims: tuple[str, ...] = ()
    required: bool = True
    description: str = ""

    @field_validator("path")
    @classmethod
    def normalize_path(cls, value: str) -> str:
        return value.replace("\\", "/").strip("/")

    @field_validator("dtype")
    @classmethod
    def normalize_dtype(cls, value: str | None) -> str | None:
        if value is None:
            return None
        normalized = value.lower()
        aliases = {"float": "float32", "int": "int32"}
        if normalized in aliases:
            return aliases[normalized]
        if "64" in normalized or normalized == "complex128":
            raise ValueError(f"Schema dtype '{value}' must use a 32-bit variant.")
        return normalized


class JsonConfigValueSpec(BaseModel):
    """One config value that may live in HDF5 or a sidecar JSON file."""

    model_config = ConfigDict(frozen=True)

    key: str
    json_key: str
    h5_path: str | None = None
    section: str | None = None
    default: Any = None
    description: str = ""

    @field_validator("h5_path")
    @classmethod
    def normalize_h5_path(cls, value: str | None) -> str | None:
        return None if value is None else value.replace("\\", "/").strip("/")

    def read_json_config(self, config: dict[str, object]) -> Any:
        source = config.get(self.section, {}) if self.section else config
        if not isinstance(source, dict):
            return self.default
        return source.get(self.json_key, self.default)


class H5SourceSchema(BaseModel):
    """Declared HDF5 source contract for a companion application."""

    model_config = ConfigDict(frozen=True)

    label: str
    companion_suffix: str
    h5_folder_name: str
    h5_filename_template: str
    config_dir_name: str | None = None
    config_filename: str | None = None
    datasets: dict[str, H5DatasetSpec] = Field(default_factory=dict)
    config_values: dict[str, JsonConfigValueSpec] = Field(default_factory=dict)

    def dataset(self, key: str) -> H5DatasetSpec:
        try:
            return self.datasets[key]
        except KeyError as exc:
            raise KeyError(f"{self.label} schema has no dataset '{key}'.") from exc

    def dataset_path(self, key: str) -> str:
        return self.dataset(key).path

    def config_value(self, key: str) -> JsonConfigValueSpec:
        try:
            return self.config_values[key]
        except KeyError as exc:
            message = f"{self.label} schema has no config value '{key}'."
            raise KeyError(message) from exc
