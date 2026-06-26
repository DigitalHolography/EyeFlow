"""Runtime containers for in-memory DopplerView reconstruction."""

from dataclasses import dataclass


@dataclass
class DopplerViewStepContext:
    """In-memory bridge for copied DopplerView steps, not an HDF5 access layer."""

    cache: dict[str, object]
    holodoppler_config: dict[str, object]
    dopplerview_config: dict[str, object]

    def require(self, key: str):
        if key not in self.cache:
            raise RuntimeError(f"Missing required context key: '{key}'")
        return self.cache[key]

    def set(self, key: str, value) -> None:
        self.cache[key] = value

    def hd_config_value(self, key: str, default=None):
        return self.holodoppler_config.get(key, default)

    def dv_config_section(self, section: str) -> dict[str, object]:
        value = self.dopplerview_config.get(section, {})
        return value if isinstance(value, dict) else {}

    def dv_config_value(self, section: str, key: str, default=None):
        return self.dv_config_section(section).get(key, default)
