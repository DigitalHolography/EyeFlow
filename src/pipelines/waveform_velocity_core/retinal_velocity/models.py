"""Runtime context for in-memory retinal velocity analysis."""

from dataclasses import dataclass


@dataclass
class VelocityAnalysisContext:
    """In-memory bridge between retinal velocity calculation stages."""

    cache: dict[str, object]
    holodoppler_config: dict[str, object]
    analysis_config: dict[str, object]

    def require(self, key: str):
        if key not in self.cache:
            raise RuntimeError(f"Missing required context key: '{key}'")
        return self.cache[key]

    def set(self, key: str, value) -> None:
        self.cache[key] = value

    def hd_config_value(self, key: str, default=None):
        return self.holodoppler_config.get(key, default)

    def analysis_config_section(self, section: str) -> dict[str, object]:
        value = self.analysis_config.get(section, {})
        return value if isinstance(value, dict) else {}

    def analysis_config_value(self, section: str, key: str, default=None):
        return self.analysis_config_section(section).get(key, default)
