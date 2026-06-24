"""Run-log file ownership for EyeFlow."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import ClassVar

from app_settings import LAST_RUN_LOG_FILENAME


@dataclass
class RunLogger:
    _current: ClassVar["RunLogger | None"] = None

    settings_path: Path
    last_saved_path: Path | None = None

    @classmethod
    def configure(cls, settings_path: Path) -> "RunLogger":
        logger = cls(settings_path)
        cls._current = logger
        return logger

    @classmethod
    def current(cls) -> "RunLogger":
        if cls._current is None:
            raise RuntimeError("RunLogger has not been configured.")
        return cls._current

    @classmethod
    def reset_current(cls) -> None:
        cls._current = None

    @property
    def path(self) -> Path:
        return self.settings_path.with_name(LAST_RUN_LOG_FILENAME)

    def write_snapshot(self, text: str) -> Path:
        log_path = self.path
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(text, encoding="utf-8")
        self.last_saved_path = log_path
        return log_path

    def ensure_file(self) -> Path:
        log_path = self.path
        if not log_path.exists():
            self.write_snapshot("")
        return log_path


def configure_logger(settings_path: Path) -> RunLogger:
    return RunLogger.configure(settings_path)


def current_logger() -> RunLogger:
    return RunLogger.current()
