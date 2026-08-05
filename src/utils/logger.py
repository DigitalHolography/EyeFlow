"""Process-wide run logging for EyeFlow."""

from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Callable, ClassVar

from app_settings import LAST_RUN_LOG_FILENAME


class LogLevel(str, Enum):
    """Levels used by the run-log event stream."""

    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"


_LEVEL_MARKER = re.compile(r"^\[(DEBUG|INFO|WARNING|ERROR|FAIL)\]\s*", re.IGNORECASE)


def _coerce_level(level: LogLevel | str) -> LogLevel:
    if isinstance(level, LogLevel):
        return level
    try:
        return LogLevel(str(level).upper())
    except ValueError as exc:
        valid_levels = ", ".join(item.value for item in LogLevel)
        raise ValueError(
            f"Unknown log level {level!r}; expected one of {valid_levels}."
        ) from exc


def format_log_message(
    message: str,
    level: LogLevel | str | None = None,
) -> str:
    """Normalize a log message without requiring embedded level markers.

    Existing ``[FAIL]``/``[WARNING]``/``[DEBUG]`` messages are accepted for
    compatibility and normalized. New call sites should use one of the
    level-specific helpers below instead of embedding a marker in the text.
    Messages without a level remain unchanged so labels such as ``[INPUT]``
    and ``[DAG]`` keep their current format.
    """

    text = str(message)
    marker = _LEVEL_MARKER.match(text)
    inferred_level: LogLevel | None = None
    if marker is not None:
        inferred_level = (
            LogLevel.ERROR
            if marker.group(1).upper() == "FAIL"
            else LogLevel(marker.group(1).upper())
        )
        text = text[marker.end() :]

    selected_level = _coerce_level(level) if level is not None else inferred_level
    if selected_level is None:
        return text
    return f"[{selected_level.value}] {text}" if text else f"[{selected_level.value}]"


@dataclass
class Logger:
    """Process-wide logger and owner of the run-log snapshot file."""

    _current: ClassVar["Logger | None"] = None

    settings_path: Path | None = None
    on_log: Callable[[str], None] | None = None
    last_saved_path: Path | None = None

    @classmethod
    def configure(
        cls,
        settings_path: Path | None = None,
        *,
        on_log: Callable[[str], None] | None = None,
    ) -> "Logger":
        if cls._current is None:
            cls._current = cls(settings_path=settings_path, on_log=on_log)
        else:
            if settings_path is not None:
                cls._current.settings_path = settings_path
            cls._current.on_log = on_log
        return cls._current

    @classmethod
    def current(cls) -> "Logger":
        if cls._current is None:
            raise RuntimeError("Logger has not been configured.")
        return cls._current

    @classmethod
    def reset_current(cls) -> None:
        cls._current = None

    @property
    def path(self) -> Path:
        if self.settings_path is None:
            raise RuntimeError("Logger has no settings path.")
        return self.settings_path.with_name(LAST_RUN_LOG_FILENAME)

    def _emit(
        self,
        message: str,
        level: LogLevel | str | None = None,
    ) -> str:
        formatted = format_log_message(message, level)
        if self.on_log is not None:
            self.on_log(formatted)
        return formatted

    @classmethod
    def log(
        cls,
        message: str,
        level: LogLevel | str | None = None,
    ) -> str:
        if cls._current is None:
            return format_log_message(message, level)
        return cls._current._emit(message, level)

    @classmethod
    def log_debug(cls, message: str) -> str:
        return cls.log(message, LogLevel.DEBUG)

    @classmethod
    def log_info(cls, message: str) -> str:
        return cls.log(message, LogLevel.INFO)

    @classmethod
    def log_warning(cls, message: str) -> str:
        return cls.log(message, LogLevel.WARNING)

    @classmethod
    def log_error(cls, message: str) -> str:
        return cls.log(message, LogLevel.ERROR)

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

# Compatibility alias for code that still refers to the old class name.
RunLogger = Logger
