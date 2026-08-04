"""Run-log file ownership for EyeFlow."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
import re
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


def emit_log(
    on_log: Callable[[str], None] | None,
    message: str,
    level: LogLevel | str | None = None,
) -> str:
    """Send a normalized log message to a callback and return the text sent."""

    formatted = format_log_message(message, level)
    if on_log is not None:
        on_log(formatted)
    return formatted


def log_debug(on_log: Callable[[str], None] | None, message: str) -> str:
    return emit_log(on_log, message, LogLevel.DEBUG)


def log_info(on_log: Callable[[str], None] | None, message: str) -> str:
    return emit_log(on_log, message, LogLevel.INFO)


def log_warning(on_log: Callable[[str], None] | None, message: str) -> str:
    return emit_log(on_log, message, LogLevel.WARNING)


def log_error(on_log: Callable[[str], None] | None, message: str) -> str:
    return emit_log(on_log, message, LogLevel.ERROR)


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
