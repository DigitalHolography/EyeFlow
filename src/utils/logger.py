"""Process-wide run logging for EyeFlow."""

from __future__ import annotations

import os
import re
import sys
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Callable, ClassVar

from app_settings import LAST_RUN_LOG_FILENAME
from input_output.writers.sidecars import write_text_file


class LogLevel(str, Enum):
    """Levels used by the run-log event stream."""

    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"
    PERF = "PERF"


_LEVEL_MARKER = re.compile(r"^\[(DEBUG|INFO|WARNING|ERROR|PERF|FAIL)\]\s*", re.IGNORECASE)
_LOCAL_LOG_SINK = ContextVar("eyeflow_log_sink", default=None)
_DIAGNOSTICS_ENV = "EYEFLOW_DIAGNOSTICS"
_PERF_ENV = "EYEFLOW_PERF_LOGGING"
_LEVEL_ENV = "EYEFLOW_LOG_LEVEL"
_LEVEL_ORDER = {
    LogLevel.DEBUG: 10,
    LogLevel.INFO: 20,
    LogLevel.WARNING: 30,
    LogLevel.ERROR: 40,
}


def _diagnostics_enabled() -> bool:
    configured = os.environ.get(_DIAGNOSTICS_ENV)
    if configured is not None:
        return configured.strip().lower() in {"1", "true", "yes", "on"}
    return not bool(getattr(sys, "frozen", False))


def _performance_enabled() -> bool:
    configured = os.environ.get(_PERF_ENV)
    if configured is not None:
        return configured.strip().lower() in {"1", "true", "yes", "on"}
    return _diagnostics_enabled()


def _visible(level: LogLevel | None) -> bool:
    if level is None:
        return True
    if level == LogLevel.PERF:
        return _performance_enabled()
    configured = os.environ.get(_LEVEL_ENV)
    if configured is None:
        configured = "DEBUG" if _diagnostics_enabled() else "INFO"
    minimum = _coerce_level(configured)
    if minimum == LogLevel.PERF:
        raise ValueError("PERF is independently filtered and cannot be a minimum log level.")
    return _LEVEL_ORDER[level] >= _LEVEL_ORDER[minimum]


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

    Existing ``[FAIL]``/``[WARNING]``/``[DEBUG]``/``[PERF]`` messages are
    accepted for compatibility and normalized. New call sites should use one of the
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
    diagnostics: bool | None = None

    @classmethod
    def configure(
        cls,
        settings_path: Path | None = None,
        *,
        on_log: Callable[[str], None] | None = None,
        diagnostics: bool | None = None,
        perf: bool | None = None,
        level: LogLevel | str | None = None,
    ) -> "Logger":
        """Configure the process-wide policy and its inheritance by spawned workers.

        Development defaults to DEBUG plus PERF. Frozen builds default to INFO;
        ``EYEFLOW_DIAGNOSTICS=1`` or ``diagnostics=True`` enables diagnostics.
        ``EYEFLOW_PERF_LOGGING`` or ``perf=`` filters PERF independently.
        ``EYEFLOW_LOG_LEVEL`` or ``level=`` selects a minimum severity separately.
        """
        if diagnostics is not None:
            os.environ[_DIAGNOSTICS_ENV] = "1" if diagnostics else "0"
            if level is None:
                os.environ[_LEVEL_ENV] = "DEBUG" if diagnostics else "INFO"
            if perf is None:
                os.environ[_PERF_ENV] = "1" if diagnostics else "0"
        elif _DIAGNOSTICS_ENV not in os.environ:
            os.environ[_DIAGNOSTICS_ENV] = (
                "0" if bool(getattr(sys, "frozen", False)) else "1"
            )
        if perf is not None:
            os.environ[_PERF_ENV] = "1" if perf else "0"
        elif _PERF_ENV not in os.environ:
            os.environ[_PERF_ENV] = "1" if _diagnostics_enabled() else "0"
        if level is not None:
            selected_level = _coerce_level(level)
            if selected_level == LogLevel.PERF:
                raise ValueError("Use diagnostics=True to enable PERF logging.")
            os.environ[_LEVEL_ENV] = selected_level.value
        elif _LEVEL_ENV not in os.environ:
            os.environ[_LEVEL_ENV] = "DEBUG" if _diagnostics_enabled() else "INFO"
        if cls._current is None:
            cls._current = cls(
                settings_path=settings_path, on_log=on_log, diagnostics=diagnostics,
            )
        else:
            if settings_path is not None:
                cls._current.settings_path = settings_path
            cls._current.on_log = on_log
            if diagnostics is not None:
                cls._current.diagnostics = diagnostics
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
        formatted = format_log_message(message, level)
        marker = _LEVEL_MARKER.match(str(message))
        selected_level = _coerce_level(level) if level is not None else (
            LogLevel.ERROR if marker and marker.group(1).upper() == "FAIL"
            else LogLevel(marker.group(1).upper()) if marker else LogLevel.INFO
        )
        if not _visible(selected_level):
            return formatted
        sink = _LOCAL_LOG_SINK.get()
        if sink is not None:
            sink(formatted)
            return formatted
        if cls._current is None:
            return formatted
        return cls._current._emit(message, level)

    @classmethod
    @contextmanager
    def capture_messages(cls, sink):
        """Route logs locally without replacing another thread's run logger."""
        token = _LOCAL_LOG_SINK.set(sink)
        try:
            yield
        finally:
            _LOCAL_LOG_SINK.reset(token)

    @classmethod
    def log_debug(cls, message: str) -> str:
        return cls.log(message, LogLevel.DEBUG)

    @classmethod
    def log_perf(cls, message: str) -> str:
        return cls.log(message, LogLevel.PERF)

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
        log_path = write_text_file(self.path, text)
        self.last_saved_path = log_path
        return log_path

    def ensure_file(self) -> Path:
        log_path = self.path
        if not log_path.exists():
            self.write_snapshot("")
        return log_path

# Compatibility alias for code that still refers to the old class name.
RunLogger = Logger
