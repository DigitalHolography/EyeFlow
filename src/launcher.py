from __future__ import annotations

import importlib
import sys
from typing import Any

from runtime_limits import configure_numeric_threads


def _call_entry(
    module_name: str,
    func_name: str,
    *args: Any,
    **kwargs: Any,
) -> Any:
    module = importlib.import_module(module_name)
    entry = getattr(module, func_name)
    return entry(*args, **kwargs)


def main(argv: list[str] | None = None) -> Any:
    """Launch the GUI with no arguments, otherwise dispatch to the CLI."""

    configure_numeric_threads()
    effective_argv = list(sys.argv[1:] if argv is None else argv)
    if effective_argv:
        return _call_entry("cli", "main", effective_argv)
    return _call_entry("eye_flow", "main")


def cli_main(argv: list[str] | None = None) -> Any:
    """Dedicated CLI entry point; never fall back to the GUI."""

    configure_numeric_threads()
    return _call_entry("cli", "main", argv)
