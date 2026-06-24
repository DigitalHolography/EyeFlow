from __future__ import annotations

import importlib
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


def main() -> Any:
    configure_numeric_threads()
    return _call_entry("eye_flow", "main")


def cli_main(argv: list[str] | None = None) -> Any:
    configure_numeric_threads()
    return _call_entry("cli", "main", argv)
