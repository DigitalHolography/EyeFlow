"""Presentation-only formatting for physical units used in figures."""

from __future__ import annotations


def display_unit(unit: str) -> str:
    """Return a compact MathText unit without changing stored unit metadata."""

    normalized = str(unit).strip().replace("Â³", "^3").replace("³", "^3")
    return {
        "mm/s": r"$\mathrm{mm.s^{-1}}$",
        "mm^3/s": r"$\mathrm{mm^3.s^{-1}}$",
    }.get(normalized, str(unit).strip())


__all__ = ["display_unit"]
