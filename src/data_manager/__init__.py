"""Stable scientific products, with physical persistence supplied by ProductStorage."""

from .handles import InvocationResource, ProductHandle
from .ownership import DataManager

__all__ = ["DataManager", "ProductHandle", "InvocationResource"]
