"""Report artifact writers and image encoders."""

from .html_report import write_html_report
from .pdf_report import write_pdf_report

__all__ = ["write_html_report", "write_pdf_report"]
