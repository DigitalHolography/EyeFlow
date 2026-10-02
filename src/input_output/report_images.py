"""Logical image references shared by PNG exporters and the PDF report."""

from pathlib import Path

REPORT_IMAGES_STATE = "report_image_paths"
ReportImagePaths = dict[tuple[str, str], Path]
