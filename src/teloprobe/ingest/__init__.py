"""Module 1: Input ingestion and preprocessing."""

from .reader import read_input, detect_format
from .platform_detect import detect_platform, get_platform_preset
from .sample_sheet import parse_sample_sheet

__all__ = [
    "read_input",
    "detect_format",
    "detect_platform",
    "get_platform_preset",
    "parse_sample_sheet",
]
