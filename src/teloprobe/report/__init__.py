"""Module 8: HTML report generation with Plotly visualizations."""

from .builder import build_report
from .plots import (
    plot_tl_distribution,
    plot_arm_boxplots,
    plot_tvr_heatmap,
    plot_qc_waterfall,
    plot_confidence_histogram,
)

__all__ = [
    "build_report",
    "plot_tl_distribution",
    "plot_arm_boxplots",
    "plot_tvr_heatmap",
    "plot_qc_waterfall",
    "plot_confidence_histogram",
]
