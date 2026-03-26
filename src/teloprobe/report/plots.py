"""Plotly figure generators for the HTML report."""

import json
from typing import Optional

import numpy as np
import pandas as pd

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


def plot_tl_distribution(read_level: pd.DataFrame) -> str:
    """Telomere length distribution histogram with KDE overlay."""
    if not HAS_PLOTLY:
        return ""

    good = read_level[
        (read_level["qc_status"] == "Good") &
        (read_level["telomere_length_bp"] > 0)
    ]

    if good.empty:
        return _empty_plot("No telomere length data")

    fig = go.Figure()

    fig.add_trace(go.Histogram(
        x=good["telomere_length_bp"],
        nbinsx=50,
        name="Telomere Length",
        marker_color="steelblue",
        opacity=0.7,
    ))

    median_tl = good["telomere_length_bp"].median()
    fig.add_vline(
        x=median_tl, line_dash="dash", line_color="red",
        annotation_text=f"Median: {median_tl:.0f} bp",
    )

    fig.update_layout(
        title="Telomere Length Distribution",
        xaxis_title="Telomere Length (bp)",
        yaxis_title="Count",
        template="plotly_white",
        height=400,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_arm_boxplots(arm_level: pd.DataFrame) -> str:
    """Per-chromosome-arm boxplot of telomere lengths."""
    if not HAS_PLOTLY or arm_level.empty:
        return _empty_plot("No arm-level data")

    fig = go.Figure()

    # Sort arms naturally
    from natsort import natsorted
    arms = natsorted(arm_level["chromosome_arm"].unique())

    for arm in arms:
        arm_data = arm_level[arm_level["chromosome_arm"] == arm]
        for _, row in arm_data.iterrows():
            fig.add_trace(go.Box(
                y=[row["median_tl"]],
                name=arm,
                marker_color="steelblue",
                boxpoints=False,
            ))

    fig.update_layout(
        title="Telomere Length by Chromosome Arm",
        yaxis_title="Telomere Length (bp)",
        xaxis_title="Chromosome Arm",
        template="plotly_white",
        height=500,
        showlegend=False,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_arm_violin(read_level: pd.DataFrame) -> str:
    """Per-chromosome-arm violin plots from read-level data."""
    if not HAS_PLOTLY:
        return ""

    good = read_level[
        (read_level["qc_status"] == "Good") &
        (read_level["chromosome_arm"] != "") &
        (read_level["telomere_length_bp"] > 0)
    ]

    if good.empty:
        return _empty_plot("No arm-assigned reads")

    fig = px.violin(
        good,
        x="chromosome_arm",
        y="telomere_length_bp",
        color="haplotype",
        box=True,
        title="Telomere Length by Chromosome Arm",
        template="plotly_white",
    )
    fig.update_layout(height=500)

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_tvr_heatmap(read_level: pd.DataFrame) -> str:
    """TVR composition heatmap across chromosome arms."""
    if not HAS_PLOTLY:
        return ""

    good = read_level[
        (read_level["qc_status"] == "Good") &
        (read_level["chromosome_arm"] != "")
    ]

    if good.empty or "tvr_fraction" not in good.columns:
        return _empty_plot("No TVR data")

    # Aggregate mean TVR fraction per arm
    tvr_by_arm = good.groupby("chromosome_arm")["tvr_fraction"].mean()

    fig = go.Figure(data=go.Bar(
        x=tvr_by_arm.index.tolist(),
        y=tvr_by_arm.values,
        marker_color="coral",
    ))

    fig.update_layout(
        title="TVR Fraction by Chromosome Arm",
        xaxis_title="Chromosome Arm",
        yaxis_title="Mean TVR Fraction",
        template="plotly_white",
        height=400,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_qc_waterfall(read_level: pd.DataFrame) -> str:
    """QC filter waterfall chart showing read attrition."""
    if not HAS_PLOTLY:
        return ""

    counts = read_level["qc_status"].value_counts()

    fig = go.Figure(data=go.Bar(
        x=counts.index.tolist(),
        y=counts.values,
        marker_color=["green" if s == "Good" else "salmon" for s in counts.index],
    ))

    fig.update_layout(
        title="QC Filter Summary",
        xaxis_title="QC Status",
        yaxis_title="Read Count",
        template="plotly_white",
        height=400,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_confidence_histogram(read_level: pd.DataFrame) -> str:
    """Distribution of confidence scores."""
    if not HAS_PLOTLY:
        return ""

    good = read_level[read_level["qc_status"] == "Good"]
    if good.empty:
        return _empty_plot("No data")

    fig = go.Figure(data=go.Histogram(
        x=good["confidence_score"],
        nbinsx=30,
        marker_color="mediumpurple",
    ))

    fig.update_layout(
        title="Confidence Score Distribution",
        xaxis_title="Confidence Score",
        yaxis_title="Count",
        template="plotly_white",
        height=350,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_downsampling_stability(stability_data: pd.DataFrame) -> str:
    """Downsampling stability plot."""
    if not HAS_PLOTLY or stability_data.empty:
        return ""

    fig = px.box(
        stability_data,
        x="fraction",
        y="median_tl",
        title="Downsampling Stability",
        template="plotly_white",
    )
    fig.update_layout(
        xaxis_title="Downsampling Fraction",
        yaxis_title="Median Telomere Length (bp)",
        height=350,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def _empty_plot(message: str) -> str:
    """Return placeholder HTML for empty plots."""
    return f'<div class="empty-plot"><p>{message}</p></div>'
