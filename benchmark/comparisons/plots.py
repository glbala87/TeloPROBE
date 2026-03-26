"""Benchmark comparison plots using Plotly."""

import logging
from typing import Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


def plot_tl_correlation(
    result_a: pd.DataFrame,
    result_b: pd.DataFrame,
    tool_a: str,
    tool_b: str,
) -> Optional[str]:
    """Scatter plot of per-read telomere lengths between two tools."""
    if not HAS_PLOTLY or result_a is None or result_b is None:
        return None

    merged = pd.merge(
        result_a[["read_id", "telomere_length"]],
        result_b[["read_id", "telomere_length"]],
        on="read_id", suffixes=(f"_{tool_a}", f"_{tool_b}"),
        how="inner",
    )
    if merged.empty:
        return None

    col_a = f"telomere_length_{tool_a}"
    col_b = f"telomere_length_{tool_b}"

    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=merged[col_a], y=merged[col_b],
        mode="markers", marker=dict(size=3, opacity=0.4, color="steelblue"),
        name="Reads",
    ))

    # Identity line
    max_val = max(merged[col_a].max(), merged[col_b].max())
    fig.add_trace(go.Scatter(
        x=[0, max_val], y=[0, max_val],
        mode="lines", line=dict(dash="dash", color="red"),
        name="y = x",
    ))

    fig.update_layout(
        title=f"Per-Read TL: {tool_a} vs {tool_b} (n={len(merged)})",
        xaxis_title=f"{tool_a} Telomere Length (bp)",
        yaxis_title=f"{tool_b} Telomere Length (bp)",
        template="plotly_white", height=500, width=550,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_bland_altman(
    result_a: pd.DataFrame,
    result_b: pd.DataFrame,
    tool_a: str,
    tool_b: str,
) -> Optional[str]:
    """Bland-Altman plot for agreement analysis."""
    if not HAS_PLOTLY or result_a is None or result_b is None:
        return None

    merged = pd.merge(
        result_a[["read_id", "telomere_length"]],
        result_b[["read_id", "telomere_length"]],
        on="read_id", suffixes=("_a", "_b"),
        how="inner",
    )
    if merged.empty:
        return None

    mean_tl = (merged["telomere_length_a"] + merged["telomere_length_b"]) / 2
    diff = merged["telomere_length_a"] - merged["telomere_length_b"]
    mean_diff = diff.mean()
    std_diff = diff.std()

    fig = go.Figure()
    fig.add_trace(go.Scattergl(
        x=mean_tl, y=diff,
        mode="markers", marker=dict(size=3, opacity=0.4, color="mediumpurple"),
        name="Reads",
    ))

    # Mean difference line
    fig.add_hline(y=mean_diff, line_dash="solid", line_color="red",
                  annotation_text=f"Mean: {mean_diff:.0f}")
    # LoA
    fig.add_hline(y=mean_diff + 1.96 * std_diff, line_dash="dash", line_color="gray",
                  annotation_text=f"+1.96 SD: {mean_diff + 1.96 * std_diff:.0f}")
    fig.add_hline(y=mean_diff - 1.96 * std_diff, line_dash="dash", line_color="gray",
                  annotation_text=f"-1.96 SD: {mean_diff - 1.96 * std_diff:.0f}")

    fig.update_layout(
        title=f"Bland-Altman: {tool_a} vs {tool_b}",
        xaxis_title="Mean TL (bp)",
        yaxis_title=f"Difference ({tool_a} - {tool_b}) (bp)",
        template="plotly_white", height=450, width=550,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_arm_comparison(
    results: dict[str, pd.DataFrame],
) -> Optional[str]:
    """Grouped bar chart of median TL per arm across tools."""
    if not HAS_PLOTLY:
        return None

    fig = go.Figure()

    for tool_name, arm_df in results.items():
        if arm_df is None or arm_df.empty:
            continue
        sorted_df = arm_df.sort_values("chromosome_arm")
        fig.add_trace(go.Bar(
            x=sorted_df["chromosome_arm"],
            y=sorted_df["median_tl"],
            name=tool_name,
            opacity=0.8,
        ))

    fig.update_layout(
        title="Median Telomere Length by Chromosome Arm",
        xaxis_title="Chromosome Arm",
        yaxis_title="Median TL (bp)",
        barmode="group",
        template="plotly_white",
        height=500,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_tl_distributions(
    results: dict[str, pd.DataFrame],
) -> Optional[str]:
    """Overlaid TL distribution histograms across tools."""
    if not HAS_PLOTLY:
        return None

    fig = go.Figure()
    colors = ["steelblue", "coral", "mediumseagreen", "mediumpurple"]

    for i, (tool_name, read_df) in enumerate(results.items()):
        if read_df is None or read_df.empty:
            continue
        tl = read_df["telomere_length"]
        tl = tl[tl > 0]
        fig.add_trace(go.Histogram(
            x=tl, name=tool_name, opacity=0.5,
            marker_color=colors[i % len(colors)],
            nbinsx=50,
        ))

    fig.update_layout(
        title="Telomere Length Distributions Across Tools",
        xaxis_title="Telomere Length (bp)",
        yaxis_title="Count",
        barmode="overlay",
        template="plotly_white",
        height=450,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_runtime_comparison(
    runtimes: dict[str, float],
) -> Optional[str]:
    """Bar chart of runtime per tool."""
    if not HAS_PLOTLY or not runtimes:
        return None

    tools = list(runtimes.keys())
    times = [runtimes[t] for t in tools]

    fig = go.Figure(data=go.Bar(
        x=tools, y=times,
        marker_color=["steelblue", "coral", "mediumseagreen", "mediumpurple"][:len(tools)],
    ))

    fig.update_layout(
        title="Runtime Comparison",
        yaxis_title="Seconds",
        template="plotly_white",
        height=350,
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)


def plot_metrics_heatmap(
    metrics_table: pd.DataFrame,
) -> Optional[str]:
    """Heatmap of comparison metrics across tool pairs."""
    if not HAS_PLOTLY or metrics_table.empty:
        return None

    numeric_cols = metrics_table.select_dtypes(include=[np.number]).columns
    if len(numeric_cols) == 0:
        return None

    fig = go.Figure(data=go.Heatmap(
        z=metrics_table[numeric_cols].values,
        x=numeric_cols.tolist(),
        y=metrics_table.get("comparison", metrics_table.index).tolist(),
        colorscale="RdYlGn",
        text=np.round(metrics_table[numeric_cols].values, 3),
        texttemplate="%{text}",
        hovertemplate="Metric: %{x}<br>Pair: %{y}<br>Value: %{z:.4f}<extra></extra>",
    ))

    fig.update_layout(
        title="Concordance Metrics Heatmap",
        template="plotly_white",
        height=300 + 30 * len(metrics_table),
    )

    return fig.to_html(full_html=False, include_plotlyjs=False)
