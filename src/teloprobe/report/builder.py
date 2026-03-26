"""HTML report builder using Jinja2 templates."""

import json
import logging
import warnings
from datetime import datetime
from pathlib import Path
from typing import Optional

import pandas as pd

from ..config import Config
from . import plots

logger = logging.getLogger(__name__)

try:
    from jinja2 import Environment, FileSystemLoader, BaseLoader
    HAS_JINJA2 = True
except ImportError:
    HAS_JINJA2 = False

try:
    import plotly.offline as _plotly_offline
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


def _plotly_script_tag(offline: bool = False) -> str:
    """Return the <script> tag for Plotly JS.

    Args:
        offline: If True, embed the full Plotly JS library inline.
                 If False, use the CDN link.

    Returns:
        An HTML <script> tag string.
    """
    if offline:
        if not HAS_PLOTLY:
            warnings.warn(
                "Offline mode requested but plotly is not installed. "
                "Falling back to CDN. Install plotly with: "
                "pip install plotly",
                stacklevel=3,
            )
            return '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>'
        js_code = _plotly_offline.get_plotlyjs()
        return f"<script>{js_code}</script>"
    return '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>'


# Inline template (fallback if template files not found)
REPORT_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>TeloPROBE Report - {{ sample_id }}</title>
    {{ plotly_script }}
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: #f8f9fa;
            color: #333;
        }
        h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
        h2 { color: #34495e; margin-top: 30px; }
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin: 20px 0;
        }
        .summary-card {
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            text-align: center;
        }
        .summary-card .value {
            font-size: 2em;
            font-weight: bold;
            color: #3498db;
        }
        .summary-card .label {
            color: #7f8c8d;
            font-size: 0.9em;
            margin-top: 5px;
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        th, td {
            padding: 10px 15px;
            text-align: left;
            border-bottom: 1px solid #ecf0f1;
        }
        th {
            background: #2c3e50;
            color: white;
            font-weight: 600;
        }
        tr:hover { background: #f5f6fa; }
        .plot-container {
            background: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            margin: 15px 0;
        }
        .qc-pass { color: #27ae60; font-weight: bold; }
        .qc-fail { color: #e74c3c; font-weight: bold; }
        .qc-warn { color: #f39c12; font-weight: bold; }
        .section { margin-bottom: 40px; }
        .footer {
            text-align: center;
            color: #95a5a6;
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ecf0f1;
        }
        .empty-plot {
            text-align: center;
            color: #95a5a6;
            padding: 40px;
            background: #f5f6fa;
            border-radius: 8px;
        }
    </style>
</head>
<body>
    <h1>TeloPROBE Report</h1>
    <p>Sample: <strong>{{ sample_id }}</strong> | Platform: {{ platform }} | Mode: {{ mode }} | Generated: {{ timestamp }}</p>

    <div class="section">
        <h2>Summary</h2>
        <div class="summary-grid">
            <div class="summary-card">
                <div class="value">{{ total_reads }}</div>
                <div class="label">Total Reads</div>
            </div>
            <div class="summary-card">
                <div class="value">{{ informative_reads }}</div>
                <div class="label">Informative Reads</div>
            </div>
            <div class="summary-card">
                <div class="value">{{ median_tl }} bp</div>
                <div class="label">Median Telomere Length</div>
            </div>
            <div class="summary-card">
                <div class="value">{{ ci_lower }}-{{ ci_upper }} bp</div>
                <div class="label">95% Confidence Interval</div>
            </div>
            {% if arms_assigned > 0 %}
            <div class="summary-card">
                <div class="value">{{ arms_assigned }}</div>
                <div class="label">Chromosome Arms</div>
            </div>
            {% endif %}
        </div>
    </div>

    <div class="section">
        <h2>Telomere Length Distribution</h2>
        <div class="plot-container">
            {{ tl_distribution_plot }}
        </div>
    </div>

    <div class="section">
        <h2>QC Summary</h2>
        <div class="plot-container">
            {{ qc_waterfall_plot }}
        </div>
        {{ qc_table }}
    </div>

    {% if arm_boxplot %}
    <div class="section">
        <h2>Per-Chromosome Arm Analysis</h2>
        <div class="plot-container">
            {{ arm_boxplot }}
        </div>
        {% if arm_table %}
        {{ arm_table }}
        {% endif %}
    </div>
    {% endif %}

    {% if tvr_heatmap %}
    <div class="section">
        <h2>Telomere Variant Repeats</h2>
        <div class="plot-container">
            {{ tvr_heatmap }}
        </div>
    </div>
    {% endif %}

    <div class="section">
        <h2>Confidence Scores</h2>
        <div class="plot-container">
            {{ confidence_plot }}
        </div>
    </div>

    {% if qc_checks %}
    <div class="section">
        <h2>QC Checks</h2>
        <table>
            <tr><th>Check</th><th>Status</th><th>Value</th><th>Details</th></tr>
            {% for check in qc_checks %}
            <tr>
                <td>{{ check.name }}</td>
                <td class="{{ 'qc-pass' if check.passed else 'qc-fail' }}">
                    {{ 'PASS' if check.passed else 'FAIL' }}
                </td>
                <td>{{ check.value }}</td>
                <td>{{ check.message }}</td>
            </tr>
            {% endfor %}
        </table>
    </div>
    {% endif %}

    {% if downsampling_plot %}
    <div class="section">
        <h2>Downsampling Stability</h2>
        <div class="plot-container">
            {{ downsampling_plot }}
        </div>
    </div>
    {% endif %}

    <div class="footer">
        <p>Generated by TeloPROBE v{{ version }} | {{ timestamp }}</p>
    </div>
</body>
</html>
"""


def build_report(
    read_level: pd.DataFrame,
    arm_level: pd.DataFrame,
    sample_level: pd.DataFrame,
    output_path: str | Path,
    config: Config,
    qc_checks: list = None,
    stability_data: pd.DataFrame = None,
    offline: bool = False,
) -> None:
    """Build the HTML report from analysis results.

    Args:
        read_level: Read-level DataFrame.
        arm_level: Arm-level DataFrame.
        sample_level: Sample-level DataFrame.
        output_path: Output HTML file path.
        config: Pipeline configuration.
        qc_checks: List of QCCheck objects.
        stability_data: Downsampling stability data.
        offline: If True, embed the full Plotly JS library inline
                 instead of loading from CDN. Requires plotly to be
                 installed. Falls back to CDN with a warning if plotly
                 is not available.
    """
    # Allow config to override the offline parameter
    if getattr(config, "offline", False):
        offline = True
    from teloprobe import __version__

    # Extract sample summary
    if not sample_level.empty:
        row = sample_level.iloc[0]
        sample_id = row.get("sample", config.sample_id)
        total_reads = row.get("total_reads", 0)
        informative = row.get("informative_reads", 0)
        median_tl = row.get("median_tl", 0)
        ci_lower = row.get("global_ci_lower", 0)
        ci_upper = row.get("global_ci_upper", 0)
        arms = row.get("arms_assigned", 0)
    else:
        sample_id = config.sample_id
        total_reads = len(read_level)
        informative = (read_level["qc_status"] == "Good").sum()
        good_tl = read_level[read_level["qc_status"] == "Good"]["telomere_length_bp"]
        median_tl = good_tl.median() if len(good_tl) > 0 else 0
        ci_lower = ci_upper = median_tl
        arms = 0

    # Generate plots
    tl_dist = plots.plot_tl_distribution(read_level)
    qc_waterfall = plots.plot_qc_waterfall(read_level)
    confidence = plots.plot_confidence_histogram(read_level)

    arm_box = ""
    tvr_heat = ""
    if not arm_level.empty:
        arm_box = plots.plot_arm_violin(read_level)
        tvr_heat = plots.plot_tvr_heatmap(read_level)

    ds_plot = ""
    if stability_data is not None and not stability_data.empty:
        ds_plot = plots.plot_downsampling_stability(stability_data)

    # QC summary table
    from ..qc.read_qc import qc_summary_table
    qc_df = qc_summary_table(read_level)
    qc_table_html = qc_df.to_html(index=False, classes="qc-table")

    # Arm-level table
    arm_table_html = ""
    if not arm_level.empty:
        arm_table_html = arm_level.to_html(index=False, classes="arm-table")

    # Plotly script tag (CDN or inline)
    plotly_script = _plotly_script_tag(offline=offline)

    # Render template
    if HAS_JINJA2:
        env = Environment(loader=BaseLoader())
        template = env.from_string(REPORT_TEMPLATE)
        html = template.render(
            sample_id=sample_id,
            platform=config.platform.value,
            mode=config.mode.value,
            plotly_script=plotly_script,
            timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            version=__version__,
            total_reads=total_reads,
            informative_reads=informative,
            median_tl=f"{median_tl:.0f}",
            ci_lower=f"{ci_lower:.0f}",
            ci_upper=f"{ci_upper:.0f}",
            arms_assigned=arms,
            tl_distribution_plot=tl_dist,
            qc_waterfall_plot=qc_waterfall,
            qc_table=qc_table_html,
            arm_boxplot=arm_box,
            arm_table=arm_table_html,
            tvr_heatmap=tvr_heat,
            confidence_plot=confidence,
            qc_checks=qc_checks or [],
            downsampling_plot=ds_plot,
        )
    else:
        html = f"<html><body><h1>TeloPROBE Report - {sample_id}</h1><p>Install jinja2 for full report.</p></body></html>"

    Path(output_path).write_text(html)
    logger.info("Report written to %s", output_path)
