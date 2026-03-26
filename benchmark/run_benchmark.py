#!/usr/bin/env python3
"""TeloPROBE Benchmark Runner.

Runs TeloPROBE against Telogator2, Topsicle, and wf-teloseq on the same
input data and generates a comprehensive comparison report.

Usage:
    # Run all tools on real data
    python run_benchmark.py --input reads.bam --output benchmark_results/ \
        --platform ont --reference subtelo_ref.fa --threads 8

    # Run on simulated data with ground truth evaluation
    python run_benchmark.py --simulate --output benchmark_results/ \
        --platform ont --n-reads 1000 --threads 8

    # Run only specific tools
    python run_benchmark.py --input reads.bam --output benchmark_results/ \
        --tools teloprobe,telogator2

    # Parse existing results (skip running tools)
    python run_benchmark.py --parse-only --output benchmark_results/
"""

import argparse
import json
import logging
import sys
import time
from datetime import datetime
from itertools import combinations
from pathlib import Path

import pandas as pd

from wrappers.base import ToolResult
from wrappers.teloprobe_wrapper import TeloPROBEWrapper
from wrappers.telogator2_wrapper import Telogator2Wrapper
from wrappers.topsicle_wrapper import TopsicleWrapper
from wrappers.wf_teloseq_wrapper import WfTeloseqWrapper
from comparisons.metrics import (
    compare_read_level,
    compare_arm_level,
    compare_sample_level,
    rank_tools,
)
from comparisons.plots import (
    plot_tl_correlation,
    plot_bland_altman,
    plot_arm_comparison,
    plot_tl_distributions,
    plot_runtime_comparison,
    plot_metrics_heatmap,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
)
logger = logging.getLogger("benchmark")

TOOL_REGISTRY = {
    "teloprobe": TeloPROBEWrapper,
    "telogator2": Telogator2Wrapper,
    "topsicle": TopsicleWrapper,
    "wf-teloseq": WfTeloseqWrapper,
}


def main():
    parser = argparse.ArgumentParser(
        description="TeloPROBE Benchmark: compare against Telogator2, Topsicle, and wf-teloseq",
    )
    parser.add_argument("-i", "--input", help="Input BAM/CRAM/FASTQ file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--platform", default="ont", choices=["ont", "pacbio"])
    parser.add_argument("-r", "--reference", help="Subtelomeric reference FASTA")
    parser.add_argument("-t", "--threads", type=int, default=4)
    parser.add_argument(
        "--tools", default="teloprobe,telogator2,topsicle,wf-teloseq",
        help="Comma-separated list of tools to run",
    )
    parser.add_argument("--simulate", action="store_true", help="Run on simulated data")
    parser.add_argument("--n-reads", type=int, default=1000, help="Reads to simulate")
    parser.add_argument("--sim-tl-mean", type=float, default=5000, help="Mean simulated TL")
    parser.add_argument("--sim-tl-std", type=float, default=2000, help="Std of simulated TL")
    parser.add_argument("--sim-error-rate", type=float, default=0.02, help="Error rate")
    parser.add_argument("--sim-tvr-fraction", type=float, default=0.05, help="TVR fraction")
    parser.add_argument("--parse-only", action="store_true", help="Parse existing results only")
    parser.add_argument("--wf-teloseq-dir", help="Path to wf-teloseq installation")
    parser.add_argument("-v", "--verbose", action="store_true")

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    out_dir = Path(args.output)
    out_dir.mkdir(parents=True, exist_ok=True)

    requested_tools = [t.strip() for t in args.tools.split(",")]

    # Step 1: Prepare input
    input_path = args.input
    ground_truth_path = None

    if args.simulate:
        logger.info("Generating simulated reads...")
        from simulation.simulator import SimulationParams, simulate_reads
        from simulation.simulator import write_simulated_fastq, write_ground_truth

        sim_params = SimulationParams(
            n_reads=args.n_reads,
            tl_mean=args.sim_tl_mean,
            tl_std=args.sim_tl_std,
            platform=args.platform,
            error_rate=args.sim_error_rate,
            tvr_fraction=args.sim_tvr_fraction,
        )

        reads = simulate_reads(sim_params)
        input_path = str(write_simulated_fastq(reads, out_dir / "simulated_reads.fastq"))
        ground_truth_path = str(write_ground_truth(reads, out_dir / "ground_truth.tsv"))
        logger.info("Simulated %d reads → %s", len(reads), input_path)

        # Save simulation params
        with open(out_dir / "simulation_params.json", "w") as f:
            json.dump(sim_params.__dict__, f, indent=2)

    elif not args.parse_only and not input_path:
        parser.error("--input is required unless --simulate or --parse-only is used")

    # Step 2: Run tools
    results: dict[str, ToolResult] = {}

    if not args.parse_only:
        for tool_name in requested_tools:
            if tool_name not in TOOL_REGISTRY:
                logger.warning("Unknown tool: %s", tool_name)
                continue

            wrapper = TOOL_REGISTRY[tool_name]()
            tool_out = str(out_dir / f"results_{tool_name.replace('-', '_')}")

            if not wrapper.check_installed():
                logger.warning("%s is not installed, skipping", tool_name)
                results[tool_name] = ToolResult(
                    tool_name=tool_name,
                    error_message="Not installed",
                )
                continue

            logger.info("Running %s...", tool_name)
            extra_args = {}
            if tool_name == "wf-teloseq" and args.wf_teloseq_dir:
                extra_args["workflow_dir"] = args.wf_teloseq_dir

            result = wrapper.run(
                input_path=input_path,
                output_dir=tool_out,
                platform=args.platform,
                reference=args.reference,
                threads=args.threads,
                extra_args=extra_args,
            )

            results[tool_name] = result

            if result.success:
                logger.info(
                    "%s completed in %.1f s — %d reads",
                    tool_name, result.runtime_seconds,
                    result.sample_summary.get("n_reads", 0),
                )
            else:
                logger.error("%s failed: %s", tool_name, result.error_message[:200])
    else:
        # Parse existing results
        for tool_name in requested_tools:
            if tool_name not in TOOL_REGISTRY:
                continue
            wrapper = TOOL_REGISTRY[tool_name]()
            tool_out = str(out_dir / f"results_{tool_name.replace('-', '_')}")
            if Path(tool_out).exists():
                logger.info("Parsing existing results for %s", tool_name)
                results[tool_name] = wrapper.parse_results(tool_out)

    # Step 3: Ground truth evaluation (simulation mode)
    truth_metrics = {}
    if ground_truth_path and Path(ground_truth_path).exists():
        logger.info("Evaluating against ground truth...")
        from simulation.ground_truth import evaluate_against_truth, evaluate_length_dependent_accuracy

        gt = pd.read_csv(ground_truth_path, sep="\t")

        for tool_name, result in results.items():
            if result.success and result.read_level is not None:
                m = evaluate_against_truth(result.read_level, gt, tool_name)
                truth_metrics[tool_name] = m
                logger.info(
                    "%s truth eval — sensitivity: %.1f%%, FPR: %.3f, MAE: %.0f bp, r²: %.3f",
                    tool_name,
                    m.get("sensitivity", 0) * 100,
                    m.get("false_positive_rate", 0),
                    m.get("mae", 0),
                    m.get("r_squared", 0),
                )

    # Step 4: Pairwise comparisons
    logger.info("Computing pairwise comparisons...")
    successful_tools = {k: v for k, v in results.items() if v.success}
    tool_names = list(successful_tools.keys())

    pairwise_read = {}
    pairwise_arm = {}
    pairwise_sample = {}

    for ta, tb in combinations(tool_names, 2):
        pair = f"{ta}_vs_{tb}"

        pairwise_read[pair] = compare_read_level(
            successful_tools[ta].read_level,
            successful_tools[tb].read_level,
            ta, tb,
        )

        pairwise_arm[pair] = compare_arm_level(
            successful_tools[ta].arm_level,
            successful_tools[tb].arm_level,
            ta, tb,
        )

        pairwise_sample[pair] = compare_sample_level(
            successful_tools[ta].sample_summary,
            successful_tools[tb].sample_summary,
            ta, tb,
        )

    # Step 5: Generate report
    logger.info("Generating benchmark report...")
    _generate_report(
        out_dir, results, truth_metrics,
        pairwise_read, pairwise_arm, pairwise_sample,
        args,
    )

    # Step 6: Save metrics to TSV
    _save_metrics(out_dir, results, truth_metrics, pairwise_read, pairwise_arm)

    logger.info("Benchmark complete. Results in %s", out_dir)


def _generate_report(
    out_dir: Path,
    results: dict[str, ToolResult],
    truth_metrics: dict,
    pairwise_read: dict,
    pairwise_arm: dict,
    pairwise_sample: dict,
    args,
):
    """Generate the HTML benchmark report."""
    successful = {k: v for k, v in results.items() if v.success}

    # Plots
    plot_htmls = {}

    # TL distributions
    read_data = {k: v.read_level for k, v in successful.items() if v.read_level is not None}
    plot_htmls["tl_distributions"] = plot_tl_distributions(read_data) or ""

    # Arm comparison
    arm_data = {k: v.arm_level for k, v in successful.items() if v.arm_level is not None}
    plot_htmls["arm_comparison"] = plot_arm_comparison(arm_data) or ""

    # Runtime
    runtimes = {k: v.runtime_seconds for k, v in successful.items()}
    plot_htmls["runtime"] = plot_runtime_comparison(runtimes) or ""

    # Pairwise correlation and Bland-Altman
    for pair_name in pairwise_read:
        ta, tb = pair_name.split("_vs_")
        if ta in successful and tb in successful:
            rl_a = successful[ta].read_level
            rl_b = successful[tb].read_level
            plot_htmls[f"correlation_{pair_name}"] = plot_tl_correlation(rl_a, rl_b, ta, tb) or ""
            plot_htmls[f"bland_altman_{pair_name}"] = plot_bland_altman(rl_a, rl_b, ta, tb) or ""

    # Metrics heatmap
    if pairwise_read:
        read_metrics_df = pd.DataFrame([
            {"comparison": k, **{mk: mv for mk, mv in v.items() if isinstance(mv, (int, float))}}
            for k, v in pairwise_read.items()
            if "error" not in v
        ])
        plot_htmls["metrics_heatmap"] = plot_metrics_heatmap(read_metrics_df) or ""

    # Build HTML
    html = _build_report_html(results, truth_metrics, pairwise_read,
                               pairwise_sample, plot_htmls, args)
    report_path = out_dir / "benchmark_report.html"
    report_path.write_text(html)
    logger.info("Report written to %s", report_path)


def _build_report_html(results, truth_metrics, pairwise_read,
                        pairwise_sample, plot_htmls, args) -> str:
    """Build the full HTML report."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Summary table
    summary_rows = ""
    for tool, result in results.items():
        status = "PASS" if result.success else f"FAIL: {result.error_message[:80]}"
        n_reads = result.sample_summary.get("n_reads", 0) if result.success else 0
        median_tl = result.sample_summary.get("median_tl", 0) if result.success else 0
        runtime = f"{result.runtime_seconds:.1f}" if result.runtime_seconds > 0 else "N/A"
        summary_rows += f"""
        <tr>
            <td><strong>{tool}</strong></td>
            <td>{result.version}</td>
            <td>{status}</td>
            <td>{n_reads}</td>
            <td>{median_tl:.0f}</td>
            <td>{runtime}</td>
        </tr>"""

    # Truth metrics table
    truth_rows = ""
    if truth_metrics:
        for tool, m in truth_metrics.items():
            truth_rows += f"""
            <tr>
                <td><strong>{tool}</strong></td>
                <td>{m.get('sensitivity', 0)*100:.1f}%</td>
                <td>{m.get('false_positive_rate', 0)*100:.2f}%</td>
                <td>{m.get('mae', 0):.0f}</td>
                <td>{m.get('rmse', 0):.0f}</td>
                <td>{m.get('mean_bias', 0):.0f}</td>
                <td>{m.get('pearson_r', 0):.3f}</td>
                <td>{m.get('r_squared', 0):.3f}</td>
                <td>{m.get('within_500bp_pct', 0):.1f}%</td>
            </tr>"""

    # Pairwise table
    pairwise_rows = ""
    for pair, m in pairwise_read.items():
        if "error" in m:
            continue
        pairwise_rows += f"""
        <tr>
            <td>{pair}</td>
            <td>{m.get('n_shared_reads', 0)}</td>
            <td>{m.get('pearson_r', 0):.3f}</td>
            <td>{m.get('lin_ccc', 0):.3f}</td>
            <td>{m.get('mae', 0):.0f}</td>
            <td>{m.get('agree_within_500bp', 0):.1f}%</td>
            <td>{m.get('ba_mean_diff', 0):.0f}</td>
        </tr>"""

    # Assemble plots
    plot_sections = ""
    plot_labels = {
        "tl_distributions": "TL Distribution Comparison",
        "arm_comparison": "Arm-Level Comparison",
        "runtime": "Runtime",
        "metrics_heatmap": "Concordance Metrics",
    }
    for key, html in plot_htmls.items():
        if not html:
            continue
        label = plot_labels.get(key, key.replace("_", " ").title())
        plot_sections += f'<div class="plot-container"><h3>{label}</h3>{html}</div>\n'

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>TeloPROBE Benchmark Report</title>
<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
       max-width: 1300px; margin: 0 auto; padding: 20px; background: #f8f9fa; }}
h1 {{ color: #2c3e50; border-bottom: 3px solid #e74c3c; padding-bottom: 10px; }}
h2 {{ color: #34495e; margin-top: 30px; }}
table {{ width: 100%; border-collapse: collapse; margin: 15px 0; background: white;
         border-radius: 8px; overflow: hidden; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
th, td {{ padding: 10px 15px; text-align: left; border-bottom: 1px solid #ecf0f1; }}
th {{ background: #2c3e50; color: white; }}
tr:hover {{ background: #f5f6fa; }}
.plot-container {{ background: white; padding: 20px; border-radius: 8px;
                   box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin: 15px 0; }}
.pass {{ color: #27ae60; font-weight: bold; }}
.fail {{ color: #e74c3c; }}
.meta {{ color: #7f8c8d; font-size: 0.9em; }}
.footer {{ text-align: center; color: #95a5a6; margin-top: 40px; padding-top: 20px;
           border-top: 1px solid #ecf0f1; }}
</style>
</head>
<body>
<h1>TeloPROBE Benchmark Report</h1>
<p class="meta">Generated: {timestamp} | Platform: {args.platform} |
   Input: {args.input or 'simulated'} | Threads: {args.threads}</p>

<h2>1. Tool Summary</h2>
<table>
<tr><th>Tool</th><th>Version</th><th>Status</th><th>Reads</th><th>Median TL (bp)</th><th>Runtime (s)</th></tr>
{summary_rows}
</table>

{"<h2>2. Ground Truth Evaluation</h2>" if truth_rows else ""}
{"<table><tr><th>Tool</th><th>Sensitivity</th><th>FPR</th><th>MAE (bp)</th><th>RMSE (bp)</th><th>Bias (bp)</th><th>Pearson r</th><th>R²</th><th>Within 500bp</th></tr>" + truth_rows + "</table>" if truth_rows else ""}

<h2>{"3" if truth_rows else "2"}. Pairwise Read-Level Concordance</h2>
{"<table><tr><th>Pair</th><th>Shared Reads</th><th>Pearson r</th><th>Lin CCC</th><th>MAE (bp)</th><th>Within 500bp</th><th>BA Mean Diff</th></tr>" + pairwise_rows + "</table>" if pairwise_rows else "<p>No pairwise comparisons available (need >= 2 successful tools).</p>"}

<h2>{"4" if truth_rows else "3"}. Visualizations</h2>
{plot_sections}

<div class="footer">
<p>Generated by TeloPROBE Benchmark Suite v2.0.0</p>
</div>
</body>
</html>"""


def _save_metrics(out_dir, results, truth_metrics, pairwise_read, pairwise_arm):
    """Save all metrics to TSV files for downstream analysis."""
    # Tool summary
    tool_rows = []
    for tool, result in results.items():
        row = {"tool": tool, "success": result.success,
               "version": result.version, "runtime_s": result.runtime_seconds}
        row.update(result.sample_summary)
        tool_rows.append(row)
    pd.DataFrame(tool_rows).to_csv(out_dir / "tool_summary.tsv", sep="\t", index=False)

    # Truth metrics
    if truth_metrics:
        pd.DataFrame(truth_metrics.values()).to_csv(
            out_dir / "truth_evaluation.tsv", sep="\t", index=False
        )

    # Pairwise read metrics
    if pairwise_read:
        rows = [{"pair": k, **{mk: mv for mk, mv in v.items() if not isinstance(mv, dict)}}
                for k, v in pairwise_read.items()]
        pd.DataFrame(rows).to_csv(out_dir / "pairwise_read_metrics.tsv", sep="\t", index=False)

    # Pairwise arm metrics
    if pairwise_arm:
        rows = [{"pair": k, **{mk: mv for mk, mv in v.items() if not isinstance(mv, dict)}}
                for k, v in pairwise_arm.items()]
        pd.DataFrame(rows).to_csv(out_dir / "pairwise_arm_metrics.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
