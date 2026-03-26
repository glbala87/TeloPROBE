"""Quantitative comparison metrics between tool outputs.

Computes concordance, correlation, and error metrics at read-level,
arm-level, and sample-level between pairs of tools.
"""

import logging
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


def compare_read_level(
    result_a: pd.DataFrame,
    result_b: pd.DataFrame,
    tool_a: str = "A",
    tool_b: str = "B",
) -> dict:
    """Compare per-read telomere lengths between two tools.

    Joins on read_id and computes concordance metrics.
    Returns dict of metrics.
    """
    if result_a is None or result_b is None:
        return {"error": "Missing read-level data"}

    # Join on read_id
    merged = pd.merge(
        result_a[["read_id", "telomere_length"]],
        result_b[["read_id", "telomere_length"]],
        on="read_id",
        suffixes=(f"_{tool_a}", f"_{tool_b}"),
        how="inner",
    )

    if merged.empty:
        return {
            "n_shared_reads": 0,
            "n_reads_a": len(result_a),
            "n_reads_b": len(result_b),
            "error": "No shared read IDs",
        }

    tl_a = merged[f"telomere_length_{tool_a}"].values.astype(float)
    tl_b = merged[f"telomere_length_{tool_b}"].values.astype(float)

    metrics = {
        "n_shared_reads": len(merged),
        "n_reads_a": len(result_a),
        "n_reads_b": len(result_b),
        "jaccard_reads": len(merged) / max(1, len(set(result_a["read_id"]) | set(result_b["read_id"]))),
    }

    # Correlation
    if len(tl_a) >= 3:
        r, p = stats.pearsonr(tl_a, tl_b)
        metrics["pearson_r"] = float(r)
        metrics["pearson_p"] = float(p)

        rho, rho_p = stats.spearmanr(tl_a, tl_b)
        metrics["spearman_rho"] = float(rho)
        metrics["spearman_p"] = float(rho_p)

    # Error metrics
    diff = tl_a - tl_b
    abs_diff = np.abs(diff)

    metrics["mean_diff"] = float(np.mean(diff))
    metrics["median_diff"] = float(np.median(diff))
    metrics["mae"] = float(np.mean(abs_diff))
    metrics["rmse"] = float(np.sqrt(np.mean(diff ** 2)))
    metrics["median_abs_diff"] = float(np.median(abs_diff))

    # Relative error
    mean_tl = (tl_a + tl_b) / 2
    valid = mean_tl > 0
    if valid.sum() > 0:
        metrics["mape"] = float(np.mean(abs_diff[valid] / mean_tl[valid]) * 100)

    # Concordance correlation coefficient (Lin's CCC)
    metrics["lin_ccc"] = _lins_ccc(tl_a, tl_b)

    # Bland-Altman metrics
    ba = bland_altman(tl_a, tl_b)
    metrics["ba_mean_diff"] = ba["mean_diff"]
    metrics["ba_lower_loa"] = ba["lower_loa"]
    metrics["ba_upper_loa"] = ba["upper_loa"]

    # Agreement within tolerance
    for tol in [100, 500, 1000]:
        metrics[f"agree_within_{tol}bp"] = float(
            (abs_diff <= tol).mean() * 100
        )

    return metrics


def compare_arm_level(
    result_a: pd.DataFrame,
    result_b: pd.DataFrame,
    tool_a: str = "A",
    tool_b: str = "B",
) -> dict:
    """Compare per-arm telomere length estimates between two tools."""
    if result_a is None or result_b is None or result_a.empty or result_b.empty:
        return {"error": "Missing arm-level data"}

    # Normalize arm names and merge
    a = result_a.copy()
    b = result_b.copy()
    a["arm_key"] = a["chromosome_arm"].str.lower().str.replace("_", "")
    b["arm_key"] = b["chromosome_arm"].str.lower().str.replace("_", "")

    merged = pd.merge(
        a[["arm_key", "median_tl", "n_reads"]],
        b[["arm_key", "median_tl", "n_reads"]],
        on="arm_key",
        suffixes=(f"_{tool_a}", f"_{tool_b}"),
        how="inner",
    )

    if merged.empty:
        return {
            "n_shared_arms": 0,
            "n_arms_a": len(a),
            "n_arms_b": len(b),
            "error": "No shared chromosome arms",
        }

    tl_a = merged[f"median_tl_{tool_a}"].values.astype(float)
    tl_b = merged[f"median_tl_{tool_b}"].values.astype(float)

    metrics = {
        "n_shared_arms": len(merged),
        "n_arms_a": a["arm_key"].nunique(),
        "n_arms_b": b["arm_key"].nunique(),
        "arm_jaccard": len(merged) / max(1, len(set(a["arm_key"]) | set(b["arm_key"]))),
    }

    if len(tl_a) >= 3:
        r, p = stats.pearsonr(tl_a, tl_b)
        metrics["pearson_r"] = float(r)
        rho, _ = stats.spearmanr(tl_a, tl_b)
        metrics["spearman_rho"] = float(rho)

    diff = tl_a - tl_b
    metrics["mae"] = float(np.mean(np.abs(diff)))
    metrics["rmse"] = float(np.sqrt(np.mean(diff ** 2)))
    metrics["lin_ccc"] = _lins_ccc(tl_a, tl_b)

    return metrics


def compare_sample_level(
    summary_a: dict,
    summary_b: dict,
    tool_a: str = "A",
    tool_b: str = "B",
) -> dict:
    """Compare sample-level summary metrics."""
    metrics = {}

    for key in ["median_tl", "mean_tl", "n_reads"]:
        val_a = summary_a.get(key, 0)
        val_b = summary_b.get(key, 0)
        metrics[f"{key}_{tool_a}"] = val_a
        metrics[f"{key}_{tool_b}"] = val_b
        if val_a and val_b:
            metrics[f"{key}_diff"] = float(val_a) - float(val_b)
            avg = (float(val_a) + float(val_b)) / 2
            if avg > 0:
                metrics[f"{key}_pct_diff"] = (
                    abs(float(val_a) - float(val_b)) / avg * 100
                )

    return metrics


def bland_altman(
    x: np.ndarray,
    y: np.ndarray,
) -> dict:
    """Compute Bland-Altman agreement statistics.

    Returns mean difference and 95% limits of agreement.
    """
    diff = x - y
    mean_diff = float(np.mean(diff))
    std_diff = float(np.std(diff, ddof=1)) if len(diff) > 1 else 0.0

    return {
        "mean_diff": mean_diff,
        "std_diff": std_diff,
        "lower_loa": mean_diff - 1.96 * std_diff,
        "upper_loa": mean_diff + 1.96 * std_diff,
    }


def _lins_ccc(x: np.ndarray, y: np.ndarray) -> float:
    """Lin's Concordance Correlation Coefficient.

    Measures agreement between two continuous measurements,
    combining precision (Pearson r) and accuracy (bias correction).
    CCC = 1.0 means perfect agreement.
    """
    if len(x) < 3:
        return 0.0

    mean_x = np.mean(x)
    mean_y = np.mean(y)
    var_x = np.var(x, ddof=1)
    var_y = np.var(y, ddof=1)
    cov_xy = np.cov(x, y, ddof=1)[0, 1]

    numerator = 2 * cov_xy
    denominator = var_x + var_y + (mean_x - mean_y) ** 2

    if denominator == 0:
        return 0.0
    return float(numerator / denominator)


def rank_tools(
    all_comparisons: dict[str, dict],
    primary_metric: str = "lin_ccc",
) -> pd.DataFrame:
    """Rank tools by a primary concordance metric.

    Args:
        all_comparisons: Dict of {tool_pair: comparison_metrics}.
        primary_metric: Metric to rank by.

    Returns:
        DataFrame with tool rankings.
    """
    rows = []
    for pair_name, metrics in all_comparisons.items():
        if "error" in metrics:
            continue
        rows.append({
            "comparison": pair_name,
            primary_metric: metrics.get(primary_metric, 0.0),
            "mae": metrics.get("mae", 0.0),
            "rmse": metrics.get("rmse", 0.0),
            "pearson_r": metrics.get("pearson_r", 0.0),
            "n_shared": metrics.get("n_shared_reads", metrics.get("n_shared_arms", 0)),
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(primary_metric, ascending=False)
    return df
