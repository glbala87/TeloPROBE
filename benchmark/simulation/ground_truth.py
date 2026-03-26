"""Ground truth evaluation against simulated data."""

import logging

import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)


def evaluate_against_truth(
    tool_result: pd.DataFrame,
    ground_truth: pd.DataFrame,
    tool_name: str = "",
) -> dict:
    """Compare tool's telomere length estimates against known ground truth.

    Args:
        tool_result: DataFrame with columns [read_id, telomere_length].
        ground_truth: DataFrame with columns [read_id, true_tl, read_type].
        tool_name: Name of the tool for labeling.

    Returns:
        Dict of accuracy metrics.
    """
    if tool_result is None or tool_result.empty:
        return {"tool": tool_name, "error": "No results"}

    # Merge on read_id
    merged = pd.merge(
        ground_truth, tool_result[["read_id", "telomere_length"]],
        on="read_id", how="left",
    )

    # Split by read type
    telo = merged[merged["read_type"] == "telomeric"].copy()
    its = merged[merged["read_type"] == "its"].copy()
    random_reads = merged[merged["read_type"] == "random"].copy()

    metrics = {"tool": tool_name}

    # --- Detection sensitivity ---
    telo_detected = telo["telomere_length"].notna() & (telo["telomere_length"] > 0)
    metrics["sensitivity"] = float(telo_detected.mean()) if len(telo) > 0 else 0.0
    metrics["n_true_telo"] = len(telo)
    metrics["n_detected_telo"] = int(telo_detected.sum())

    # --- False positive rate (ITS and random incorrectly called telomeric) ---
    its_fp = its["telomere_length"].notna() & (its["telomere_length"] > 0)
    random_fp = random_reads["telomere_length"].notna() & (random_reads["telomere_length"] > 0)
    total_neg = len(its) + len(random_reads)
    total_fp = int(its_fp.sum()) + int(random_fp.sum())
    metrics["false_positive_rate"] = total_fp / max(1, total_neg)
    metrics["its_false_positives"] = int(its_fp.sum())
    metrics["random_false_positives"] = int(random_fp.sum())

    # --- TL accuracy (detected telomeric reads only) ---
    detected = telo[telo_detected].copy()
    if len(detected) >= 3:
        true_tl = detected["true_tl"].values.astype(float)
        est_tl = detected["telomere_length"].values.astype(float)
        diff = est_tl - true_tl

        metrics["mae"] = float(np.mean(np.abs(diff)))
        metrics["rmse"] = float(np.sqrt(np.mean(diff ** 2)))
        metrics["mean_bias"] = float(np.mean(diff))
        metrics["median_bias"] = float(np.median(diff))

        # Relative error
        metrics["mape"] = float(np.mean(np.abs(diff) / np.maximum(true_tl, 1)) * 100)

        # Correlation
        r, p = stats.pearsonr(true_tl, est_tl)
        metrics["pearson_r"] = float(r)
        metrics["pearson_p"] = float(p)

        rho, _ = stats.spearmanr(true_tl, est_tl)
        metrics["spearman_rho"] = float(rho)

        # Agreement bins
        abs_diff = np.abs(diff)
        for tol in [100, 500, 1000, 2000]:
            metrics[f"within_{tol}bp_pct"] = float((abs_diff <= tol).mean() * 100)

        # R-squared
        ss_res = np.sum(diff ** 2)
        ss_tot = np.sum((true_tl - true_tl.mean()) ** 2)
        metrics["r_squared"] = float(1 - ss_res / ss_tot) if ss_tot > 0 else 0.0

    return metrics


def evaluate_length_dependent_accuracy(
    tool_result: pd.DataFrame,
    ground_truth: pd.DataFrame,
    bins: list[tuple[float, float]] = None,
) -> pd.DataFrame:
    """Evaluate accuracy stratified by true telomere length bins.

    Returns DataFrame with accuracy metrics per length bin.
    """
    if bins is None:
        bins = [
            (0, 1000), (1000, 3000), (3000, 5000),
            (5000, 8000), (8000, 12000), (12000, 20000),
        ]

    merged = pd.merge(
        ground_truth[ground_truth["read_type"] == "telomeric"],
        tool_result[["read_id", "telomere_length"]],
        on="read_id", how="inner",
    )

    rows = []
    for low, high in bins:
        mask = (merged["true_tl"] >= low) & (merged["true_tl"] < high)
        subset = merged[mask]

        if len(subset) < 3:
            continue

        diff = subset["telomere_length"].values - subset["true_tl"].values
        rows.append({
            "bin": f"{low}-{high}",
            "n_reads": len(subset),
            "mae": float(np.mean(np.abs(diff))),
            "mean_bias": float(np.mean(diff)),
            "rmse": float(np.sqrt(np.mean(diff ** 2))),
            "mape": float(np.mean(np.abs(diff) / np.maximum(subset["true_tl"].values, 1)) * 100),
        })

    return pd.DataFrame(rows)
