"""Per-chromosome-arm telomere length aggregation with bootstrap CIs."""

import logging
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats

from ..models.arm_estimate import ArmEstimate
from .confidence import bootstrap_ci

logger = logging.getLogger(__name__)


def compute_arm_level(
    read_level: pd.DataFrame,
    bootstrap_n: int = 1000,
    bootstrap_alpha: float = 0.05,
    min_reads: int = 5,
) -> pd.DataFrame:
    """Aggregate read-level data to chromosome-arm level.

    Groups by (chromosome_arm, haplotype, allele_cluster) and computes
    summary statistics with bootstrap confidence intervals.

    Args:
        read_level: DataFrame from compute_read_level.
        bootstrap_n: Number of bootstrap resamples.
        bootstrap_alpha: Significance level for CIs.
        min_reads: Minimum reads per group to include.

    Returns:
        DataFrame with one row per arm/haplotype/allele combination.
    """
    # Filter to good reads with arm assignments
    good = read_level[
        (read_level["qc_status"] == "Good") &
        (read_level["chromosome_arm"] != "") &
        (read_level["telomere_length_bp"] > 0)
    ].copy()

    if good.empty:
        logger.warning("No reads with arm assignments for arm-level aggregation")
        return pd.DataFrame(columns=[
            "sample", "chromosome_arm", "haplotype", "allele_id",
            "n_reads", "median_tl", "mean_tl", "trimmed_mean_tl",
            "ci_lower", "ci_upper", "cv", "anchor_confidence",
        ])

    # Group by arm, haplotype, allele
    group_cols = ["chromosome_arm", "haplotype"]
    if "allele_cluster" in good.columns:
        group_cols.append("allele_cluster")

    rows = []
    for name, group in good.groupby(group_cols, dropna=False):
        if len(group) < min_reads:
            continue

        tl = group["telomere_length_bp"].values.astype(float)
        sample_id = group["platform"].iloc[0]  # Use sample if available

        # Get sample_id from the data
        if "sample" in group.columns:
            sample_id = group["sample"].iloc[0]

        estimate = _aggregate_arm(
            tl, name, sample_id, bootstrap_n, bootstrap_alpha
        )
        rows.append(estimate.to_row())

    df = pd.DataFrame(rows)
    logger.info("Computed arm-level estimates for %d groups", len(df))
    return df


def _aggregate_arm(
    tl_values: np.ndarray,
    group_name: tuple,
    sample_id: str,
    bootstrap_n: int,
    bootstrap_alpha: float,
) -> ArmEstimate:
    """Compute summary statistics for a single arm group."""
    chr_arm = group_name[0] if len(group_name) > 0 else ""
    haplotype = group_name[1] if len(group_name) > 1 else "unknown"
    allele_id = group_name[2] if len(group_name) > 2 else None

    median_tl = float(np.median(tl_values))
    mean_tl = float(np.mean(tl_values))
    std_tl = float(np.std(tl_values))
    cv = std_tl / mean_tl if mean_tl > 0 else 0.0

    # Trimmed mean (10% trim)
    trimmed = float(stats.trim_mean(tl_values, 0.1))

    # Bootstrap CI for median
    ci_low, ci_high = bootstrap_ci(
        tl_values,
        statistic=np.median,
        n_boot=bootstrap_n,
        alpha=bootstrap_alpha,
    )

    # Quartiles
    q1 = float(np.percentile(tl_values, 25))
    q3 = float(np.percentile(tl_values, 75))

    return ArmEstimate(
        sample_id=sample_id,
        chromosome_arm=chr_arm,
        haplotype=haplotype,
        allele_id=int(allele_id) if allele_id is not None else None,
        n_reads=len(tl_values),
        median_tl=median_tl,
        mean_tl=mean_tl,
        trimmed_mean_tl=trimmed,
        std_tl=std_tl,
        cv=cv,
        q1=q1,
        q3=q3,
        ci_lower=ci_low,
        ci_upper=ci_high,
    )


def boxplot_stats(values: np.ndarray) -> dict:
    """Compute boxplot statistics for a distribution."""
    if len(values) == 0:
        return {"min": 0, "q1": 0, "median": 0, "q3": 0, "max": 0}

    q1 = float(np.percentile(values, 25))
    q3 = float(np.percentile(values, 75))
    iqr = q3 - q1
    whisker_low = float(max(np.min(values), q1 - 1.5 * iqr))
    whisker_high = float(min(np.max(values), q3 + 1.5 * iqr))

    return {
        "min": whisker_low,
        "q1": q1,
        "median": float(np.median(values)),
        "q3": q3,
        "max": whisker_high,
    }
