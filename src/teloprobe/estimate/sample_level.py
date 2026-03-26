"""Sample-level telomere length aggregation."""

import logging

import numpy as np
import pandas as pd
from scipy import stats

from ..constants import Platform, Mode
from ..models.sample_result import SampleResult
from .confidence import bootstrap_ci

logger = logging.getLogger(__name__)


def compute_sample_level(
    read_level: pd.DataFrame,
    arm_level: pd.DataFrame,
    sample_id: str = "sample",
    platform: Platform = Platform.ONT,
    mode: Mode = Mode.TELOSEQ,
    bootstrap_n: int = 1000,
    bootstrap_alpha: float = 0.05,
) -> pd.DataFrame:
    """Aggregate to sample-level telomere statistics.

    Args:
        read_level: Read-level DataFrame.
        arm_level: Arm-level DataFrame.
        sample_id: Sample identifier.
        platform: Sequencing platform.
        mode: Analysis mode.
        bootstrap_n: Bootstrap resamples.
        bootstrap_alpha: CI significance level.

    Returns:
        Single-row DataFrame with sample-level statistics.
    """
    total_reads = len(read_level)
    good_mask = read_level["qc_status"] == "Good"
    good_reads = read_level[good_mask]

    tl_values = good_reads["telomere_length_bp"].values.astype(float)
    tl_positive = tl_values[tl_values > 0]

    # QC summary: count by status
    qc_counts = read_level["qc_status"].value_counts().to_dict()

    result = SampleResult(
        sample_id=sample_id,
        platform=platform,
        mode=mode,
        total_reads=total_reads,
        telomeric_reads=int(good_mask.sum()),
        passing_reads=len(tl_positive),
        informative_reads=len(tl_positive),
        qc_summary=qc_counts,
    )

    if len(tl_positive) > 0:
        result.median_tl = float(np.median(tl_positive))
        result.mean_tl = float(np.mean(tl_positive))
        result.trimmed_mean_tl = float(stats.trim_mean(tl_positive, 0.1))

        std = float(np.std(tl_positive))
        result.cv = std / result.mean_tl if result.mean_tl > 0 else 0.0

        ci_low, ci_high = bootstrap_ci(
            tl_positive,
            statistic=np.median,
            n_boot=bootstrap_n,
            alpha=bootstrap_alpha,
        )
        result.ci_lower = ci_low
        result.ci_upper = ci_high

    if not arm_level.empty:
        result.arms_assigned = int(arm_level["chromosome_arm"].nunique())

    return pd.DataFrame([result.to_row()])
