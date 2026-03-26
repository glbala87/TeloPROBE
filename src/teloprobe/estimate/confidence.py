"""Confidence scoring and bootstrap confidence intervals.

Combines three evidence layers (motif, boundary, anchor) into
a single confidence score, weighted by analysis mode.
"""

import logging
from typing import Callable, Optional

import numpy as np

from ..constants import Mode, EVIDENCE_WEIGHTS
from ..config import Config

logger = logging.getLogger(__name__)


def bootstrap_ci(
    values: np.ndarray,
    statistic: Callable = np.median,
    n_boot: int = 1000,
    alpha: float = 0.05,
    random_seed: int = 42,
) -> tuple[float, float]:
    """Compute BCa (bias-corrected accelerated) bootstrap confidence interval.

    Args:
        values: Data array.
        statistic: Summary statistic function.
        n_boot: Number of bootstrap resamples.
        alpha: Significance level (default 0.05 for 95% CI).
        random_seed: For reproducibility.

    Returns:
        Tuple of (ci_lower, ci_upper).
    """
    if len(values) < 3:
        val = statistic(values) if len(values) > 0 else 0.0
        return float(val), float(val)

    rng = np.random.RandomState(random_seed)
    n = len(values)
    observed = float(statistic(values))

    # Generate bootstrap distribution
    boot_stats = np.array([
        statistic(rng.choice(values, size=n, replace=True))
        for _ in range(n_boot)
    ])

    # BCa correction: bias correction factor
    z0 = _norm_ppf(np.mean(boot_stats < observed))

    # Acceleration factor (jackknife)
    jackknife_stats = np.array([
        statistic(np.delete(values, i)) for i in range(n)
    ])
    jack_mean = jackknife_stats.mean()
    num = ((jack_mean - jackknife_stats) ** 3).sum()
    den = 6.0 * ((jack_mean - jackknife_stats) ** 2).sum() ** 1.5
    a = num / den if den != 0 else 0.0

    # Adjusted percentiles
    z_alpha_low = _norm_ppf(alpha / 2)
    z_alpha_high = _norm_ppf(1 - alpha / 2)

    alpha_low = _norm_cdf(z0 + (z0 + z_alpha_low) / (1 - a * (z0 + z_alpha_low)))
    alpha_high = _norm_cdf(z0 + (z0 + z_alpha_high) / (1 - a * (z0 + z_alpha_high)))

    # Clamp to valid percentile range
    alpha_low = max(0.001, min(0.999, alpha_low))
    alpha_high = max(0.001, min(0.999, alpha_high))

    ci_lower = float(np.percentile(boot_stats, alpha_low * 100))
    ci_upper = float(np.percentile(boot_stats, alpha_high * 100))

    return ci_lower, ci_upper


def _norm_ppf(p: float) -> float:
    """Normal distribution percent point function (inverse CDF)."""
    from scipy.stats import norm
    p = max(1e-10, min(1 - 1e-10, p))
    return float(norm.ppf(p))


def _norm_cdf(z: float) -> float:
    """Normal distribution CDF."""
    from scipy.stats import norm
    return float(norm.cdf(z))


def compute_confidence_score(
    call,  # TelomereCall
    config: Config,
) -> float:
    """Compute overall confidence score from three evidence layers.

    The score is a weighted combination of:
    - Motif evidence: motif density support
    - Boundary evidence: agreement between HMM and PELT
    - Anchor evidence: subtelomeric alignment quality

    Weights depend on analysis mode (WGS vs TeloPROBE enriched).

    Returns:
        Confidence score between 0.0 and 1.0.
    """
    weights = EVIDENCE_WEIGHTS.get(
        Mode(config.mode) if isinstance(config.mode, str) else config.mode,
        EVIDENCE_WEIGHTS[Mode.TELOSEQ],
    )

    motif_score = call.evidence.motif_score
    boundary_score = call.evidence.boundary_score
    anchor_score = call.evidence.anchor_score

    # Weighted sum
    score = (
        weights["motif"] * motif_score +
        weights["boundary"] * boundary_score +
        weights["anchor"] * anchor_score
    )

    # Penalty for failed QC
    from ..constants import QCStatus
    if call.qc_status != QCStatus.GOOD:
        score *= 0.1

    return round(max(0.0, min(1.0, score)), 4)


def evidence_weights(mode: Mode) -> dict[str, float]:
    """Return evidence layer weights for a given mode."""
    return dict(EVIDENCE_WEIGHTS.get(mode, EVIDENCE_WEIGHTS[Mode.TELOSEQ]))
