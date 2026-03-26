"""Change-point detection for telomere boundary inference.

Uses the ruptures library for statistically principled change-point
detection (PELT, BinSeg) on the motif density signal.

This replaces the convolution-based edge filter from the original
wf-teloprobe, providing:
- Automatic penalty selection via BIC
- Multiple change-point detection
- Statistical robustness to noise
"""

import logging
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

try:
    import ruptures as rpt
    HAS_RUPTURES = True
except ImportError:
    HAS_RUPTURES = False
    logger.warning("ruptures not installed; change-point detection disabled")


def detect_changepoints(
    density_vector: np.ndarray,
    method: str = "pelt",
    penalty: Optional[float] = None,
    min_size: int = 60,
    n_bkps: int = 5,
) -> list[int]:
    """Detect change-points in the motif density signal.

    Args:
        density_vector: Per-base motif density from compute_motif_density.
        method: Detection method - "pelt" (default), "binseg", or "bottomup".
        penalty: Penalty value for PELT. None = auto-select via BIC.
        min_size: Minimum segment length between changepoints.
        n_bkps: Maximum number of breakpoints (for binseg/bottomup).

    Returns:
        List of changepoint positions (in base-pair indices).
        Last element is always the signal length (ruptures convention).
    """
    if not HAS_RUPTURES:
        return _fallback_changepoint(density_vector, min_size)

    signal = density_vector.astype(np.float64)

    if len(signal) < min_size * 2:
        return [len(signal)]

    try:
        if method == "pelt":
            if penalty is None:
                penalty = auto_penalty(signal)
            algo = rpt.Pelt(model="rbf", min_size=min_size).fit(signal)
            bkps = algo.predict(pen=penalty)
        elif method == "binseg":
            algo = rpt.Binseg(model="rbf", min_size=min_size).fit(signal)
            bkps = algo.predict(n_bkps=n_bkps)
        elif method == "bottomup":
            algo = rpt.BottomUp(model="rbf", min_size=min_size).fit(signal)
            bkps = algo.predict(n_bkps=n_bkps)
        else:
            raise ValueError(f"Unknown changepoint method: {method}")

        return bkps

    except Exception as e:
        logger.warning("Change-point detection failed (%s): %s", method, e)
        return _fallback_changepoint(density_vector, min_size)


def auto_penalty(density_vector: np.ndarray) -> float:
    """Automatically select PELT penalty using BIC-like criterion.

    The penalty controls sensitivity: higher = fewer changepoints.
    We scale it with signal variance and length.
    """
    n = len(density_vector)
    var = np.var(density_vector)

    # BIC-inspired penalty: log(n) * variance * scaling factor
    penalty = np.log(n) * max(var, 0.01) * 2.0

    # Clamp to reasonable range
    return max(0.1, min(penalty, 10.0))


def select_boundary_changepoint(
    changepoints: list[int],
    density: np.ndarray,
    min_density_drop: float = 0.3,
) -> Optional[int]:
    """From detected changepoints, select the telomere-subtelomere boundary.

    The boundary is the changepoint with the largest density drop
    (from high to low), consistent with a telomere ending.

    Args:
        changepoints: List of changepoint positions from detect_changepoints.
        density: Motif density array.
        min_density_drop: Minimum density difference to qualify as boundary.

    Returns:
        The selected boundary position, or None if no valid boundary found.
    """
    if not changepoints or len(density) == 0:
        return None

    # Remove the terminal position (signal length)
    interior_bkps = [cp for cp in changepoints if cp < len(density)]
    if not interior_bkps:
        return None

    best_pos = None
    best_drop = 0.0

    for cp in interior_bkps:
        # Compute density before and after the changepoint
        window = min(60, cp, len(density) - cp)
        if window < 10:
            continue

        before = np.mean(density[max(0, cp - window):cp])
        after = np.mean(density[cp:min(len(density), cp + window)])
        drop = before - after

        if drop > best_drop and drop >= min_density_drop:
            best_drop = drop
            best_pos = cp

    return best_pos


def _fallback_changepoint(
    density_vector: np.ndarray,
    min_size: int = 60,
) -> list[int]:
    """Simple threshold-based fallback when ruptures is not available.

    Finds the position where density drops below 0.5 for the first time
    after an initial high-density region.
    """
    n = len(density_vector)
    if n < min_size * 2:
        return [n]

    # Look for sustained drop below threshold
    threshold = 0.5
    in_telo = True
    boundary = None

    for i in range(min_size, n - min_size):
        window = density_vector[i:i + min_size // 2]
        if in_telo and np.mean(window) < threshold:
            boundary = i
            break

    if boundary is not None:
        return [boundary, n]
    return [n]


def detect_changepoints_bayesian(
    density_vector: np.ndarray,
    min_size: int = 60,
    prior_scale: float = 1.0,
) -> list[int]:
    """Bayesian online change-point detection (optional alternative).

    Uses a simple Bayesian approach where each position is evaluated
    as a potential changepoint based on the likelihood ratio of
    the data being from two different distributions.
    """
    n = len(density_vector)
    if n < min_size * 2:
        return [n]

    # Compute cumulative statistics
    cumsum = np.cumsum(density_vector)
    cumsum2 = np.cumsum(density_vector ** 2)

    log_likelihood_ratios = np.zeros(n)

    for t in range(min_size, n - min_size):
        # Before segment: [0, t)
        n1 = t
        mean1 = cumsum[t - 1] / n1
        var1 = max(cumsum2[t - 1] / n1 - mean1 ** 2, 1e-6)

        # After segment: [t, n)
        n2 = n - t
        mean2 = (cumsum[n - 1] - cumsum[t - 1]) / n2
        var2_sum = cumsum2[n - 1] - cumsum2[t - 1]
        var2 = max(var2_sum / n2 - mean2 ** 2, 1e-6)

        # Full segment variance
        mean_all = cumsum[n - 1] / n
        var_all = max(cumsum2[n - 1] / n - mean_all ** 2, 1e-6)

        # Log-likelihood ratio
        llr = (n / 2) * np.log(var_all) - (n1 / 2) * np.log(var1) - (n2 / 2) * np.log(var2)
        log_likelihood_ratios[t] = llr

    # Find peaks above threshold
    threshold = prior_scale * np.log(n)
    peaks = []
    for i in range(min_size, n - min_size):
        if (log_likelihood_ratios[i] > threshold and
            log_likelihood_ratios[i] >= log_likelihood_ratios[max(0, i - 1)] and
            log_likelihood_ratios[i] >= log_likelihood_ratios[min(n - 1, i + 1)]):
            peaks.append(i)

    peaks.append(n)
    return peaks
