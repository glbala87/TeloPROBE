"""Sample-level QC checks and diagnostics."""

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class QCCheck:
    """A single QC check result."""
    name: str
    passed: bool
    value: float
    threshold: float
    message: str


def sample_qc_checks(
    read_level: pd.DataFrame,
    arm_level: pd.DataFrame,
    min_reads: int = 20,
    min_arm_reads: int = 10,
    max_cv: float = 1.0,
) -> list[QCCheck]:
    """Run sample-level QC checks.

    Checks:
    1. Minimum informative reads
    2. Coefficient of variation
    3. Arm assignment rate
    4. Expected chromosome arm count
    """
    checks = []
    good = read_level[read_level["qc_status"] == "Good"]
    tl = good["telomere_length_bp"].values.astype(float)
    tl_positive = tl[tl > 0]

    # 1. Minimum informative reads
    checks.append(QCCheck(
        name="min_reads",
        passed=len(tl_positive) >= min_reads,
        value=len(tl_positive),
        threshold=min_reads,
        message=f"Informative reads: {len(tl_positive)} (min: {min_reads})",
    ))

    # 2. CV check
    if len(tl_positive) > 1:
        mean_val = float(np.mean(tl_positive))
        cv = float(np.std(tl_positive) / mean_val) if mean_val > 0 else 0.0
        checks.append(QCCheck(
            name="cv",
            passed=cv <= max_cv,
            value=cv,
            threshold=max_cv,
            message=f"CV: {cv:.3f} (max: {max_cv})",
        ))

    # 3. Arm assignment rate
    if not arm_level.empty:
        assigned = good["chromosome_arm"].apply(
            lambda x: x != "" if isinstance(x, str) else False
        ).sum()
        rate = assigned / len(good) if len(good) > 0 else 0.0
        checks.append(QCCheck(
            name="arm_assignment_rate",
            passed=rate > 0.1,
            value=rate,
            threshold=0.1,
            message=f"Arm assignment rate: {rate:.1%}",
        ))

        # 4. Number of arms
        n_arms = arm_level["chromosome_arm"].nunique()
        checks.append(QCCheck(
            name="arm_count",
            passed=n_arms >= 10,
            value=n_arms,
            threshold=10,
            message=f"Chromosome arms detected: {n_arms}",
        ))

    return checks


def flag_low_coverage_arms(
    arm_level: pd.DataFrame,
    min_reads: int = 10,
) -> pd.DataFrame:
    """Flag arms with insufficient reads for reliable estimation."""
    result = arm_level.copy()
    result["low_coverage"] = result["n_reads"] < min_reads
    n_flagged = result["low_coverage"].sum()
    if n_flagged > 0:
        logger.warning(
            "%d arms have fewer than %d reads (low coverage)",
            n_flagged, min_reads,
        )
    return result


def downsampling_stability(
    tl_values: np.ndarray,
    fractions: list[float] = None,
    n_replicates: int = 10,
    random_seed: int = 42,
) -> pd.DataFrame:
    """Test stability of median TL estimate under downsampling.

    Returns DataFrame with columns: fraction, replicate, median_tl.
    """
    if fractions is None:
        fractions = [0.1, 0.25, 0.5, 0.75, 1.0]

    rng = np.random.RandomState(random_seed)
    rows = []

    for frac in fractions:
        n_sample = max(1, int(len(tl_values) * frac))
        for rep in range(n_replicates):
            subsample = rng.choice(tl_values, size=n_sample, replace=False)
            rows.append({
                "fraction": frac,
                "replicate": rep,
                "median_tl": float(np.median(subsample)),
                "n_reads": n_sample,
            })

    return pd.DataFrame(rows)
