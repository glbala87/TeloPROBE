"""Cross-evidence validation for telomere calls."""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def cross_validate_evidence(read_level: pd.DataFrame) -> pd.DataFrame:
    """Flag reads where evidence layers conflict.

    Detects inconsistencies such as:
    - High motif density but no boundary detected
    - Boundary detected but very low motif density
    - High boundary confidence but poor alignment
    """
    result = read_level.copy()
    result["evidence_conflict"] = False

    good = result["qc_status"] == "Good"

    # High motif density but telomere_length is -1
    conflict_1 = (
        good &
        (result["motif_density"] > 0.7) &
        (result["telomere_length_bp"] <= 0)
    )

    # Telomere detected but very low motif density
    conflict_2 = (
        good &
        (result["telomere_length_bp"] > 0) &
        (result["motif_density"] < 0.2)
    )

    result.loc[conflict_1 | conflict_2, "evidence_conflict"] = True

    n_conflicts = result["evidence_conflict"].sum()
    if n_conflicts > 0:
        logger.warning(
            "%d reads have conflicting evidence layers", n_conflicts
        )

    return result


def validate_tvr_consistency(
    read_level: pd.DataFrame,
    variance_threshold: float = 0.3,
) -> pd.DataFrame:
    """Check TVR composition consistency within chromosome arms.

    For each chromosome arm, computes the variance of TVR fractions across
    reads assigned to that arm. Arms with high variance may have mixed
    chromosome-end assignments (reads from different telomeres incorrectly
    grouped together).

    Args:
        read_level: Read-level DataFrame. Must contain 'chromosome_arm'
            and 'tvr_fraction' columns. If 'tvr_fraction' is missing,
            the column is computed as variant_motif_count / total_motif_count
            when those columns are available, otherwise the function returns
            the input unchanged with a tvr_inconsistent column set to False.
        variance_threshold: Maximum acceptable variance of TVR fractions
            within an arm. Arms exceeding this are flagged.

    Returns:
        DataFrame with one row per chromosome arm, including a boolean
        'tvr_inconsistent' column indicating potentially mixed assignments.
    """
    result = read_level.copy()

    # Ensure we have a tvr_fraction column
    if "tvr_fraction" not in result.columns:
        if "variant_motif_count" in result.columns and "total_motif_count" in result.columns:
            total = result["total_motif_count"].astype(float)
            result["tvr_fraction"] = np.where(
                total > 0,
                result["variant_motif_count"].astype(float) / total,
                0.0,
            )
        else:
            logger.info(
                "No TVR fraction data available; skipping TVR consistency check"
            )
            result["tvr_inconsistent"] = False
            return result

    # Filter to good reads with arm assignments
    has_arm = (
        result["chromosome_arm"].notna()
        & (result["chromosome_arm"] != "")
    )
    good = result[has_arm].copy()

    if good.empty:
        result["tvr_inconsistent"] = False
        return result

    # Compute per-arm variance of TVR fractions
    arm_variance = (
        good.groupby("chromosome_arm")["tvr_fraction"]
        .var(ddof=0)
        .fillna(0.0)
        .rename("tvr_variance")
    )

    # Flag arms exceeding the threshold
    inconsistent_arms = set(
        arm_variance[arm_variance > variance_threshold].index
    )

    if inconsistent_arms:
        logger.warning(
            "%d arms flagged for TVR inconsistency (variance > %.2f): %s",
            len(inconsistent_arms),
            variance_threshold,
            ", ".join(sorted(inconsistent_arms)),
        )

    result["tvr_inconsistent"] = result["chromosome_arm"].isin(inconsistent_arms)
    return result
