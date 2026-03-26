"""Consensus boundary detection combining HMM and change-point evidence.

This module merges the results from the HMM labeler and PELT change-point
detector to produce a robust boundary call with confidence scoring.
"""

import logging
from typing import Optional

import numpy as np

from ..constants import QCStatus, MOTIF_LEN
from ..config import Config

logger = logging.getLogger(__name__)


def consensus_boundary(
    hmm_bnd: Optional[int],
    pelt_bnd: Optional[int],
    motif_array: np.ndarray,
    tolerance: int = 60,
) -> tuple[Optional[int], str, float]:
    """Merge HMM and PELT boundary estimates into a consensus.

    Args:
        hmm_bnd: Boundary from HMM Viterbi decoding (or None).
        pelt_bnd: Boundary from PELT change-point detection (or None).
        motif_array: Binary motif hit array for additional validation.
        tolerance: Maximum distance (bp) for boundaries to be considered agreeing.

    Returns:
        Tuple of (boundary_position, method_used, agreement_score).
        - agreement_score: 1.0 if both agree, lower if single-source or disagreement.
    """
    if hmm_bnd is not None and pelt_bnd is not None:
        distance = abs(hmm_bnd - pelt_bnd)

        if distance <= tolerance:
            # Both agree — high confidence, use the mean
            boundary = int((hmm_bnd + pelt_bnd) / 2)
            agreement = 1.0 - (distance / tolerance) * 0.3
            return boundary, "consensus", agreement
        else:
            # Disagree — use the one with better density contrast
            hmm_contrast = _density_contrast(motif_array, hmm_bnd)
            pelt_contrast = _density_contrast(motif_array, pelt_bnd)

            if hmm_contrast >= pelt_contrast:
                return hmm_bnd, "hmm", 0.5
            else:
                return pelt_bnd, "pelt", 0.5

    elif hmm_bnd is not None:
        return hmm_bnd, "hmm", 0.6

    elif pelt_bnd is not None:
        return pelt_bnd, "pelt", 0.6

    else:
        # Neither method found a boundary — try simple density threshold
        fallback = _threshold_boundary(motif_array)
        if fallback is not None:
            return fallback, "threshold", 0.3
        return None, "none", 0.0


def _density_contrast(motif_array: np.ndarray, boundary: int, window: int = 60) -> float:
    """Compute density contrast across a boundary position."""
    if boundary <= window or boundary >= len(motif_array) - window:
        return 0.0

    before = motif_array[boundary - window:boundary].mean()
    after = motif_array[boundary:boundary + window].mean()
    return before - after


def _threshold_boundary(
    motif_array: np.ndarray,
    threshold: float = 0.5,
    min_telo_length: int = 120,
) -> Optional[int]:
    """Simple threshold-based boundary as last resort.

    Finds the first position after an initial telomeric block where
    motif density drops below threshold.
    """
    from ..candidate.motif_scanner import compute_motif_density_fast

    if len(motif_array) < min_telo_length:
        return None

    density = compute_motif_density_fast(motif_array, window_size=60)

    # Check that the start is telomeric
    if density[:min_telo_length].mean() < threshold:
        return None

    # Find where density drops below threshold
    below = np.where(density < threshold)[0]
    candidates = below[below >= min_telo_length]

    if len(candidates) > 0:
        return int(candidates[0])
    return None


def classify_boundary(
    boundary: Optional[int],
    read_length: int,
    motif_array: np.ndarray,
    qualities: Optional[np.ndarray],
    motif_count: int,
    config: Config,
) -> QCStatus:
    """Apply the QC classification cascade to a boundary call.

    This implements the 9-filter cascade from the original wf-teloprobe,
    now informed by HMM state path when available.

    Args:
        boundary: Detected boundary position.
        read_length: Total read length.
        motif_array: Binary motif hit array.
        qualities: Per-base quality array (or None).
        motif_count: Total number of motif matches.
        config: Pipeline configuration.

    Returns:
        QCStatus indicating pass/fail and reason.
    """
    filter_width_bp = config.filter_width * MOTIF_LEN

    # 1. Read too short
    if read_length < config.min_read_length_boundary:
        return QCStatus.TOO_SHORT

    # 2. Too few repeats
    if motif_count < config.min_repeats:
        return QCStatus.TOO_FEW_REPEATS

    # 3. No boundary detected
    if boundary is None:
        return QCStatus.NO_BOUNDARY

    # 4. Boundary too close to start
    if boundary <= filter_width_bp:
        return QCStatus.TOO_CLOSE_START

    # 5. Boundary too close to end
    if boundary + filter_width_bp // 2 > read_length:
        return QCStatus.TOO_CLOSE_END

    # 6. Start not repeats — first portion must be mostly telomeric
    from ..candidate.read_filter import check_start_repeats
    if check_start_repeats(
        motif_array, boundary,
        config.start_window_frac, config.start_repeats_frac
    ):
        return QCStatus.START_NOT_REPEATS

    # 7. Low subtelomeric quality
    from ..candidate.read_filter import check_post_boundary_quality
    if check_post_boundary_quality(qualities, boundary, config.min_qual_non_telo):
        return QCStatus.LOW_SUBTELO_QUAL

    # 8. Telomere only (no subtelomeric flank)
    from ..candidate.read_filter import check_post_boundary_ccc
    sequence_after = ""  # Need full sequence for this check
    # This check will be done at a higher level where sequence is available

    # 9. Error motif check done at higher level

    return QCStatus.GOOD


def determine_telomere_side(
    motif_array: np.ndarray,
    boundary: Optional[int],
    read_length: int,
    min_terminal_fraction: float = 0.6,
) -> str:
    """Determine which end of the read has the telomeric tract.

    For vertebrate telomeres, strand/orientation patterns help
    discriminate genuine terminal telomeres from internal repeats.

    Returns: "5prime", "3prime", "both", or "internal"
    """
    if boundary is None or len(motif_array) == 0:
        return "internal"

    # Check 5' end (start of read)
    five_prime_density = motif_array[:min(boundary, 120)].mean() if boundary > 0 else 0.0

    # Check 3' end (end of read)
    end_start = max(read_length - 120, boundary)
    three_prime_density = motif_array[end_start:].mean() if end_start < read_length else 0.0

    has_5prime = five_prime_density >= min_terminal_fraction
    has_3prime = three_prime_density >= min_terminal_fraction

    if has_5prime and has_3prime:
        return "both"
    elif has_5prime:
        return "5prime"
    elif has_3prime:
        return "3prime"
    else:
        return "internal"
