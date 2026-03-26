"""Fast telomeric motif scanning using compiled regex and vectorized operations.

Two-pass candidate detection:
  Pass 1: k-mer scan for telomere-like motifs (this module)
  Pass 2: HMM/probabilistic labeler (segment module)
"""

import re
import logging
from typing import Optional

import numpy as np

from ..constants import (
    MOTIF_LEN,
    FORWARD_MOTIFS,
    REVERSE_MOTIFS,
    VARIANT_MOTIFS,
    ALL_TELOMERE_MOTIFS,
)

logger = logging.getLogger(__name__)

# Precompile regex patterns for all telomere motifs
_CANONICAL_PATTERN = re.compile(
    "|".join(re.escape(m) for m in sorted(FORWARD_MOTIFS | REVERSE_MOTIFS))
)
_ALL_PATTERN = re.compile(
    "|".join(re.escape(m) for m in sorted(ALL_TELOMERE_MOTIFS))
)
_VARIANT_PATTERN = re.compile(
    "|".join(re.escape(m) for m in sorted(VARIANT_MOTIFS))
)


def scan_motifs(
    sequence: str,
    include_variants: bool = True,
) -> np.ndarray:
    """Scan sequence for telomeric motifs, returning a binary hit array.

    For each position in the sequence, marks 1 if a motif starts at that
    position (or overlaps the 6bp window), 0 otherwise.

    Args:
        sequence: DNA sequence string.
        include_variants: Whether to include TVR motifs in scanning.

    Returns:
        Binary numpy array of length len(sequence), where 1 indicates
        a telomeric motif covers that position.
    """
    seq_len = len(sequence)
    if seq_len < MOTIF_LEN:
        return np.zeros(seq_len, dtype=np.int8)

    hits = np.zeros(seq_len, dtype=np.int8)
    pattern = _ALL_PATTERN if include_variants else _CANONICAL_PATTERN

    for match in pattern.finditer(sequence):
        start = match.start()
        end = match.end()
        hits[start:end] = 1

    return hits


def scan_motifs_by_type(
    sequence: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Scan for canonical and variant motifs separately.

    Returns:
        Tuple of (canonical_hits, variant_hits) binary arrays.
    """
    seq_len = len(sequence)
    canonical = np.zeros(seq_len, dtype=np.int8)
    variant = np.zeros(seq_len, dtype=np.int8)

    for match in _CANONICAL_PATTERN.finditer(sequence):
        canonical[match.start():match.end()] = 1

    for match in _VARIANT_PATTERN.finditer(sequence):
        variant[match.start():match.end()] = 1

    return canonical, variant


def compute_motif_density(
    motif_array: np.ndarray,
    window_size: int = 120,
) -> np.ndarray:
    """Compute sliding-window motif density across the read.

    Args:
        motif_array: Binary array from scan_motifs.
        window_size: Window size in bp for density calculation.

    Returns:
        Array of motif density values (0.0 to 1.0), same length as input.
    """
    seq_len = len(motif_array)
    if seq_len < window_size:
        total = motif_array.sum()
        return np.full(seq_len, total / seq_len if seq_len > 0 else 0.0)

    # Use cumulative sum for efficient windowed computation
    cumsum = np.cumsum(motif_array, dtype=np.float64)
    density = np.zeros(seq_len, dtype=np.float64)

    half_w = window_size // 2

    for i in range(seq_len):
        left = max(0, i - half_w)
        right = min(seq_len, i + half_w)
        window_sum = cumsum[right - 1] - (cumsum[left - 1] if left > 0 else 0)
        density[i] = window_sum / (right - left)

    return density


def compute_motif_density_fast(
    motif_array: np.ndarray,
    window_size: int = 120,
) -> np.ndarray:
    """Vectorized motif density computation using uniform filter.

    Faster than compute_motif_density for large arrays.
    """
    from scipy.ndimage import uniform_filter1d
    return uniform_filter1d(
        motif_array.astype(np.float64), size=window_size, mode="nearest"
    )


def quick_filter(
    sequence: str,
    min_repeats: int = 100,
    include_variants: bool = True,
) -> bool:
    """Fast rejection filter: check if sequence has enough telomeric repeats.

    This is a quick first pass to avoid expensive HMM/change-point
    computation on clearly non-telomeric reads.
    """
    pattern = _ALL_PATTERN if include_variants else _CANONICAL_PATTERN
    count = len(pattern.findall(sequence))
    return count >= min_repeats


def longest_telomeric_run(motif_array: np.ndarray) -> int:
    """Find the longest consecutive run of telomeric positions."""
    if len(motif_array) == 0:
        return 0

    max_run = 0
    current_run = 0
    for val in motif_array:
        if val:
            current_run += 1
            max_run = max(max_run, current_run)
        else:
            current_run = 0
    return max_run


def telomeric_fraction(motif_array: np.ndarray) -> float:
    """Fraction of the read covered by telomeric motifs."""
    if len(motif_array) == 0:
        return 0.0
    return float(motif_array.sum()) / len(motif_array)


def entropy_score(sequence: str, window_size: int = 120) -> float:
    """Compute Shannon entropy of the sequence as a complexity measure.

    Low entropy = low complexity (repetitive), high entropy = complex.
    Telomeric sequence should have low entropy.
    """
    from collections import Counter
    if len(sequence) < window_size:
        window_size = len(sequence)
    if window_size == 0:
        return 0.0

    # Use first window_size bases
    counts = Counter(sequence[:window_size])
    total = sum(counts.values())
    entropy = 0.0
    for count in counts.values():
        if count > 0:
            p = count / total
            entropy -= p * np.log2(p)
    return entropy
