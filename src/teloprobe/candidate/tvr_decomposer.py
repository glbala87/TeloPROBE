"""Telomere Variant Repeat (TVR) decomposition and annotation.

TVRs are biologically real variants of the canonical TTAGGG repeat.
Their composition varies across chromosome arms and alleles, providing
signatures for arm/allele identification.
"""

import re
from collections import Counter
from typing import Optional

from ..constants import FORWARD_MOTIFS, REVERSE_MOTIFS, VARIANT_MOTIFS, MOTIF_LEN


# Map each motif to its canonical name for reporting
_MOTIF_NAMES = {}
for m in FORWARD_MOTIFS | REVERSE_MOTIFS:
    _MOTIF_NAMES[m] = "TTAGGG"
for m in VARIANT_MOTIFS:
    _MOTIF_NAMES[m] = m  # TVRs keep their own identity

# All motifs as a single regex (longest first for proper matching)
_ALL_SORTED = sorted(
    FORWARD_MOTIFS | REVERSE_MOTIFS | VARIANT_MOTIFS,
    key=len, reverse=True,
)
_TVR_PATTERN = re.compile("|".join(re.escape(m) for m in _ALL_SORTED))


def decompose_tvr(
    sequence: str,
    boundary: Optional[int] = None,
) -> dict[str, int]:
    """Decompose the telomeric region into motif counts.

    Args:
        sequence: Full read sequence.
        boundary: Position of telomere-subtelomere boundary.
                  If provided, only the region [0:boundary] is analyzed.
                  If None, the entire sequence is analyzed.

    Returns:
        Dict mapping motif name -> count.
        Canonical TTAGGG and its rotations are merged under "TTAGGG".
    """
    region = sequence[:boundary] if boundary else sequence
    counts = Counter()

    for match in _TVR_PATTERN.finditer(region):
        motif = match.group()
        name = _MOTIF_NAMES.get(motif, motif)
        counts[name] += 1

    return dict(counts)


def tvr_fraction(composition: dict[str, int]) -> float:
    """Compute the fraction of non-canonical (variant) repeats.

    Returns 0.0 if all repeats are canonical TTAGGG, higher values
    indicate more TVR diversity.
    """
    total = sum(composition.values())
    if total == 0:
        return 0.0
    canonical = composition.get("TTAGGG", 0)
    return (total - canonical) / total


def tvr_profile_vector(composition: dict[str, int]) -> list[float]:
    """Convert TVR composition to a fixed-length normalized vector.

    Used as feature input for allele clustering.
    """
    # Fixed motif order for consistent vectorization
    motif_order = sorted(VARIANT_MOTIFS)
    total = sum(composition.values())
    if total == 0:
        return [0.0] * len(motif_order)

    return [composition.get(m, 0) / total for m in motif_order]


def tvr_pattern_string(composition: dict[str, int]) -> str:
    """Human-readable TVR pattern summary.

    Example: "TTAGGG:450 TGAGGG:12 TCAGGG:3"
    """
    if not composition:
        return "none"
    sorted_items = sorted(composition.items(), key=lambda x: -x[1])
    return " ".join(f"{motif}:{count}" for motif, count in sorted_items)
