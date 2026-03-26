"""Chromosome-arm assignment from subtelomeric alignment.

Parses reference names, applies quality filters, and assigns
reads to chromosome arms with confidence levels.
"""

import logging
import re
from typing import Optional

import pysam
import pandas as pd

logger = logging.getLogger(__name__)

# Patterns for parsing reference sequence names
# Supports: chr1_PATERNAL_P, chr1_mat_p, chr1p, etc.
_ARM_PATTERNS = [
    # T2T-CHM13 style: chr1_PATERNAL_P or chr1_pat_p
    re.compile(
        r"^(chr\d+|chrX|chrY)[-_](pat(?:ernal)?|mat(?:ernal)?|hap[12])[-_]([pq])$",
        re.IGNORECASE,
    ),
    # Simple style: chr1p or chr1q
    re.compile(
        r"^(chr\d+|chrX|chrY)([pq])$",
        re.IGNORECASE,
    ),
    # Underscore style: chr1_p or chr1_q
    re.compile(
        r"^(chr\d+|chrX|chrY)[-_]([pq])$",
        re.IGNORECASE,
    ),
]


def parse_reference_name(ref_name: str) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """Parse a reference sequence name into (chromosome, haplotype, arm).

    Returns:
        Tuple of (chromosome, haplotype, arm_end).
        E.g., ("chr1", "pat", "p") or ("chrX", None, "q").
        Returns (None, None, None) if unparseable.
    """
    for pattern in _ARM_PATTERNS:
        m = pattern.match(ref_name)
        if m:
            groups = m.groups()
            if len(groups) == 3:
                chrom = groups[0].lower()
                haplotype = _normalize_haplotype(groups[1])
                arm = groups[2].lower()
                return chrom, haplotype, arm
            elif len(groups) == 2:
                chrom = groups[0].lower()
                arm = groups[1].lower()
                return chrom, None, arm

    return None, None, None


def _normalize_haplotype(hap_str: str) -> str:
    """Normalize haplotype string to 'pat', 'mat', or 'unknown'."""
    hap = hap_str.lower()
    if hap in ("pat", "paternal", "hap1"):
        return "pat"
    elif hap in ("mat", "maternal", "hap2"):
        return "mat"
    return "unknown"


def assign_arm(
    record: pysam.AlignedSegment,
    identity_threshold: float = 0.8,
    mapq_threshold: int = 20,
) -> dict:
    """Assign a chromosome arm to an aligned read.

    Args:
        record: pysam AlignedSegment from subtelomeric alignment.
        identity_threshold: Minimum alignment identity.
        mapq_threshold: Minimum mapping quality.

    Returns:
        Dict with keys: read_id, chr_arm, haplotype, identity,
        mapq, anchor_confidence.
    """
    result = {
        "read_id": record.query_name,
        "chr_arm": None,
        "haplotype": None,
        "alignment_identity": None,
        "mapping_quality": record.mapping_quality,
        "anchor_confidence": "none",
    }

    if record.is_unmapped:
        return result

    # Compute alignment identity
    identity = _compute_identity(record)
    result["alignment_identity"] = identity
    result["mapping_quality"] = record.mapping_quality

    # Apply quality filters
    if identity < identity_threshold or record.mapping_quality < mapq_threshold:
        result["anchor_confidence"] = "low"
        return result

    # Parse reference name
    ref_name = record.reference_name
    chrom, haplotype, arm = parse_reference_name(ref_name)

    if chrom and arm:
        result["chr_arm"] = f"{chrom}{arm}"
        result["haplotype"] = haplotype or "unknown"

        # Confidence based on mapping quality and identity
        if record.mapping_quality >= 60 and identity >= 0.95:
            result["anchor_confidence"] = "high"
        elif record.mapping_quality >= mapq_threshold and identity >= identity_threshold:
            result["anchor_confidence"] = "medium"
        else:
            result["anchor_confidence"] = "low"
    else:
        result["anchor_confidence"] = "low"

    return result


def _compute_identity(record: pysam.AlignedSegment) -> float:
    """Compute gap-compressed alignment identity.

    Uses the de:f tag if available (minimap2), otherwise computes
    from CIGAR and NM tag.
    """
    # Try de:f tag (gap-compressed divergence from minimap2)
    try:
        de = record.get_tag("de")
        return 1.0 - de
    except KeyError:
        pass

    # Compute from NM tag and alignment length
    try:
        nm = record.get_tag("NM")
        aligned_length = record.query_alignment_length
        if aligned_length > 0:
            return 1.0 - (nm / aligned_length)
    except KeyError:
        pass

    # Fallback: compute from CIGAR
    if record.cigartuples:
        matches = sum(l for op, l in record.cigartuples if op in (0, 7))
        total = sum(l for op, l in record.cigartuples if op in (0, 1, 2, 7, 8))
        if total > 0:
            return matches / total

    return 0.0


def process_alignments(
    bam_path: str,
    identity_threshold: float = 0.8,
    mapq_threshold: int = 20,
) -> pd.DataFrame:
    """Process all alignments and return arm assignment DataFrame."""
    assignments = []

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for record in bam:
            assignment = assign_arm(record, identity_threshold, mapq_threshold)
            assignments.append(assignment)

    df = pd.DataFrame(assignments)
    logger.info(
        "Processed %d alignments: %d with arm assignment",
        len(df),
        df["chr_arm"].notna().sum(),
    )
    return df
