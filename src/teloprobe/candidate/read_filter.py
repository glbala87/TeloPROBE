"""Read filtering: adapter trimming, quality checks, error motif detection.

Ported and extended from the original wf-teloprobe process_reads.py.
"""

import re
import logging
from typing import Optional

import numpy as np

try:
    import edlib
except ImportError:
    edlib = None

from ..constants import (
    QCStatus,
    TELOSEQ_BARCODES,
    TELOMERE_MARKER,
    ERROR_MOTIFS,
    MOTIF_LEN,
)
from ..models.read_record import ReadRecord
from ..config import Config

logger = logging.getLogger(__name__)


def trim_adapters(
    record: ReadRecord,
    barcodes: Optional[list[str]] = None,
    prefix_length: int = 200,
    max_mismatches: int = 3,
) -> tuple[ReadRecord, Optional[str]]:
    """Trim TeloPROBE barcode adapters from read start.

    Searches the first prefix_length bases for barcode sequences
    using edlib edit-distance alignment.

    Returns:
        Tuple of (trimmed ReadRecord, detected barcode or None).
    """
    if edlib is None:
        logger.warning("edlib not installed, skipping adapter trimming")
        return record, None

    if barcodes is None:
        barcodes = TELOSEQ_BARCODES

    seq = record.sequence
    prefix = seq[:prefix_length]
    detected_barcode = None
    trim_pos = 0

    # Try each barcode
    for bc in barcodes:
        # Build query: barcode + telomere marker
        query = bc + TELOMERE_MARKER
        result = edlib.align(
            query, prefix, mode="HW", task="locations",
            k=max_mismatches,
        )
        if result["editDistance"] >= 0 and result["locations"]:
            end_pos = result["locations"][0][1] + 1
            if end_pos > trim_pos:
                trim_pos = end_pos
                detected_barcode = bc

    # Fallback: search for telomere marker alone
    if trim_pos == 0:
        result = edlib.align(
            TELOMERE_MARKER, prefix, mode="HW", task="locations", k=1,
        )
        if result["editDistance"] >= 0 and result["locations"]:
            trim_pos = result["locations"][0][1] + 1

    if trim_pos > 0:
        trimmed_seq = seq[trim_pos:]
        trimmed_quals = (
            record.qualities[trim_pos:] if record.qualities is not None else None
        )
        return ReadRecord(
            read_id=record.read_id,
            sequence=trimmed_seq,
            qualities=trimmed_quals,
            read_length=len(trimmed_seq),
            platform=record.platform,
            mode=record.mode,
            sample_id=record.sample_id,
            barcode=detected_barcode,
            mean_quality=record.mean_quality,
        ), detected_barcode

    return record, None


def apply_preflight_filters(
    record: ReadRecord,
    config: Config,
) -> QCStatus:
    """Apply basic read-level quality filters before boundary detection.

    Returns QCStatus.GOOD if the read passes all filters.
    """
    if record.read_length < config.min_read_length_boundary:
        return QCStatus.TOO_SHORT

    if config.min_quality > 0 and record.mean_quality < config.min_quality:
        return QCStatus.LOW_SUBTELO_QUAL

    return QCStatus.GOOD


def check_error_motifs(
    sequence: str,
    boundary: int,
    max_errors: int = 5,
    error_distance: int = 500,
) -> bool:
    """Check for clustering of basecalling error motifs near the boundary.

    Returns True if too many error motifs are clustered (read should be
    flagged as TooErrorful).
    """
    positions = []
    for motif in ERROR_MOTIFS:
        for match in re.finditer(re.escape(motif), sequence):
            positions.append(match.start())

    if len(positions) < max_errors:
        return False

    return _largest_cluster(positions, error_distance) >= max_errors


def _largest_cluster(positions: list[int], max_distance: int) -> int:
    """Find the largest cluster of positions within max_distance of each other."""
    if not positions:
        return 0

    positions = sorted(positions)
    max_cluster = 1
    current_cluster = 1

    for i in range(1, len(positions)):
        if positions[i] - positions[i - 1] <= max_distance:
            current_cluster += 1
            max_cluster = max(max_cluster, current_cluster)
        else:
            current_cluster = 1

    return max_cluster


def check_start_repeats(
    motif_array: np.ndarray,
    boundary: int,
    start_window_frac: float = 0.3,
    start_repeats_frac: float = 0.8,
) -> bool:
    """Check that the first portion of the read before boundary is mostly telomeric.

    Returns True if the start fails the repeat check (should be flagged).
    """
    window_end = int(boundary * start_window_frac)
    if window_end <= 0:
        return True  # boundary too close to start

    window = motif_array[:window_end]
    density = window.sum() / len(window) if len(window) > 0 else 0.0
    return density < start_repeats_frac


def check_post_boundary_quality(
    qualities: Optional[np.ndarray],
    boundary: int,
    min_qual: float = 9.0,
) -> bool:
    """Check that bases after the boundary have sufficient quality.

    Returns True if quality is too low (should be flagged).
    """
    if qualities is None:
        return False  # no quality info, pass

    post_quals = qualities[boundary:]
    if len(post_quals) == 0:
        return True  # nothing after boundary

    return float(np.median(post_quals)) < min_qual


def check_post_boundary_ccc(
    sequence: str,
    boundary: int,
    threshold: float = 0.25,
) -> bool:
    """Check if post-boundary region has excessive CCC content.

    High CCC after boundary suggests the read is entirely telomeric
    with no subtelomeric flank. Returns True if should be flagged.
    """
    post_seq = sequence[boundary:]
    if len(post_seq) == 0:
        return True

    ccc_count = post_seq.count("CCC") + post_seq.count("ccc")
    ccc_fraction = (ccc_count * 3) / len(post_seq)
    return ccc_fraction > threshold
