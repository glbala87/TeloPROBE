"""Per-read telomere length estimation.

Combines motif scanning, boundary detection, and anchoring into
a unified per-read measurement with confidence scoring.
"""

import logging
from typing import Optional

import numpy as np
import pandas as pd

from ..config import Config
from ..constants import QCStatus
from ..models.read_record import ReadRecord
from ..models.telomere_call import TelomereCall, EvidenceLayers
from ..candidate.motif_scanner import (
    scan_motifs,
    compute_motif_density_fast,
    quick_filter,
    telomeric_fraction,
)
from ..candidate.tvr_decomposer import decompose_tvr, tvr_fraction
from ..candidate.read_filter import (
    trim_adapters,
    apply_preflight_filters,
    check_error_motifs,
    check_post_boundary_ccc,
)
from ..segment.hmm_labeler import build_hmm, label_windows, hmm_boundary
from ..segment.changepoint import detect_changepoints, select_boundary_changepoint
from ..segment.boundary import (
    consensus_boundary,
    classify_boundary,
    determine_telomere_side,
)
from .confidence import compute_confidence_score

logger = logging.getLogger(__name__)


def process_single_read(
    record: ReadRecord,
    config: Config,
    hmm=None,
) -> TelomereCall:
    """Process a single read through the full telomere detection pipeline.

    Steps:
    1. Adapter trimming (Mode B only)
    2. Preflight QC filters
    3. Motif scanning (Pass 1: fast k-mer scan)
    4. Quick rejection filter
    5. Motif density computation
    6. HMM labeling (Pass 2: probabilistic)
    7. Change-point detection
    8. Consensus boundary
    9. Boundary QC classification
    10. TVR decomposition
    11. Confidence scoring
    """
    call = TelomereCall(
        read_id=record.read_id,
        read_length=record.read_length,
        platform=record.platform.value,
        sample_id=record.sample_id,
    )

    # Step 1: Adapter trimming (TeloPROBE mode)
    if config.mode.value == "teloprobe":
        record, barcode = trim_adapters(record)
        call.read_length = record.read_length

    # Step 2: Preflight filters
    status = apply_preflight_filters(record, config)
    if status != QCStatus.GOOD:
        call.qc_status = status
        return call

    # Step 3: Quick motif filter
    if not quick_filter(record.sequence, min_repeats=config.min_repeats):
        call.qc_status = QCStatus.TOO_FEW_REPEATS
        return call

    # Step 4: Full motif scan
    motif_array = scan_motifs(record.sequence, include_variants=True)
    call.motif_count = int(motif_array.sum() // 6)  # approximate motif count
    call.motif_density = telomeric_fraction(motif_array)

    # Step 5: Compute motif density vector
    density = compute_motif_density_fast(motif_array, window_size=config.motif_window_size)

    # Step 6: HMM labeling
    hmm_bnd = None
    state_path = None
    if hmm is not None:
        state_path = label_windows(density, hmm, window_size=config.hmm_window_size)
        hmm_bnd = hmm_boundary(state_path, window_size=config.hmm_window_size)
        call.hmm_boundary = hmm_bnd

    # Step 7: Change-point detection
    changepoints = detect_changepoints(
        density,
        method=config.changepoint_method,
        penalty=config.changepoint_penalty,
        min_size=config.changepoint_min_size,
    )
    pelt_bnd = select_boundary_changepoint(changepoints, density)
    call.pelt_boundary = pelt_bnd
    call.changepoint_positions = [cp for cp in changepoints if cp < len(density)]

    # Step 8: Consensus boundary
    boundary, method, agreement = consensus_boundary(
        hmm_bnd, pelt_bnd, motif_array
    )
    call.boundary_pos = boundary
    call.boundary_method = method

    # Step 9: QC classification
    qc = classify_boundary(
        boundary, record.read_length, motif_array,
        record.qualities, call.motif_count, config,
    )

    # Additional QC checks that need the full sequence
    if qc == QCStatus.GOOD and boundary is not None:
        if check_error_motifs(
            record.sequence, boundary,
            config.max_errors, config.error_distance
        ):
            qc = QCStatus.TOO_ERRORFUL

        if check_post_boundary_ccc(
            record.sequence, boundary,
            config.post_boundary_ccc_threshold
        ):
            qc = QCStatus.TELOMERE_ONLY

    call.qc_status = qc

    # Step 10: Telomere length and side
    if boundary is not None and qc == QCStatus.GOOD:
        call.telomere_length = boundary
        call.telomere_side = determine_telomere_side(
            motif_array, boundary, record.read_length
        )

    # Step 11: TVR decomposition
    if boundary is not None and boundary > 0:
        tvr_comp = decompose_tvr(record.sequence, boundary)
        call.tvr_composition = tvr_comp
        call.tvr_fraction = tvr_fraction(tvr_comp)

    # Step 12: Evidence layers and confidence
    call.evidence = EvidenceLayers(
        motif_score=min(1.0, call.motif_density / 0.8),
        boundary_score=agreement,
        anchor_score=0.0,  # updated after alignment
    )
    call.confidence_score = compute_confidence_score(call, config)

    return call


def compute_read_level(
    records: list[ReadRecord],
    config: Config,
) -> pd.DataFrame:
    """Process all reads and return read-level DataFrame.

    Args:
        records: List of ReadRecord objects from ingestion.
        config: Pipeline configuration.

    Returns:
        DataFrame with one row per read, columns matching read_level.tsv schema.
    """
    # Build HMM once
    hmm = build_hmm(
        n_states=config.hmm_n_states,
        random_seed=config.random_seed,
    )

    rows = []
    total = len(records)

    for i, record in enumerate(records):
        if (i + 1) % 1000 == 0:
            logger.info("Processing read %d / %d", i + 1, total)

        call = process_single_read(record, config, hmm=hmm)
        rows.append(call.to_row())

    df = pd.DataFrame(rows)
    logger.info(
        "Processed %d reads: %d Good, %d filtered",
        len(df),
        (df["qc_status"] == "Good").sum(),
        (df["qc_status"] != "Good").sum(),
    )
    return df
