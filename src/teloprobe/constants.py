"""Telomere motifs, platform presets, and enumeration types."""

from enum import Enum


class Platform(str, Enum):
    ONT = "ont"
    PACBIO_HIFI = "pacbio"
    UNKNOWN = "unknown"


class Mode(str, Enum):
    WGS = "wgs"
    TELOSEQ = "teloprobe"


class QCStatus(str, Enum):
    GOOD = "Good"
    TOO_SHORT = "TooShort"
    TOO_FEW_REPEATS = "TooFewRepeats"
    TOO_CLOSE_START = "TooCloseStart"
    TOO_CLOSE_END = "TooCloseEnd"
    START_NOT_REPEATS = "StartNotRepeats"
    LOW_SUBTELO_QUAL = "LowSubTeloQual"
    TELOMERE_ONLY = "TelomereOnly"
    TOO_ERRORFUL = "TooErrorful"
    BAD_ALIGN = "BadAlign"
    LOW_CONFIDENCE = "LowConfidence"
    NO_BOUNDARY = "NoBoundary"


# Motif length (hexanucleotide repeats)
MOTIF_LEN = 6

# Canonical telomere motif (reverse complement of TTAGGG, as found in reads)
CANONICAL_MOTIF = "TAACCC"

# All rotations of the canonical motif for forward strand
FORWARD_MOTIFS = frozenset({
    "TTAGGG", "TAGGGT", "AGGGTT", "GGGTTA", "GGTTAG", "GTTAGG",
})

# All rotations for reverse complement strand
REVERSE_MOTIFS = frozenset({
    "CCCTAA", "CCTAAC", "CTAACC", "TAACCC", "AACCCT", "ACCCTA",
})

# Telomere variant repeat (TVR) motifs — biologically real variants
VARIANT_MOTIFS = frozenset({
    "CACCCT", "ACCCCT", "CCCAAA", "CCCCGA",
    "CCCTGA", "CCCTCA", "CCCTAC", "CCCTAT", "CCCTAG",
    # Additional known TVRs
    "TTGGGG", "TGAGGG", "TCAGGG", "TTCGGG", "ATAGGG",
    "GTAGGG", "CTAGGG",
})

# Combined set for scanning
ALL_TELOMERE_MOTIFS = FORWARD_MOTIFS | REVERSE_MOTIFS | VARIANT_MOTIFS

# Basecalling error motifs (known ONT artifacts)
ERROR_MOTIFS = frozenset({
    "GTATAG", "CGCGCGCG", "CCACCG", "AGCGACAG",
    "ATAAGT", "CCTCGTCC", "TATAGT", "AGTACT",
    "GAGTCC", "TATACA", "TGGTCC", "CTCTCCTCT",
})

# TeloPROBE barcodes
TELOSEQ_BARCODES = [
    "CACAAAGACGATAGCACTTCGCAGTCTCTATTTG",
    "ACAGACGACTACAAACGGAATCGAGATACAGAGC",
    "CCTGGTAACTGGGACACAAGACTCTATATAGAGC",
    "TAGGGAAACACGATAGAATCCGAATAAGGATACC",
    "AAGGTTACACAAACCCTGGACAAGAGAAGACCAG",
    "ACTACTGCTGCCTAAGACCCTAACATTTCCCTCG",
    "ATGCTTCTTTCGATGCTCAGTACGTCAGAATACC",
    "AACAGGGGATTTCGATCTATAGATTGAAGAGAGG",
    "TACAGTCCGAGCCTCATGTGATCTATAGCTACTG",
    "AGAAGGCATCGAGCGGAGTACTATATCAGCTCCG",
    "GTTCATAACTCGTAATGGATTGCACTAAGCTGCG",
    "AGCATATGATCGAGGCTTCTAGAACTCAAATCGC",
]

# Telomere marker in barcode adapter
TELOMERE_MARKER = "CCTAACC"

# Platform-specific presets
PLATFORM_PRESETS = {
    Platform.ONT: {
        "min_quality": 9,
        "min_length": 100,
        "motif_density_threshold": 0.6,
        "motif_mismatch_tolerance": 0.25,
        "error_motif_check": True,
        "filter_width": 10,
        "min_repeats": 100,
        "start_window_frac": 0.3,
        "start_repeats_frac": 0.8,
        "min_qual_non_telo": 9,
        "post_boundary_ccc_threshold": 0.25,
        "max_errors": 5,
        "error_distance": 500,
        "minimap2_preset": "lr:hq",
        "identity_threshold": 0.8,
        "mapq_threshold": 20,
        "changepoint_penalty": None,  # auto
        "hmm_density_telo": 0.85,
        "hmm_density_subtelo": 0.05,
    },
    Platform.PACBIO_HIFI: {
        "min_quality": 20,
        "min_length": 100,
        "motif_density_threshold": 0.75,
        "motif_mismatch_tolerance": 0.15,
        "error_motif_check": False,
        "filter_width": 10,
        "min_repeats": 100,
        "start_window_frac": 0.3,
        "start_repeats_frac": 0.85,
        "min_qual_non_telo": 20,
        "post_boundary_ccc_threshold": 0.20,
        "max_errors": 3,
        "error_distance": 500,
        "minimap2_preset": "map-hifi",
        "identity_threshold": 0.9,
        "mapq_threshold": 30,
        "changepoint_penalty": None,
        "hmm_density_telo": 0.92,
        "hmm_density_subtelo": 0.02,
    },
}

# Evidence weights for confidence scoring by mode
EVIDENCE_WEIGHTS = {
    Mode.TELOSEQ: {"motif": 0.4, "boundary": 0.4, "anchor": 0.2},
    Mode.WGS: {"motif": 0.2, "boundary": 0.3, "anchor": 0.5},
}

# Default chromosome arms (T2T-CHM13 compatible)
CHROMOSOME_ARMS = [
    f"chr{c}{arm}"
    for c in list(range(1, 23)) + ["X", "Y"]
    for arm in ["p", "q"]
]
