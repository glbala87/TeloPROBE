"""Per-read telomere measurement data model."""

from dataclasses import dataclass, field
from typing import Optional

from ..constants import QCStatus


@dataclass
class EvidenceLayers:
    """Scores from three evidence layers."""

    motif_score: float = 0.0      # 0-1: motif density support
    boundary_score: float = 0.0   # 0-1: boundary detection clarity
    anchor_score: float = 0.0     # 0-1: subtelomeric anchoring quality

    @property
    def as_dict(self) -> dict:
        return {
            "motif_score": self.motif_score,
            "boundary_score": self.boundary_score,
            "anchor_score": self.anchor_score,
        }


@dataclass
class TelomereCall:
    """Complete telomere measurement for a single read."""

    read_id: str
    read_length: int
    platform: str = ""
    sample_id: str = ""

    # Boundary detection
    boundary_pos: Optional[int] = None
    telomere_length: Optional[int] = None
    telomere_side: str = "5prime"  # 5prime, 3prime, both, internal
    boundary_method: str = ""  # hmm, pelt, consensus

    # Quality
    qc_status: QCStatus = QCStatus.GOOD
    confidence_score: float = 0.0
    evidence: EvidenceLayers = field(default_factory=EvidenceLayers)

    # Motif analysis
    motif_density: float = 0.0
    motif_count: int = 0
    tvr_composition: dict = field(default_factory=dict)
    tvr_fraction: float = 0.0

    # HMM/changepoint details
    hmm_boundary: Optional[int] = None
    pelt_boundary: Optional[int] = None
    changepoint_positions: list = field(default_factory=list)

    # Anchoring
    chr_arm: Optional[str] = None
    haplotype: Optional[str] = None
    allele_cluster: Optional[int] = None
    alignment_identity: Optional[float] = None
    mapping_quality: Optional[int] = None
    anchor_confidence: str = "none"  # high, medium, low, none

    def to_row(self) -> dict:
        """Convert to flat dict for DataFrame output."""
        return {
            "read_id": self.read_id,
            "platform": self.platform,
            "read_length": self.read_length,
            "telomere_length_bp": self.telomere_length if self.telomere_length else -1,
            "telomere_side": self.telomere_side,
            "motif_score": round(self.evidence.motif_score, 4),
            "motif_density": round(self.motif_density, 4),
            "boundary_confidence": round(self.evidence.boundary_score, 4),
            "boundary_method": self.boundary_method,
            "subtelomere_anchor": round(self.evidence.anchor_score, 4),
            "chromosome_arm": self.chr_arm or "",
            "haplotype": self.haplotype or "",
            "allele_cluster": self.allele_cluster if self.allele_cluster is not None else "",
            "tvr_fraction": round(self.tvr_fraction, 4),
            "confidence_score": round(self.confidence_score, 4),
            "qc_status": self.qc_status.value,
            "alignment_identity": round(self.alignment_identity, 4) if self.alignment_identity else "",
            "mapping_quality": self.mapping_quality if self.mapping_quality is not None else "",
            "anchor_confidence": self.anchor_confidence,
        }
