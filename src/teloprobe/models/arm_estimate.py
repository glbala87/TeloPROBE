"""Per-chromosome-arm telomere estimate data model."""

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class ArmEstimate:
    """Telomere length estimate for a chromosome arm."""

    sample_id: str
    chromosome_arm: str
    haplotype: str = "unknown"
    allele_id: Optional[int] = None
    n_reads: int = 0
    median_tl: float = 0.0
    mean_tl: float = 0.0
    trimmed_mean_tl: float = 0.0
    std_tl: float = 0.0
    cv: float = 0.0
    q1: float = 0.0
    q3: float = 0.0
    ci_lower: float = 0.0
    ci_upper: float = 0.0
    anchor_confidence: str = "none"
    tvr_consensus: dict = field(default_factory=dict)

    def to_row(self) -> dict:
        """Convert to flat dict for DataFrame output."""
        return {
            "sample": self.sample_id,
            "chromosome_arm": self.chromosome_arm,
            "haplotype": self.haplotype,
            "allele_id": self.allele_id if self.allele_id is not None else "",
            "n_reads": self.n_reads,
            "median_tl": round(self.median_tl, 1),
            "mean_tl": round(self.mean_tl, 1),
            "trimmed_mean_tl": round(self.trimmed_mean_tl, 1),
            "ci_lower": round(self.ci_lower, 1),
            "ci_upper": round(self.ci_upper, 1),
            "cv": round(self.cv, 4),
            "anchor_confidence": self.anchor_confidence,
        }
