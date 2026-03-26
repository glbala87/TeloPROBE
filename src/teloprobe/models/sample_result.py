"""Sample-level result data model."""

from dataclasses import dataclass, field
from typing import Optional

from ..constants import Platform, Mode


@dataclass
class SampleResult:
    """Aggregated telomere length estimate for a sample."""

    sample_id: str
    platform: Platform = Platform.UNKNOWN
    mode: Mode = Mode.TELOSEQ
    total_reads: int = 0
    telomeric_reads: int = 0
    passing_reads: int = 0
    informative_reads: int = 0
    median_tl: float = 0.0
    mean_tl: float = 0.0
    trimmed_mean_tl: float = 0.0
    ci_lower: float = 0.0
    ci_upper: float = 0.0
    cv: float = 0.0
    arms_assigned: int = 0
    qc_summary: dict = field(default_factory=dict)

    def to_row(self) -> dict:
        """Convert to flat dict for DataFrame output."""
        return {
            "sample": self.sample_id,
            "platform": self.platform.value,
            "mode": self.mode.value,
            "total_reads": self.total_reads,
            "telomeric_reads": self.telomeric_reads,
            "informative_reads": self.informative_reads,
            "median_tl": round(self.median_tl, 1),
            "trimmed_mean_tl": round(self.trimmed_mean_tl, 1),
            "global_ci_lower": round(self.ci_lower, 1),
            "global_ci_upper": round(self.ci_upper, 1),
            "cv": round(self.cv, 4),
            "arms_assigned": self.arms_assigned,
            "qc_summary": str(self.qc_summary),
        }
