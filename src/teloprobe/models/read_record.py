"""Per-read input data model."""

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from ..constants import Platform, Mode


@dataclass
class ReadRecord:
    """A single sequencing read with metadata."""

    read_id: str
    sequence: str
    qualities: Optional[np.ndarray]
    read_length: int
    platform: Platform = Platform.UNKNOWN
    mode: Mode = Mode.TELOSEQ
    sample_id: str = ""
    barcode: Optional[str] = None
    mean_quality: float = 0.0

    @classmethod
    def from_pysam(cls, record, platform: Platform = Platform.UNKNOWN,
                   mode: Mode = Mode.TELOSEQ, sample_id: str = ""):
        """Create ReadRecord from a pysam AlignedSegment or FastxRecord."""
        if hasattr(record, "query_sequence"):
            # pysam AlignedSegment (BAM)
            seq = record.query_sequence or ""
            quals = record.query_qualities
            name = record.query_name
        else:
            # pysam FastxRecord (FASTQ)
            seq = record.sequence
            if record.quality:
                qa = getattr(record, "get_quality_array", None)
                if qa and callable(qa):
                    quals = np.array(qa())
                elif hasattr(record, "quality_array"):
                    quals = np.array(record.quality_array)
                else:
                    # Parse quality string manually
                    quals = np.array([ord(c) - 33 for c in record.quality], dtype=np.int32)
            else:
                quals = None
            name = record.name

        mean_q = float(np.mean(quals)) if quals is not None and len(quals) > 0 else 0.0

        return cls(
            read_id=name,
            sequence=seq,
            qualities=quals,
            read_length=len(seq),
            platform=platform,
            mode=mode,
            sample_id=sample_id,
            mean_quality=mean_q,
        )
