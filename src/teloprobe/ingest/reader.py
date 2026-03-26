"""BAM/CRAM/FASTQ reading via pysam."""

import logging
from pathlib import Path
from typing import Iterator

import pysam

from ..constants import Platform, Mode
from ..models.read_record import ReadRecord

logger = logging.getLogger(__name__)


def detect_format(path: Path) -> str:
    """Detect input file format from extension and magic bytes."""
    path = Path(path)
    suffixes = [s.lower() for s in path.suffixes]

    if ".bam" in suffixes:
        return "bam"
    elif ".cram" in suffixes:
        return "cram"
    elif any(s in suffixes for s in [".fastq", ".fq"]):
        return "fastq"
    elif ".sam" in suffixes:
        return "sam"

    # Try magic bytes
    try:
        with open(path, "rb") as f:
            header = f.read(4)
        if header[:3] == b"\x1f\x8b\x08":
            # Gzipped — could be FASTQ or BAM
            # BAM has specific magic after gzip header
            return "bam"  # pysam can handle both
        if header[:4] == b"@HD\t" or header[:4] == b"@SQ\t":
            return "sam"
    except Exception:
        pass

    return "fastq"  # default assumption


def read_input(
    path: str | Path,
    format: str = "auto",
    platform: Platform = Platform.UNKNOWN,
    mode: Mode = Mode.TELOSEQ,
    sample_id: str = "",
    min_length: int = 0,
    min_quality: float = 0.0,
) -> Iterator[ReadRecord]:
    """Read sequences from BAM/CRAM/FASTQ, yielding ReadRecord objects.

    Applies basic length and quality filtering during ingestion.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    if format == "auto":
        format = detect_format(path)

    logger.info("Reading %s as %s format", path, format)

    if format in ("bam", "cram", "sam"):
        yield from _read_xam(path, format, platform, mode, sample_id,
                             min_length, min_quality)
    elif format == "fastq":
        yield from _read_fastq(path, platform, mode, sample_id,
                               min_length, min_quality)
    else:
        raise ValueError(f"Unsupported format: {format}")


def _read_xam(
    path: Path,
    format: str,
    platform: Platform,
    mode: Mode,
    sample_id: str,
    min_length: int,
    min_quality: float,
) -> Iterator[ReadRecord]:
    """Read BAM/CRAM/SAM file."""
    open_mode = {"bam": "rb", "cram": "rc", "sam": "r"}[format]

    with pysam.AlignmentFile(str(path), open_mode, check_sq=False) as af:
        total = 0
        passed = 0
        for record in af:
            total += 1
            if record.is_secondary or record.is_supplementary:
                continue

            seq = record.query_sequence
            if not seq:
                continue
            if len(seq) < min_length:
                continue

            read = ReadRecord.from_pysam(
                record, platform=platform, mode=mode, sample_id=sample_id
            )

            if min_quality > 0 and read.mean_quality < min_quality:
                continue

            passed += 1
            yield read

    logger.info("Read %d records, %d passed filters from %s", total, passed, path)


def _read_fastq(
    path: Path,
    platform: Platform,
    mode: Mode,
    sample_id: str,
    min_length: int,
    min_quality: float,
) -> Iterator[ReadRecord]:
    """Read FASTQ file (plain or gzipped)."""
    total = 0
    passed = 0

    with pysam.FastxFile(str(path)) as fq:
        for record in fq:
            total += 1
            if len(record.sequence) < min_length:
                continue

            read = ReadRecord.from_pysam(
                record, platform=platform, mode=mode, sample_id=sample_id
            )

            if min_quality > 0 and read.mean_quality < min_quality:
                continue

            passed += 1
            yield read

    logger.info("Read %d records, %d passed filters from %s", total, passed, path)


def read_directory(
    dir_path: str | Path,
    platform: Platform = Platform.UNKNOWN,
    mode: Mode = Mode.TELOSEQ,
    sample_id: str = "",
    min_length: int = 0,
    min_quality: float = 0.0,
) -> Iterator[ReadRecord]:
    """Read all BAM/FASTQ files from a directory."""
    dir_path = Path(dir_path)
    extensions = {".bam", ".cram", ".fastq", ".fq", ".fastq.gz", ".fq.gz"}

    files = sorted(
        f for f in dir_path.iterdir()
        if f.is_file() and any(str(f).endswith(ext) for ext in extensions)
    )

    if not files:
        raise FileNotFoundError(f"No sequence files found in {dir_path}")

    for f in files:
        yield from read_input(
            f, platform=platform, mode=mode, sample_id=sample_id,
            min_length=min_length, min_quality=min_quality,
        )
