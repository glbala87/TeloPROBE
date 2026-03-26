"""Minimap2 wrapper for subtelomeric flank alignment.

Extracts the non-telomeric flank from each read and aligns it
to a subtelomeric reference for chromosome-arm assignment.
"""

import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

import pysam

from ..constants import Platform

logger = logging.getLogger(__name__)


def check_external_tools() -> None:
    """Verify that minimap2 and samtools are installed and on PATH.

    Raises:
        EnvironmentError: If either tool is not found.
    """
    missing = []
    for tool in ("minimap2", "samtools"):
        if shutil.which(tool) is None:
            missing.append(tool)
    if missing:
        raise EnvironmentError(
            f"Required external tool(s) not found on PATH: {', '.join(missing)}. "
            "Please install them before running the alignment step."
        )


def get_minimap2_preset(platform: Platform) -> str:
    """Return minimap2 preset for the given platform."""
    presets = {
        Platform.ONT: "lr:hq",
        Platform.PACBIO_HIFI: "map-hifi",
        Platform.UNKNOWN: "lr:hq",
    }
    return presets[platform]


def extract_flank_fastq(
    reads_bam: str | Path,
    boundaries_tsv: str | Path,
    output_fastq: str | Path,
    min_flank_length: int = 200,
) -> int:
    """Extract subtelomeric flanks from reads, write as FASTQ for alignment.

    For each read with a detected boundary, extracts the non-telomeric
    portion (after the boundary) and writes it as a FASTQ entry.

    Args:
        reads_bam: Path to BAM with candidate reads.
        boundaries_tsv: TSV with columns: read_id, boundary_pos, qc_status.
        output_fastq: Output FASTQ path.
        min_flank_length: Minimum flank length to include.

    Returns:
        Number of flanks written.
    """
    import pandas as pd

    boundaries = pd.read_csv(boundaries_tsv, sep="\t")
    boundary_map = dict(zip(boundaries["read_id"], boundaries["boundary_pos"]))
    good_reads = set(
        boundaries[boundaries["qc_status"] == "Good"]["read_id"]
    )

    count = 0
    with pysam.AlignmentFile(str(reads_bam), "rb", check_sq=False) as bam, \
         open(output_fastq, "w") as out:

        for record in bam:
            read_id = record.query_name
            if read_id not in good_reads:
                continue

            boundary = boundary_map.get(read_id)
            if boundary is None or not isinstance(boundary, (int, float)):
                continue
            boundary = int(boundary)

            seq = record.query_sequence
            if not seq or boundary >= len(seq):
                continue

            flank_seq = seq[boundary:]
            if len(flank_seq) < min_flank_length:
                continue

            # Build quality string from BAM quality scores when available.
            quals = record.query_qualities
            if quals is not None and len(quals) > boundary:
                flank_quals = quals[boundary:]
                qual_str = "".join(chr(q + 33) for q in flank_quals)
            else:
                # No quality information stored in BAM; use placeholder.
                qual_str = "I" * len(flank_seq)

            out.write(f"@{read_id}\n{flank_seq}\n+\n{qual_str}\n")
            count += 1

    logger.info("Extracted %d subtelomeric flanks", count)
    return count


def align_subtelomeric(
    fastq_path: str | Path,
    reference: str | Path,
    output_bam: str | Path,
    threads: int = 4,
    platform: Platform = Platform.ONT,
) -> Path:
    """Align subtelomeric flanks to reference using minimap2.

    Args:
        fastq_path: FASTQ of extracted subtelomeric flanks.
        reference: Reference FASTA with subtelomeric sequences.
        output_bam: Output BAM path.
        threads: Number of threads for minimap2.
        platform: Sequencing platform for preset selection.

    Returns:
        Path to sorted, indexed BAM.
    """
    check_external_tools()

    output_bam = Path(output_bam)
    preset = get_minimap2_preset(platform)

    # Build argument lists (no shell, no injection risk).
    minimap2_cmd = [
        "minimap2",
        "-yax", preset,
        "--secondary=no",
        "-t", str(threads),
        str(reference),
        str(fastq_path),
    ]
    samtools_cmd = [
        "samtools", "sort",
        "-@", str(threads),
        "-o", str(output_bam),
    ]

    logger.info(
        "Running alignment: %s | %s",
        " ".join(minimap2_cmd),
        " ".join(samtools_cmd),
    )

    # Pipe minimap2 stdout into samtools stdin without using shell=True.
    minimap2_proc = subprocess.Popen(
        minimap2_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    samtools_proc = subprocess.Popen(
        samtools_cmd,
        stdin=minimap2_proc.stdout,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    # Allow minimap2 to receive SIGPIPE if samtools exits early.
    minimap2_proc.stdout.close()

    samtools_stdout, samtools_stderr = samtools_proc.communicate()
    minimap2_stderr = minimap2_proc.stderr.read()
    minimap2_proc.wait()

    if minimap2_proc.returncode != 0:
        raise RuntimeError(
            f"minimap2 failed (exit {minimap2_proc.returncode}): "
            f"{minimap2_stderr.decode('utf-8', errors='replace')}"
        )
    if samtools_proc.returncode != 0:
        raise RuntimeError(
            f"samtools sort failed (exit {samtools_proc.returncode}): "
            f"{samtools_stderr.decode('utf-8', errors='replace')}"
        )

    # Index the BAM
    try:
        pysam.index(str(output_bam))
    except pysam.SamtoolsError as exc:
        raise RuntimeError(
            f"Failed to index BAM file {output_bam}: {exc}"
        ) from exc

    logger.info("Alignment complete: %s", output_bam)
    return output_bam
