"""Simulate long reads with known telomere lengths for validation.

Generates synthetic FASTQ files with:
- Variable telomere lengths (known ground truth)
- ONT-like or HiFi-like error profiles
- Canonical-only or TVR-rich telomeres
- Interstitial telomeric sequences (negative controls)
- Subtelomeric flanks from reference
"""

import logging
import random
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

CANONICAL_MOTIF = "TAACCC"
VARIANT_MOTIFS = ["TGAGGG", "TCAGGG", "CACCCT", "TTGGGG", "CCCTCA"]


@dataclass
class SimulationParams:
    """Parameters for read simulation."""

    n_reads: int = 1000
    tl_mean: float = 5000.0
    tl_std: float = 2000.0
    tl_min: float = 500.0
    tl_max: float = 20000.0
    subtelo_length: int = 2000
    platform: str = "ont"  # ont or pacbio
    tvr_fraction: float = 0.05
    error_rate: float = 0.0  # 0-1; 0 = no errors
    include_interstitial: int = 50  # number of ITS reads
    include_non_telomeric: int = 100  # random reads
    seed: int = 42


@dataclass
class SimulatedRead:
    """A simulated read with known ground truth."""

    read_id: str
    sequence: str
    quality_string: str
    true_tl: float
    true_arm: str = ""
    read_type: str = "telomeric"  # telomeric, its, random


def simulate_reads(params: SimulationParams) -> list[SimulatedRead]:
    """Generate simulated reads with known telomere lengths."""
    rng = np.random.RandomState(params.seed)
    reads = []

    # Telomeric reads
    for i in range(params.n_reads):
        tl = rng.normal(params.tl_mean, params.tl_std)
        tl = max(params.tl_min, min(params.tl_max, tl))
        tl = int(tl)

        telo_seq = _generate_telomere(tl, params.tvr_fraction, rng)
        subtelo_seq = _generate_random_seq(params.subtelo_length, rng)
        full_seq = telo_seq + subtelo_seq

        if params.error_rate > 0:
            full_seq = _introduce_errors(full_seq, params.error_rate, params.platform, rng)

        qual = _generate_quality(len(full_seq), params.platform, rng)

        reads.append(SimulatedRead(
            read_id=f"sim_telo_{i:06d}",
            sequence=full_seq,
            quality_string=qual,
            true_tl=float(tl),
            read_type="telomeric",
        ))

    # Interstitial telomeric sequence (ITS) reads — negative controls
    for i in range(params.include_interstitial):
        pre_len = rng.randint(500, 3000)
        its_len = rng.randint(100, 800)
        post_len = rng.randint(500, 3000)

        pre = _generate_random_seq(pre_len, rng)
        its = _generate_telomere(its_len, 0.0, rng)
        post = _generate_random_seq(post_len, rng)
        full_seq = pre + its + post

        if params.error_rate > 0:
            full_seq = _introduce_errors(full_seq, params.error_rate, params.platform, rng)

        qual = _generate_quality(len(full_seq), params.platform, rng)

        reads.append(SimulatedRead(
            read_id=f"sim_its_{i:06d}",
            sequence=full_seq,
            quality_string=qual,
            true_tl=0.0,
            read_type="its",
        ))

    # Non-telomeric reads — negative controls
    for i in range(params.include_non_telomeric):
        seq_len = rng.randint(1000, 10000)
        full_seq = _generate_random_seq(seq_len, rng)

        if params.error_rate > 0:
            full_seq = _introduce_errors(full_seq, params.error_rate, params.platform, rng)

        qual = _generate_quality(len(full_seq), params.platform, rng)

        reads.append(SimulatedRead(
            read_id=f"sim_random_{i:06d}",
            sequence=full_seq,
            quality_string=qual,
            true_tl=0.0,
            read_type="random",
        ))

    rng.shuffle(reads)
    return reads


def write_simulated_fastq(
    reads: list[SimulatedRead],
    output_path: str | Path,
) -> Path:
    """Write simulated reads to FASTQ file."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        for read in reads:
            f.write(f"@{read.read_id}\n")
            f.write(f"{read.sequence}\n")
            f.write("+\n")
            f.write(f"{read.quality_string}\n")

    logger.info("Wrote %d simulated reads to %s", len(reads), output_path)
    return output_path


def write_ground_truth(
    reads: list[SimulatedRead],
    output_path: str | Path,
) -> Path:
    """Write ground truth telomere lengths to TSV."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        f.write("read_id\ttrue_tl\tread_type\ttrue_arm\n")
        for read in reads:
            f.write(f"{read.read_id}\t{read.true_tl}\t{read.read_type}\t{read.true_arm}\n")

    return output_path


def _generate_telomere(length: int, tvr_fraction: float, rng: np.random.RandomState) -> str:
    """Generate telomeric sequence with optional TVRs."""
    n_repeats = length // 6 + 1
    motifs = []
    for _ in range(n_repeats):
        if rng.random() < tvr_fraction:
            motifs.append(rng.choice(VARIANT_MOTIFS))
        else:
            motifs.append(CANONICAL_MOTIF)
    return "".join(motifs)[:length]


def _generate_random_seq(length: int, rng: np.random.RandomState) -> str:
    """Generate random non-telomeric sequence."""
    # Avoid generating telomere-like sequence by using biased base frequencies
    bases = list("ACGT")
    return "".join(rng.choice(bases, size=length))


def _introduce_errors(
    seq: str,
    error_rate: float,
    platform: str,
    rng: np.random.RandomState,
) -> str:
    """Introduce sequencing errors into a sequence."""
    seq_list = list(seq)
    bases = list("ACGT")

    for i in range(len(seq_list)):
        if rng.random() < error_rate:
            error_type = rng.random()

            if platform == "ont":
                # ONT: higher deletion rate in homopolymers
                if error_type < 0.4:
                    # Substitution
                    seq_list[i] = rng.choice([b for b in bases if b != seq_list[i]])
                elif error_type < 0.7:
                    # Deletion (mark for removal)
                    seq_list[i] = ""
                else:
                    # Insertion
                    seq_list[i] = seq_list[i] + rng.choice(bases)
            else:
                # PacBio HiFi: mostly substitutions
                if error_type < 0.8:
                    seq_list[i] = rng.choice([b for b in bases if b != seq_list[i]])
                elif error_type < 0.9:
                    seq_list[i] = ""
                else:
                    seq_list[i] = seq_list[i] + rng.choice(bases)

    return "".join(seq_list)


def _generate_quality(length: int, platform: str, rng: np.random.RandomState) -> str:
    """Generate quality string appropriate for platform."""
    if platform == "ont":
        # ONT: mean ~15, range 5-30
        quals = rng.normal(15, 3, length).clip(5, 35).astype(int)
    else:
        # PacBio HiFi: mean ~30, range 20-40
        quals = rng.normal(30, 3, length).clip(20, 40).astype(int)

    return "".join(chr(q + 33) for q in quals)
