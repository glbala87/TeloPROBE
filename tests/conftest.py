"""Shared test fixtures for teloprobe tests."""

import numpy as np
import pytest

from teloprobe.constants import Platform, Mode
from teloprobe.config import Config
from teloprobe.models.read_record import ReadRecord


@pytest.fixture
def default_config():
    """Default pipeline configuration."""
    config = Config()
    config.min_repeats = 10  # Lower for testing
    config.min_read_length_boundary = 50
    return config


@pytest.fixture
def ont_config(default_config):
    """ONT-specific configuration."""
    default_config.platform = Platform.ONT
    default_config.mode = Mode.TELOSEQ
    default_config.apply_platform_preset()
    default_config.min_repeats = 10
    default_config.min_read_length_boundary = 50
    return default_config


@pytest.fixture
def pacbio_config(default_config):
    """PacBio-specific configuration."""
    default_config.platform = Platform.PACBIO_HIFI
    default_config.mode = Mode.WGS
    default_config.apply_platform_preset()
    default_config.min_repeats = 10
    default_config.min_read_length_boundary = 50
    return default_config


def make_telomeric_sequence(telo_length: int = 1000, subtelo_length: int = 500) -> str:
    """Generate a synthetic telomeric read with known boundary.

    Returns sequence with telomeric repeats (TAACCC) followed by
    random genomic sequence.
    """
    # Telomeric region
    n_repeats = telo_length // 6
    telo = "TAACCC" * n_repeats
    telo = telo[:telo_length]

    # Subtelomeric region (random but not telomere-like)
    rng = np.random.RandomState(42)
    bases = "ACGT"
    subtelo = "".join(rng.choice(list(bases)) for _ in range(subtelo_length))

    return telo + subtelo


def make_tvr_sequence(
    telo_length: int = 600,
    tvr_fraction: float = 0.1,
    subtelo_length: int = 400,
) -> str:
    """Generate sequence with telomere variant repeats."""
    rng = np.random.RandomState(42)
    n_repeats = telo_length // 6
    tvr_motifs = ["TGAGGG", "TCAGGG", "CACCCT"]

    motifs = []
    for _ in range(n_repeats):
        if rng.random() < tvr_fraction:
            motifs.append(rng.choice(tvr_motifs))
        else:
            motifs.append("TAACCC")

    telo = "".join(motifs)[:telo_length]

    bases = "ACGT"
    subtelo = "".join(rng.choice(list(bases)) for _ in range(subtelo_length))

    return telo + subtelo


@pytest.fixture
def telomeric_read():
    """A synthetic read with clear telomeric region."""
    seq = make_telomeric_sequence(1000, 500)
    return ReadRecord(
        read_id="test_read_1",
        sequence=seq,
        qualities=np.full(len(seq), 20, dtype=np.int32),
        read_length=len(seq),
        platform=Platform.ONT,
        mode=Mode.TELOSEQ,
        sample_id="test",
        mean_quality=20.0,
    )


@pytest.fixture
def non_telomeric_read():
    """A synthetic read with no telomeric content."""
    rng = np.random.RandomState(123)
    seq = "".join(rng.choice(list("ACGT")) for _ in range(1500))
    return ReadRecord(
        read_id="test_read_2",
        sequence=seq,
        qualities=np.full(len(seq), 20, dtype=np.int32),
        read_length=len(seq),
        platform=Platform.ONT,
        mode=Mode.TELOSEQ,
        sample_id="test",
        mean_quality=20.0,
    )


@pytest.fixture
def tvr_read():
    """A synthetic read with telomere variant repeats."""
    seq = make_tvr_sequence(600, 0.15, 400)
    return ReadRecord(
        read_id="test_read_tvr",
        sequence=seq,
        qualities=np.full(len(seq), 20, dtype=np.int32),
        read_length=len(seq),
        platform=Platform.ONT,
        mode=Mode.TELOSEQ,
        sample_id="test",
        mean_quality=20.0,
    )
