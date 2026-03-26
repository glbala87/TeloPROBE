"""Integration tests for the full teloprobe pipeline."""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from teloprobe.config import Config
from teloprobe.constants import Platform, Mode, QCStatus
from teloprobe.models.read_record import ReadRecord
from teloprobe.estimate.read_level import process_single_read, compute_read_level
from teloprobe.estimate.confidence import bootstrap_ci
from tests.conftest import make_telomeric_sequence, make_tvr_sequence


class TestProcessSingleRead:
    def test_good_telomeric_read(self, ont_config):
        """A clear telomeric read should be classified as Good."""
        seq = make_telomeric_sequence(1000, 500)
        record = ReadRecord(
            read_id="telo_1", sequence=seq,
            qualities=np.full(len(seq), 20, dtype=np.int32),
            read_length=len(seq), platform=Platform.ONT,
            mode=Mode.TELOSEQ, sample_id="test", mean_quality=20,
        )
        call = process_single_read(record, ont_config)
        assert call.qc_status == QCStatus.GOOD
        assert call.telomere_length is not None
        assert 800 < call.telomere_length < 1200  # approximate

    def test_non_telomeric_read(self, ont_config):
        """A random read should not pass the motif filter."""
        rng = np.random.RandomState(456)
        seq = "".join(rng.choice(list("ACGT")) for _ in range(1500))
        record = ReadRecord(
            read_id="random_1", sequence=seq,
            qualities=np.full(len(seq), 20, dtype=np.int32),
            read_length=len(seq), platform=Platform.ONT,
            mode=Mode.TELOSEQ, sample_id="test", mean_quality=20,
        )
        call = process_single_read(record, ont_config)
        assert call.qc_status != QCStatus.GOOD

    def test_short_read_rejected(self, ont_config):
        ont_config.min_read_length_boundary = 500
        seq = "TAACCC" * 30
        record = ReadRecord(
            read_id="short_1", sequence=seq,
            qualities=np.full(len(seq), 20, dtype=np.int32),
            read_length=len(seq), platform=Platform.ONT,
            mode=Mode.TELOSEQ, sample_id="test", mean_quality=20,
        )
        call = process_single_read(record, ont_config)
        assert call.qc_status == QCStatus.TOO_SHORT

    def test_tvr_read_detected(self, ont_config):
        """TVR composition should be populated."""
        seq = make_tvr_sequence(600, 0.2, 400)
        record = ReadRecord(
            read_id="tvr_1", sequence=seq,
            qualities=np.full(len(seq), 20, dtype=np.int32),
            read_length=len(seq), platform=Platform.ONT,
            mode=Mode.TELOSEQ, sample_id="test", mean_quality=20,
        )
        call = process_single_read(record, ont_config)
        if call.qc_status == QCStatus.GOOD:
            assert call.tvr_fraction > 0
            assert len(call.tvr_composition) > 1  # more than just canonical

    def test_confidence_score_populated(self, ont_config):
        seq = make_telomeric_sequence(800, 400)
        record = ReadRecord(
            read_id="conf_1", sequence=seq,
            qualities=np.full(len(seq), 20, dtype=np.int32),
            read_length=len(seq), platform=Platform.ONT,
            mode=Mode.TELOSEQ, sample_id="test", mean_quality=20,
        )
        call = process_single_read(record, ont_config)
        assert 0 <= call.confidence_score <= 1


class TestComputeReadLevel:
    def test_multiple_reads(self, ont_config):
        """Process a batch of reads."""
        records = []
        rng = np.random.RandomState(42)

        for i in range(5):
            tl = rng.randint(500, 2000)
            seq = make_telomeric_sequence(tl, 500)
            records.append(ReadRecord(
                read_id=f"read_{i}", sequence=seq,
                qualities=np.full(len(seq), 20, dtype=np.int32),
                read_length=len(seq), platform=Platform.ONT,
                mode=Mode.TELOSEQ, sample_id="test", mean_quality=20,
            ))

        df = compute_read_level(records, ont_config)
        assert len(df) == 5
        assert "read_id" in df.columns
        assert "telomere_length_bp" in df.columns
        assert "qc_status" in df.columns
        assert "confidence_score" in df.columns


# ---------------------------------------------------------------------------
# Helper: write synthetic telomeric reads to a FASTQ file
# ---------------------------------------------------------------------------

def _write_synthetic_fastq(path: Path, n_reads: int, telo_length: int = 800,
                           subtelo_length: int = 500, seed: int = 42):
    """Write *n_reads* synthetic telomeric reads to *path* as FASTQ."""
    rng = np.random.RandomState(seed)
    with open(path, "w") as fh:
        for i in range(n_reads):
            tl = max(300, telo_length + rng.randint(-200, 200))
            seq = make_telomeric_sequence(tl, subtelo_length)
            qual = "I" * len(seq)  # Phred ~40
            fh.write(f"@synth_read_{i}\n{seq}\n+\n{qual}\n")


def _write_random_fastq(path: Path, n_reads: int, read_length: int = 1500,
                        seed: int = 99):
    """Write *n_reads* of purely random (non-telomeric) sequence."""
    rng = np.random.RandomState(seed)
    with open(path, "w") as fh:
        for i in range(n_reads):
            seq = "".join(rng.choice(list("ACGT")) for _ in range(read_length))
            qual = "I" * len(seq)
            fh.write(f"@random_read_{i}\n{seq}\n+\n{qual}\n")


# ---------------------------------------------------------------------------
# Pipeline integration tests
# ---------------------------------------------------------------------------

class TestFullPipelineNoMapping:
    """Run ``_run_pipeline`` with skip_mapping=True on synthetic reads."""

    def test_full_pipeline_no_mapping(self, tmp_path):
        """Full pipeline on synthetic telomeric FASTQ produces all outputs."""
        from teloprobe.cli import _run_pipeline

        fastq = tmp_path / "reads.fastq"
        _write_synthetic_fastq(fastq, n_reads=10)

        config = Config()
        config.input_path = str(fastq)
        config.output_dir = str(tmp_path / "out")
        config.sample_id = "synth"
        config.platform = Platform.ONT
        config.mode = Mode.TELOSEQ
        config.skip_mapping = True
        config.min_repeats = 10
        config.min_read_length_boundary = 50
        config.min_length = 100
        config.min_quality = 0
        config.bootstrap_n = 50  # fast
        config.apply_platform_preset()
        # Re-apply relaxed thresholds after preset
        config.min_repeats = 10
        config.min_read_length_boundary = 50
        config.min_length = 100
        config.min_quality = 0

        result = _run_pipeline(config, skip_clustering=True)

        sample_dir = tmp_path / "out" / "synth"
        assert (sample_dir / "read_level.tsv").exists()
        assert (sample_dir / "sample_level.tsv").exists()
        assert (sample_dir / "arm_level.tsv").exists()
        assert (sample_dir / "qc_flags.json").exists()
        assert (tmp_path / "out" / "params.json").exists()
        assert (tmp_path / "out" / "report.html").exists()

        # Verify result is a non-empty DataFrame
        assert result is not None
        assert not result.empty


class TestEmptyInput:
    """Verify graceful handling of an empty input file."""

    def test_empty_input(self, tmp_path):
        """Pipeline should return valid output when FASTQ is empty."""
        from teloprobe.cli import _run_pipeline

        fastq = tmp_path / "empty.fastq"
        fastq.write_text("")  # empty file

        config = Config()
        config.input_path = str(fastq)
        config.output_dir = str(tmp_path / "out")
        config.sample_id = "empty_sample"
        config.platform = Platform.ONT
        config.mode = Mode.TELOSEQ
        config.skip_mapping = True
        config.min_length = 100
        config.min_quality = 0

        result = _run_pipeline(config, skip_clustering=True)

        # Should return a DataFrame (with zero-read summary)
        assert result is not None
        assert isinstance(result, pd.DataFrame)

        # A sample_level.tsv should still be written
        sample_dir = tmp_path / "out" / "empty_sample"
        assert (sample_dir / "sample_level.tsv").exists()

        sample_df = pd.read_csv(sample_dir / "sample_level.tsv", sep="\t")
        assert sample_df.iloc[0]["total_reads"] == 0


class TestSingleRead:
    """Verify pipeline handles a single read without crashing."""

    def test_single_read(self, tmp_path):
        from teloprobe.cli import _run_pipeline

        fastq = tmp_path / "single.fastq"
        _write_synthetic_fastq(fastq, n_reads=1, telo_length=900)

        config = Config()
        config.input_path = str(fastq)
        config.output_dir = str(tmp_path / "out")
        config.sample_id = "single"
        config.platform = Platform.ONT
        config.mode = Mode.TELOSEQ
        config.skip_mapping = True
        config.min_repeats = 10
        config.min_read_length_boundary = 50
        config.min_length = 100
        config.min_quality = 0
        config.bootstrap_n = 50
        config.apply_platform_preset()
        config.min_repeats = 10
        config.min_read_length_boundary = 50
        config.min_length = 100
        config.min_quality = 0

        result = _run_pipeline(config, skip_clustering=True)

        assert result is not None
        assert isinstance(result, pd.DataFrame)

        # read_level.tsv or sample_level.tsv should exist
        sample_dir = tmp_path / "out" / "single"
        assert (sample_dir / "sample_level.tsv").exists()


class TestAllFilteredReads:
    """When all reads fail QC, pipeline should produce valid outputs with zeros."""

    def test_all_filtered_reads(self, tmp_path):
        from teloprobe.cli import _run_pipeline

        fastq = tmp_path / "random.fastq"
        _write_random_fastq(fastq, n_reads=5, read_length=1500)

        config = Config()
        config.input_path = str(fastq)
        config.output_dir = str(tmp_path / "out")
        config.sample_id = "allfiltered"
        config.platform = Platform.ONT
        config.mode = Mode.TELOSEQ
        config.skip_mapping = True
        config.min_repeats = 10
        config.min_read_length_boundary = 50
        config.min_length = 100
        config.min_quality = 0
        config.bootstrap_n = 50
        config.apply_platform_preset()
        config.min_repeats = 10
        config.min_read_length_boundary = 50
        config.min_length = 100
        config.min_quality = 0

        result = _run_pipeline(config, skip_clustering=True)

        assert result is not None
        assert isinstance(result, pd.DataFrame)

        sample_dir = tmp_path / "out" / "allfiltered"
        assert (sample_dir / "sample_level.tsv").exists()

        sample_df = pd.read_csv(sample_dir / "sample_level.tsv", sep="\t")
        assert sample_df.iloc[0].get("informative_reads", 0) == 0


class TestOutputSchema:
    """Verify read_level.tsv has all expected columns."""

    EXPECTED_COLUMNS = [
        "read_id",
        "platform",
        "read_length",
        "telomere_length_bp",
        "telomere_side",
        "motif_score",
        "motif_density",
        "boundary_confidence",
        "boundary_method",
        "subtelomere_anchor",
        "chromosome_arm",
        "haplotype",
        "allele_cluster",
        "tvr_fraction",
        "confidence_score",
        "qc_status",
        "alignment_identity",
        "mapping_quality",
        "anchor_confidence",
    ]

    def test_output_schema(self, ont_config):
        """read_level DataFrame must contain all expected columns."""
        rng = np.random.RandomState(77)
        records = []
        for i in range(3):
            tl = rng.randint(500, 1500)
            seq = make_telomeric_sequence(tl, 500)
            records.append(ReadRecord(
                read_id=f"schema_{i}", sequence=seq,
                qualities=np.full(len(seq), 20, dtype=np.int32),
                read_length=len(seq), platform=Platform.ONT,
                mode=Mode.TELOSEQ, sample_id="test", mean_quality=20,
            ))

        df = compute_read_level(records, ont_config)
        for col in self.EXPECTED_COLUMNS:
            assert col in df.columns, f"Missing expected column: {col}"


class TestBootstrapSmallSample:
    """Verify bootstrap CI works with very few reads (3-5)."""

    def test_bootstrap_3_reads(self):
        """Bootstrap CI should return valid bounds with 3 values."""
        values = np.array([500.0, 700.0, 900.0])
        ci_lower, ci_upper = bootstrap_ci(values, n_boot=200, random_seed=42)
        assert ci_lower <= ci_upper
        assert ci_lower >= 0
        assert np.isfinite(ci_lower)
        assert np.isfinite(ci_upper)

    def test_bootstrap_5_reads(self):
        """Bootstrap CI should return valid bounds with 5 values."""
        values = np.array([400.0, 600.0, 800.0, 1000.0, 1200.0])
        ci_lower, ci_upper = bootstrap_ci(values, n_boot=200, random_seed=42)
        assert ci_lower <= ci_upper
        assert ci_lower > 0
        assert np.isfinite(ci_lower)
        assert np.isfinite(ci_upper)
        # CI should contain or be near the median
        median = np.median(values)
        assert ci_lower <= median <= ci_upper

    def test_bootstrap_2_reads_degenerate(self):
        """With fewer than 3 reads, bootstrap_ci returns point estimate."""
        values = np.array([500.0, 700.0])
        ci_lower, ci_upper = bootstrap_ci(values, n_boot=200, random_seed=42)
        assert ci_lower == ci_upper  # degenerate case
        assert np.isfinite(ci_lower)

    def test_bootstrap_1_read(self):
        """Single value should return that value as both bounds."""
        values = np.array([1000.0])
        ci_lower, ci_upper = bootstrap_ci(values, n_boot=200, random_seed=42)
        assert ci_lower == ci_upper == 1000.0

    def test_bootstrap_empty(self):
        """Empty array should return (0, 0)."""
        values = np.array([])
        ci_lower, ci_upper = bootstrap_ci(values, n_boot=200, random_seed=42)
        assert ci_lower == ci_upper == 0.0


class TestConfigValidation:
    """Test that invalid configs raise ValueError."""

    def test_negative_bootstrap_n(self):
        config = Config()
        config.bootstrap_n = -1
        with pytest.raises(ValueError, match="bootstrap_n"):
            config.validate()

    def test_zero_bootstrap_n(self):
        config = Config()
        config.bootstrap_n = 0
        with pytest.raises(ValueError, match="bootstrap_n"):
            config.validate()

    def test_negative_min_length(self):
        config = Config()
        config.min_length = -10
        with pytest.raises(ValueError, match="min_length"):
            config.validate()

    def test_zero_min_length(self):
        config = Config()
        config.min_length = 0
        with pytest.raises(ValueError, match="min_length"):
            config.validate()

    def test_negative_threads(self):
        config = Config()
        config.threads = 0
        with pytest.raises(ValueError, match="threads"):
            config.validate()

    def test_invalid_bootstrap_alpha_zero(self):
        config = Config()
        config.bootstrap_alpha = 0.0
        with pytest.raises(ValueError, match="bootstrap_alpha"):
            config.validate()

    def test_invalid_bootstrap_alpha_one(self):
        config = Config()
        config.bootstrap_alpha = 1.0
        with pytest.raises(ValueError, match="bootstrap_alpha"):
            config.validate()

    def test_invalid_changepoint_method(self):
        config = Config()
        config.changepoint_method = "invalid_method"
        with pytest.raises(ValueError, match="changepoint_method"):
            config.validate()

    def test_invalid_allele_method(self):
        config = Config()
        config.allele_method = "dbscan"
        with pytest.raises(ValueError, match="allele_method"):
            config.validate()

    def test_negative_min_repeats(self):
        config = Config()
        config.min_repeats = -5
        with pytest.raises(ValueError, match="min_repeats"):
            config.validate()

    def test_valid_config_passes(self):
        """A default Config should pass validation without errors."""
        config = Config()
        config.validate()  # should not raise
