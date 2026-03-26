"""Tests for read filtering module."""

import numpy as np
import pytest

from teloprobe.candidate.read_filter import (
    apply_preflight_filters,
    check_error_motifs,
    check_start_repeats,
    check_post_boundary_quality,
    check_post_boundary_ccc,
    _largest_cluster,
)
from teloprobe.constants import QCStatus
from teloprobe.config import Config
from teloprobe.models.read_record import ReadRecord


class TestPreflightFilters:
    def test_too_short(self):
        config = Config()
        config.min_read_length_boundary = 160
        record = ReadRecord(
            read_id="r1", sequence="A" * 100,
            qualities=None, read_length=100, mean_quality=20,
        )
        assert apply_preflight_filters(record, config) == QCStatus.TOO_SHORT

    def test_passes(self):
        config = Config()
        config.min_read_length_boundary = 160
        record = ReadRecord(
            read_id="r1", sequence="A" * 500,
            qualities=None, read_length=500, mean_quality=20,
        )
        assert apply_preflight_filters(record, config) == QCStatus.GOOD


class TestErrorMotifs:
    def test_no_errors(self):
        seq = "TAACCC" * 200
        assert check_error_motifs(seq, 600, max_errors=5, error_distance=500) is False

    def test_clustered_errors(self):
        # Insert error motifs in a cluster
        seq = "TAACCC" * 50 + ("GTATAG" * 10) + "TAACCC" * 50
        assert check_error_motifs(seq, 300, max_errors=5, error_distance=500) is True


class TestStartRepeats:
    def test_good_start(self):
        motif = np.concatenate([np.ones(500), np.zeros(500)])
        assert check_start_repeats(motif, 500, 0.3, 0.8) == False

    def test_bad_start(self):
        motif = np.concatenate([np.zeros(200), np.ones(300), np.zeros(500)])
        assert check_start_repeats(motif, 500, 0.3, 0.8) == True


class TestPostBoundaryQuality:
    def test_good_quality(self):
        quals = np.full(1000, 20, dtype=np.int32)
        assert check_post_boundary_quality(quals, 500, 9) is False

    def test_low_quality(self):
        quals = np.full(1000, 5, dtype=np.int32)
        assert check_post_boundary_quality(quals, 500, 9) is True

    def test_no_quality(self):
        assert check_post_boundary_quality(None, 500, 9) is False


class TestPostBoundaryCCC:
    def test_normal_subtelomere(self):
        seq = "TAACCC" * 100 + "ATCGATCGATCG" * 50
        assert check_post_boundary_ccc(seq, 600, 0.25) is False

    def test_all_telomere(self):
        seq = "TAACCC" * 200
        assert check_post_boundary_ccc(seq, 600, 0.25) is True


class TestLargestCluster:
    def test_cluster(self):
        positions = [100, 110, 120, 500, 510]
        assert _largest_cluster(positions, 50) == 3

    def test_no_cluster(self):
        positions = [100, 300, 500, 700]
        assert _largest_cluster(positions, 50) == 1

    def test_empty(self):
        assert _largest_cluster([], 50) == 0
